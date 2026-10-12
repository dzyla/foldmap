# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Render jobs for the app: every expensive call (figure, files, sequence, 3D page, contacts check, size check) is a
module-level job function, run through run_job.

On your own machine a job runs in this process. Hosted, it runs in a process forked from a preloaded forkserver, with
a hard timeout, an address-space cap and a cap on jobs at once; going over a limit raises JobLimit, and the failure
is remembered (remember_failure) so the next visitor asking for the same layout is answered at once.

Errors: a job raises ValueError (the input is wrong), OSError / RuntimeError (the file is gone or unreadable),
JobLimit (over a limit) or, from a job process, JobFailed (anything else: a bug). JobFailed is a RuntimeError, so a
caller that tells them apart itself catches JobFailed before RuntimeError; explain() does it for every caller and
returns the message to show and the control to point at.

No Streamlit here: the app wraps these in its own caches."""

from __future__ import annotations

import functools
import json
import logging
import multiprocessing
import os
import re
import signal
import sys
import tempfile
import threading
import time
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path

SLOW_SECONDS = 4.0  # a figure slower than this turns automatic redrawing off and brings its files with it (Q-A6)
PAGE_EMBED_MAX = 10_000_000  # a 3D page larger than this is offered as a download only, never embedded
BUSY_WAIT = 30.0  # seconds a job waits for a free slot before the app says it is busy
ISOLATED_TIMEOUT = 90.0  # the timeout of a job given a memory cap and no timeout: an isolated job is always timed
RESULT_MAX = 100_000_000  # bytes a job process may hand back (a giant's 3D page is 20-48 MB)

_CAN_FORK = sys.platform != "win32" and "forkserver" in multiprocessing.get_all_start_methods()
_PRELOAD = [
    "foldmap.app",
    "foldmap.app_jobs",
    "foldmap.spec",
    "foldmap.cli",
    "foldmap.render",
    "foldmap.interactive",
    "foldmap.seqplot",
]
REPORT = "https://github.com/dzyla/foldmap/issues"
REPORT_SENTENCE = f" Please report it at {REPORT}."  # ends every bug message (the app links it)
_LOG = logging.getLogger("foldmap.app_jobs")
_WARNED: set[str] = set()  # warnings already given in this process


def _warn_once(key: str, message: str) -> None:
    if key not in _WARNED:
        _WARNED.add(key)
        _LOG.warning(message)


@dataclass(frozen=True)
class Limits:
    """How a job runs: timeout (s) and memory_mb (address space) > 0 run it isolated, in its own process; jobs: at
    most this many at once (0: no cap); max_residues: larger structures are refused before drawing (0: none); slow:
    seconds after which a figure counts as slow. A memory cap without a timeout gets ISOLATED_TIMEOUT: a job in
    its own process is never waited on for ever."""

    timeout: float = 0.0
    memory_mb: int = 0
    max_residues: int = 0
    jobs: int = 0
    slow: float = SLOW_SECONDS

    def __post_init__(self):
        if self.memory_mb > 0 and self.timeout <= 0:
            object.__setattr__(self, "timeout", ISOLATED_TIMEOUT)

    @property
    def isolated(self) -> bool:
        """Whether jobs run in their own process: a limit is set and forkserver exists (not on Windows)."""
        return (self.timeout > 0 or self.memory_mb > 0) and _CAN_FORK


HOSTED = Limits(timeout=90, memory_mb=2048, max_residues=30_000, jobs=2)  # Q-A11


def limits(hosted: bool) -> Limits:
    """HOSTED or no limits, each overridden by its FOLDMAP_APP_* variable (an unparsable value keeps the default)."""
    base = HOSTED if hosted else Limits()

    def env(name: str, cast, default):
        try:
            return cast(os.environ[name])
        except (KeyError, ValueError):
            return default

    return Limits(
        timeout=env("FOLDMAP_APP_TIMEOUT", float, base.timeout),
        memory_mb=env("FOLDMAP_APP_MEMORY_MB", int, base.memory_mb),
        max_residues=env("FOLDMAP_APP_MAX_RESIDUES", int, base.max_residues),
        jobs=env("FOLDMAP_APP_JOBS", int, base.jobs),
        slow=env("FOLDMAP_APP_SLOW", float, SLOW_SECONDS),
    )


class JobLimit(Exception):
    """A job over a limit. kind: 'time', 'memory', 'busy' (no free slot) or 'size' (refused before drawing)."""

    def __init__(self, kind: str, message: str = ""):
        super().__init__(kind, message)
        self.kind, self.message = kind, message

    def __str__(self) -> str:
        return self.message or self.kind


class JobFailed(RuntimeError):
    """An unexpected exception inside an isolated job, as '<Type>: <text>'. A RuntimeError: catch it first."""


# ---------------------------------------------------------------------------------------------- running a job
_CONTEXT = None
_LOCK = threading.Lock()
_SLOTS: dict[int, threading.BoundedSemaphore] = {}


_ONE_THREAD = {"OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1"}


def _context():
    """The forkserver context, made once, with its server running (started again if it died). The server imports
    the pipeline once for every job, with BLAS on one thread: OpenBLAS sets aside buffers per thread, which on a
    many-core machine take much of a job's memory allowance before it starts. This process keeps its own settings."""
    global _CONTEXT
    from multiprocessing import forkserver

    with _LOCK:
        if _CONTEXT is None:
            _CONTEXT = multiprocessing.get_context("forkserver")
            _CONTEXT.set_forkserver_preload(_PRELOAD)
        saved = {k: os.environ.get(k) for k in _ONE_THREAD}
        os.environ.update(_ONE_THREAD)  # the server inherits this environment when it starts
        try:
            forkserver.ensure_running()
        finally:
            for k, v in saved.items():
                if v is None:
                    os.environ.pop(k, None)
                else:
                    os.environ[k] = v
        return _CONTEXT


def _semaphore(jobs: int) -> threading.BoundedSemaphore:
    """The process-wide slots for `jobs` jobs at once."""
    with _LOCK:
        if jobs not in _SLOTS:
            _SLOTS[jobs] = threading.BoundedSemaphore(jobs)
        return _SLOTS[jobs]


def _own_size() -> int:
    """This process's address space in bytes (Linux /proc; 0 where there is none)."""
    try:
        import resource

        return int(Path("/proc/self/statm").read_text(encoding="utf-8").split()[0]) * resource.getpagesize()
    except (ImportError, OSError, ValueError, IndexError):
        return 0


def _set_cap(memory_mb: int) -> bool:
    """Cap this process's address space at its present size plus memory_mb; False where the OS refuses (macOS)."""
    try:
        import resource

        cap = _own_size() + memory_mb * 1024 * 1024  # memory_mb on top of what the job holds from the server
        resource.setrlimit(resource.RLIMIT_AS, (cap, cap))
    except (ImportError, ValueError, OSError):
        return False
    return True


def _child(conn, fn, args, memory_mb: int, max_bytes: int, set_cap=_set_cap) -> None:
    """In the job process: cap the address space, run, send (kind, payload, capped): ('ok', result), ('memory',
    ''), ('large', size in bytes) when the result is over max_bytes, or ('error', (type name, text, is ValueError,
    is OSError))."""
    from multiprocessing.reduction import ForkingPickler

    capped = set_cap(memory_mb) if memory_mb else False
    try:
        message = ForkingPickler.dumps(("ok", fn(*args), capped))
        if len(message) > max_bytes:
            message = ForkingPickler.dumps(("large", len(message), capped))
    except MemoryError:
        message = ForkingPickler.dumps(("memory", "", capped))
    except Exception as err:  # every kind of error goes back to the app, which says what it means
        what = (type(err).__name__, str(err), isinstance(err, ValueError), isinstance(err, OSError))
        message = ForkingPickler.dumps(("error", what, capped))
    try:
        conn.send_bytes(message)
    finally:
        conn.close()


def run_job(fn, *args, limits: Limits):
    """fn(*args): here when `limits` sets no isolation, else in a job process under the time and memory caps,
    one of at most `limits.jobs` at once. Over a limit: JobLimit. A ValueError, OSError or RuntimeError in the job
    is raised again as such (with its text); anything else, a job process ending without an answer (other than
    killed) or a result over RESULT_MAX, as JobFailed."""
    if not limits.isolated:
        if limits.timeout > 0 or limits.memory_mb > 0:
            _warn_once("fork", "jobs run in the app's own process, with no time or memory limit (no forkserver here)")
        return fn(*args)
    slots = _semaphore(limits.jobs) if limits.jobs > 0 else None
    if slots is not None and not slots.acquire(timeout=BUSY_WAIT):
        raise JobLimit("busy", "The app is busy drawing other figures; try again in a minute.")
    try:
        return _isolated(fn, args, limits)
    finally:
        if slots is not None:
            slots.release()


def _isolated(fn, args, limits: Limits):
    ctx = _context()
    here, there = ctx.Pipe(duplex=False)
    proc = ctx.Process(target=_child, args=(there, fn, args, limits.memory_mb, RESULT_MAX, _set_cap), daemon=True)
    try:
        proc.start()
    except BaseException:
        here.close()
        there.close()
        raise
    there.close()  # the job holds the only sending end: its death reads as end of file here
    answered = False
    try:
        if not here.poll(limits.timeout or None):
            raise JobLimit("time", f"stopped after {limits.timeout:g} s")
        try:
            kind, payload, capped = here.recv()
            answered = True
        except EOFError:  # ended without a word
            proc.join(5)
            if proc.exitcode == -signal.SIGKILL:  # killed: by the OS, which takes it for its memory
                raise JobLimit("memory", "the job process was killed") from None
            raise JobFailed(f"the job process ended without a result (exit code {proc.exitcode})") from None
    finally:
        if answered:
            proc.join(5)  # it has answered and is closing
        if proc.is_alive():
            proc.kill()
            proc.join(5)
        here.close()
    if limits.memory_mb and not capped:
        _warn_once("cap", "jobs run with no memory cap: this system refuses an address-space limit (RLIMIT_AS)")
    if kind == "ok":
        return payload
    if kind == "memory":
        raise JobLimit("memory", f"stopped at {limits.memory_mb} MB")
    if kind == "large":
        raise JobFailed(
            f"the result is {payload / 1e6:.0f} MB, larger than the {RESULT_MAX / 1e6:.0f} MB a job may hand back"
        )
    name, text, is_value, is_os = payload
    if is_value:
        raise ValueError(text)
    if is_os:
        raise OSError(text)
    if name == "RuntimeError":
        raise RuntimeError(text)
    raise JobFailed(f"{name}: {text}")


# ---------------------------------------------------------------------------------------------- jobs
_DRAWN: OrderedDict = OrderedDict()  # (path, settings JSON) -> spec.Drawn, the last 4 drawn in this process
_DRAWN_MAX = 4
_FIGURES = threading.Lock()  # a drawn figure is saved by one thread at a time (downloads run in threads)


def _keep_drawn(key: tuple[str, str], drawn) -> None:
    with _LOCK:
        _DRAWN[key] = drawn
        _DRAWN.move_to_end(key)
        while len(_DRAWN) > _DRAWN_MAX:
            _DRAWN.popitem(last=False)


def _drawn(path, settings_json: str):
    """The Drawn for these settings: the one kept from an earlier job in this process, else drawn now (and kept)."""
    from . import app
    from . import spec as spec_mod

    key = (str(path), settings_json)
    with _LOCK:
        drawn = _DRAWN.get(key)
        if drawn is not None:
            _DRAWN.move_to_end(key)
    if drawn is None:
        drawn = spec_mod.draw_spec(app.spec_of(path, json.loads(settings_json)))
        _keep_drawn(key, drawn)
    return drawn


def figure_bytes(fig, ext: str) -> bytes:
    """A figure as SVG (inline, as the app shows it), PNG (300 dpi) or PDF: the files the app offers."""
    from . import render

    if ext == "svg":
        return render.save_svg(fig).encode()
    with tempfile.TemporaryDirectory() as tmp:
        return render.save(fig, Path(tmp) / f"figure.{ext}").read_bytes()


def draw_job(path, settings_json: str, eager_after: float | None = None) -> dict:
    """Draw the figure: its SVG and what the app's pickers and notes need. When it took longer than
    `eager_after` seconds, its PNG and PDF come with it (a job process cannot hand back the figure itself)."""
    from . import app, render
    from . import spec as spec_mod

    t0 = time.perf_counter()
    drawn = spec_mod.draw_spec(app.spec_of(path, json.loads(settings_json)))
    svg = render.save_svg(drawn.figure)
    seconds = time.perf_counter() - t0
    files = {}
    if eager_after is not None and seconds > eager_after:
        files = {ext: figure_bytes(drawn.figure, ext) for ext in ("png", "pdf")}
    chains = list(dict.fromkeys(label.chain for label in drawn.bb.labels))
    nucleic = getattr(drawn.layout, "nucleic", None)
    na_chains = list(dict.fromkeys(s.chain for s in nucleic.strands)) if nucleic is not None else []
    out = {
        "svg": svg,
        "fallback": drawn.stuck,
        "over_labels": drawn.over_labels,
        "elements": elements_of(drawn),
        "chains": chains,
        "na_chains": na_chains,
        "summary": summary_of(drawn),
        "seconds": seconds,
        "files": files,
    }
    _keep_drawn((str(path), settings_json), drawn)  # after every save: another thread may save it from here on
    return out


def files_job(path, settings_json: str, ext: str) -> bytes:
    """The figure as an svg, png or pdf file, from the figure kept in this process when there is one."""
    drawn = _drawn(path, settings_json)
    with _FIGURES:
        return figure_bytes(drawn.figure, ext)


def sequence_job(path, settings_json: str, view_json: str, ext: str) -> bytes:
    """The sequence figure (app.sequence_figure) as an svg, png or pdf file."""
    from . import app

    return figure_bytes(app.sequence_figure(path, json.loads(settings_json), json.loads(view_json)), ext)


def page_job(path, settings_json: str) -> str:
    """The interactive page (figure beside the 3D structure) as HTML."""
    from . import app
    from . import interactive as interactive_mod

    return interactive_mod.build_page(path, **app.spec_of(path, json.loads(settings_json)).page_kwargs())


def contacts_job(path, settings_json: str) -> dict | None:
    """The contact beads of the figure checked against the coordinates (foldmap contacts): counts and one line per
    problem; None when the figure has no beads."""
    from .cli import contact_problems

    drawn = _drawn(path, settings_json)
    with _FIGURES:  # the kept figure may be saved by a download thread meanwhile
        report = drawn.verify_contacts()
    if report is None:
        return None
    return {
        "drawn": len(report.beads),
        "verified": len(report.verified),
        "false": len(report.false),
        "missing": len(report.missing),
        "unmatched": len(report.unmatched),
        "wrong_colour": len(report.wrong_colour),
        "problems": contact_problems(report),
    }


# The settings that decide the contacts check's answer: what is drawn (chains, assembly, symmetry and the subunit
# shown, which nucleic acid, the elements), whether the beads are drawn (letters or beads) and which residues have a
# colour (the colouring and the data behind it). The bead colours themselves are checked against the figure's own
# colouring, so a palette, a theme's colours or a position change gives the same answer. The app's cache of the check
# is keyed on these settings only and does not hash the settings JSON it passes to contacts_job (`_blob`): that is
# safe only while this list is complete, so a setting that may change which beads are drawn belongs here (a superset
# costs a check more, never a wrong answer).
CONTACT_KEYS = (
    "chains",
    "mode",
    "assembly",
    "symmetry",
    "symmetry_tol",
    "protomers",
    "focus",
    "highlight",
    "dna",
    "na_extent",
    "na_style",
    "nucleotide_labels",
    "helices_310",
    "color_by",
    "msa_path",
    "msa_reference",
    "domains",
    "uniprot",
    "uniprot_domains",
)


def contacts_key(settings: dict) -> str:
    """The settings of CONTACT_KEYS as JSON: the key of the app's contacts check cache."""
    return json.dumps({k: settings[k] for k in CONTACT_KEYS if k in settings}, sort_keys=True, default=str)


def size_job(path) -> dict:
    """What the file holds, read before any drawing: polymer residues and chains of its first model (as
    deposited; `chain_ids` in the file's order, nucleic acid too), and its assemblies."""
    import gemmi

    st = gemmi.read_structure(str(path))
    st.setup_entities()
    polymers = [(chain.name, len(chain.get_polymer())) for chain in st[0]] if len(st) else []
    polymers = [(name, n) for name, n in polymers if n]
    return {
        "residues": sum(n for _, n in polymers),
        "chains": len(polymers),
        "chain_ids": list(dict.fromkeys(name for name, _ in polymers)),
        "assemblies": [a.name for a in st.assemblies],
    }


def elements_of(drawn) -> list[dict]:
    """The drawn elements (no ghosts) in sequence order, as the app's pickers show them."""
    from .model import res_ref

    placed = drawn.layout.placed
    return [
        {
            "ref": res_ref(s.chain, s.first, s.last),
            "label": placed[s.id].label,
            "chain": s.chain,
            "kind": s.kind,
            "first": s.first.tag,
            "last": s.last.tag,
        }
        for s in drawn.sses
        if s.id in placed and not placed[s.id].ghost
    ]


def _count(n: int, one: str, many: str) -> str:
    return f"{n:,} {one if n == 1 else many}"


def summary_of(drawn) -> str:
    """One line about the drawn structure: chains, residues, helices, strands, nucleotides, symmetry; zero counts
    left out."""
    from .symmetry import detect_symmetry

    bb, sses = drawn.bb, drawn.sses
    try:
        sym = detect_symmetry(bb, sses)
    except ValueError:
        sym = None
    nucleic = getattr(drawn.layout, "nucleic", None)
    counts = (
        (len({label.chain for label in bb.labels}), "chain", "chains"),
        (len(bb), "residue", "residues"),
        (sum(s.kind == "H" for s in sses), "helix", "helices"),
        (sum(s.kind == "E" for s in sses), "strand", "strands"),
        (nucleic.size if nucleic is not None else 0, "nucleotide", "nucleotides"),
    )
    parts = [_count(n, one, many) for n, one, many in counts if n]
    if sym is not None:
        parts.append(f"symmetry {sym.label}")
    return " · ".join(parts)


def restricted_summary(summary: str, total: int) -> str:
    """The summary of a figure drawn with only some of the structure's `total` protein chains: '2 of 6 chains · …'
    (unchanged when it names no chain count, or all of them)."""
    m = re.match(r"([\d,]+) chains?\b", summary)
    if not m or int(m.group(1).replace(",", "")) >= total:
        return summary
    return f"{m.group(1)} of {total:,} chains{summary[m.end() :]}"


# ---------------------------------------------------------------------------------------------- downloads
@functools.lru_cache(maxsize=24)
def figure_file(path: str, blob: str, ext: str, limits: Limits, eager: tuple | None = None) -> bytes:
    """A figure download: from `eager` (the render's files as sorted (ext, bytes) items) when it holds `ext`,
    else drawn by a job."""
    files = dict(eager or ())
    if ext in files:
        return files[ext]
    return run_job(files_job, path, blob, ext, limits=limits)


@functools.lru_cache(maxsize=16)
def sequence_file(path: str, blob: str, view_json: str, ext: str, limits: Limits) -> bytes:
    """A sequence figure download."""
    return run_job(sequence_job, path, blob, view_json, ext, limits=limits)


@functools.lru_cache(maxsize=16)
def page_file(path: str, blob: str, limits: Limits) -> bytes:
    """The interactive page download (UTF-8 HTML)."""
    return run_job(page_job, path, blob, limits=limits).encode("utf-8")


# ---------------------------------------------------------------------------------------------- failures
FAILED: OrderedDict = OrderedDict()  # (path, job, layout_key, view_json) -> JobLimit, process-wide, the latest 256
_FAILED_MAX = 256


def layout_key(settings: dict) -> str:
    """The settings that decide the layout (app.LAYOUT_DEFAULTS and style.LAYOUT_KEYS), as JSON: a figure over a
    limit stays over it in any colour."""
    from . import app
    from .style import LAYOUT_KEYS

    keys = [*app.LAYOUT_DEFAULTS, *LAYOUT_KEYS]
    return json.dumps({k: settings[k] for k in keys if k in settings}, sort_keys=True, default=str)


def _failure_key(path, settings: dict, job: str, view_json: str) -> tuple:
    return (str(path), job, layout_key(settings), view_json)


def remember_failure(path, settings: dict, err: JobLimit, job: str = "figure", view_json: str = "") -> None:
    """Remember a job over a limit, so the same request is answered at once. job: 'figure', 'sequence' or 'page'
    (each remembered apart); view_json: the sequence view it was drawn with. A busy moment says nothing about the
    structure and is never remembered."""
    if err.kind == "busy":
        return
    key = _failure_key(path, settings, job, view_json)
    with _LOCK:
        FAILED[key] = err
        FAILED.move_to_end(key)
        while len(FAILED) > _FAILED_MAX:
            FAILED.popitem(last=False)


def known_failure(path, settings: dict, job: str = "figure", view_json: str = "") -> JobLimit | None:
    """The limit this job went over before for this layout (and sequence view), else None."""
    with _LOCK:
        return FAILED.get(_failure_key(path, settings, job, view_json))


def _elsewhere(name: str) -> str:
    return f"pip install foldmap, then foldmap plot {name} -o {name}.svg"


def size_message(name: str, size: dict, limits: Limits) -> str:
    return (
        f"{name} is too large for this app: {size['residues']:,} residues in {size['chains']} chain(s), and the "
        f"limit here is {limits.max_residues:,}. Draw it on your own computer: {_elsewhere(name)} (add --chains "
        "A,B to draw part of it)."
    )


def limit_message(name: str, err: JobLimit, limits: Limits) -> str:
    """What to tell the visitor about a job over a limit."""
    ending = f"Draw part of it (Chains & assembly), or on your own computer: {_elsewhere(name)}."
    if err.kind == "time":
        return f"{name} took longer than {limits.timeout:g} s to draw, the limit of this app. {ending}"
    if err.kind == "memory":
        return f"{name} needs more memory than this app allows ({limits.memory_mb} MB). {ending}"
    if err.kind == "busy":
        return "The app is busy drawing other figures; try again in a minute."
    return err.message


_CHAINS = "Chains (Chains & assembly)"  # app_edits.CONTROLS["chains"], as control_for names it


def explain(err: Exception, name: str, limits: Limits, last_changed: str | None = None) -> tuple[str, str | None]:
    """The one handler for a failed job: (message for the visitor, control to point at, as 'Label (Section)', or
    None). JobLimit: limit_message, and for time and memory the Chains control. ValueError: its text and the control
    it names (else the one changed last).
    JobFailed (checked before RuntimeError) or any exception that is not ValueError, OSError or RuntimeError (a
    bug, in a job process or in this one): 'foldmap could not draw this setting', with a link to report it.
    OSError, RuntimeError: app.trouble (the file is gone or unreadable; no server path, no setting to blame)."""
    from . import app
    from .app_edits import control_for

    if isinstance(err, JobLimit):  # time and memory: draw part of it, the message says
        return limit_message(name, err, limits), _CHAINS if err.kind in ("time", "memory") else None
    if isinstance(err, ValueError):
        return str(err), control_for(str(err), last_changed)
    if isinstance(err, JobFailed) or not isinstance(err, (OSError, RuntimeError)):
        text = str(err) if isinstance(err, JobFailed) else f"{type(err).__name__}: {err}"
        message = f"foldmap could not draw this setting ({text}).{REPORT_SENTENCE}"
        return message, control_for(text, last_changed)
    return app.trouble(err), None


def clear() -> None:
    """Forget every kept figure, download and failure (tests)."""
    for cached in (figure_file, sequence_file, page_file):
        cached.cache_clear()
    with _LOCK:
        _DRAWN.clear()
        FAILED.clear()
