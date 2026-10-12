"""Render jobs for the app: the same figure as the command line, limits that hold, failures remembered."""

import json
import os
import signal
import sys
import time
from pathlib import Path

import pytest

from foldmap import app_jobs
from foldmap.app_jobs import JobLimit, Limits, run_job

DATA = Path(__file__).parent / "data"
needs_fork = pytest.mark.skipif(not app_jobs._CAN_FORK, reason="isolated jobs need forkserver (not on Windows)")
linux = pytest.mark.skipif(sys.platform != "linux", reason="RLIMIT_AS is reliable on Linux only")
ENV = (
    "FOLDMAP_APP_TIMEOUT",
    "FOLDMAP_APP_MEMORY_MB",
    "FOLDMAP_APP_MAX_RESIDUES",
    "FOLDMAP_APP_JOBS",
    "FOLDMAP_APP_SLOW",
)


def _blob(**settings):
    return json.dumps({"theme": "publication", **settings}, sort_keys=True, default=str)


def _sleep(seconds):
    time.sleep(seconds)
    return seconds


def _hog():
    return len(bytearray(4 * 1024**3))


def _boom():
    raise KeyError("x")


def _gone():
    raise FileNotFoundError("the upload is gone")


def _headroom():
    """(MB the job may still take under its address-space cap, its BLAS thread settings)."""
    import resource

    size = int(Path("/proc/self/statm").read_text(encoding="utf-8").split()[0]) * resource.getpagesize()
    soft, _ = resource.getrlimit(resource.RLIMIT_AS)
    return (soft - size) / 2**20, os.environ.get("OPENBLAS_NUM_THREADS"), os.environ.get("OMP_NUM_THREADS")


def _exits():
    os._exit(3)  # dies without a word, but not killed for its memory


def _large():
    return "x" * 5000


def _refuse_cap(memory_mb):
    return False  # as on macOS, where the address-space cap is refused


def _killed():
    os.kill(os.getpid(), signal.SIGKILL)  # what the OS does to a process over its memory


def test_app_figure_is_the_command_line_figure():
    from foldmap.render import save_svg
    from foldmap.spec import draw_spec, figure_spec

    for name, extra, given in (("1LMB", {}, {}), ("2LZM", {"swap": "α1,α3"}, {"swap": ["α1,α3"]})):
        out = app_jobs.draw_job(DATA / f"{name}.cif", _blob(**extra))
        cli = save_svg(draw_spec(figure_spec(DATA / f"{name}.cif", theme="publication", **given)).figure)
        assert out["svg"] == cli, name


@needs_fork
def test_isolated_job_draws_the_same_figure():
    """A figure drawn in a job process is the figure drawn in this one, byte for byte (protein; nucleic acid)."""
    for name in ("2LZM", "1LMB"):
        here = app_jobs.draw_job(DATA / f"{name}.cif", _blob())
        there = run_job(app_jobs.draw_job, DATA / f"{name}.cif", _blob(), limits=Limits(timeout=120))
        assert there["svg"] == here["svg"], name
        assert there["elements"] == here["elements"] and there["summary"] == here["summary"], name


def test_draw_job_reports_what_the_pickers_need():
    out = app_jobs.draw_job(DATA / "1UBQ.cif", _blob())
    assert out["summary"] == "1 chain · 76 residues · 1 helix · 5 strands"
    assert out["chains"] == ["A"] and out["na_chains"] == []
    assert [(e["ref"], e["label"]) for e in out["elements"][:4]] == [
        ("res:A:2-7", "A"), ("res:A:12-16", "B"), ("res:A:23-34", "α1"), ("res:A:38-40", "η1"),
    ]  # fmt: skip
    assert out["elements"][0] == {
        "ref": "res:A:2-7",
        "label": "A",
        "chain": "A",
        "kind": "E",
        "first": "2",
        "last": "7",
    }
    assert out["fallback"] == [] and out["over_labels"] == [] and out["files"] == {}


def test_summary_never_reads_like_a_fault():
    out = app_jobs.draw_job(DATA / "1LMB.cif", _blob())
    assert out["summary"] == "2 chains · 179 residues · 10 helices · 40 nucleotides · symmetry C2"
    assert out["chains"] == ["3", "4"] and out["na_chains"] == ["1", "2"]


def test_slow_figures_bring_their_files():
    out = app_jobs.draw_job(DATA / "1UBQ.cif", _blob(), 0.0)
    assert out["files"]["png"][:4] == b"\x89PNG" and out["files"]["pdf"][:4] == b"%PDF"


def test_files_reuse_the_drawn_figure(monkeypatch):
    import foldmap.spec as spec

    app_jobs.clear()
    app_jobs.draw_job(DATA / "1UBQ.cif", _blob())
    calls = []
    real = spec.draw_spec
    monkeypatch.setattr(spec, "draw_spec", lambda *a, **k: (calls.append(1), real(*a, **k))[1])
    png = app_jobs.files_job(DATA / "1UBQ.cif", _blob(), "png")
    assert png[:4] == b"\x89PNG" and calls == []
    assert app_jobs.files_job(DATA / "1UBQ.cif", _blob(theme="trace"), "svg").startswith(b"<svg") and calls == [1]


def test_contacts_are_checked_on_request():
    report = app_jobs.contacts_job(DATA / "1LMB.cif", _blob())
    assert report["drawn"] == 21 and report["verified"] == 21 and report["false"] == 0
    assert report["problems"] == []
    assert app_jobs.contacts_job(DATA / "1UBQ.cif", _blob()) is None


def test_size_job_counts_before_drawing():
    size = app_jobs.size_job(DATA / "1UBQ.cif")
    assert size["residues"] == 76 and size["chains"] == 1 and isinstance(size["assemblies"], list)


def test_limits_come_from_the_environment(monkeypatch):
    for name in ENV:
        monkeypatch.delenv(name, raising=False)
    assert app_jobs.limits(hosted=True) == Limits(timeout=90, memory_mb=2048, max_residues=30_000, jobs=2)
    assert app_jobs.limits(hosted=False) == Limits()
    monkeypatch.setenv("FOLDMAP_APP_TIMEOUT", "30")
    monkeypatch.setenv("FOLDMAP_APP_MAX_RESIDUES", "nonsense")
    monkeypatch.setenv("FOLDMAP_APP_SLOW", "0")
    got = app_jobs.limits(hosted=True)
    assert got.timeout == 30.0 and got.max_residues == 30_000 and got.slow == 0.0


def test_slow_threshold_is_read_when_the_limits_are(monkeypatch):
    for name in ENV:
        monkeypatch.delenv(name, raising=False)
    monkeypatch.setattr(app_jobs, "SLOW_SECONDS", 0.0)
    assert app_jobs.limits(hosted=False).slow == 0.0 and app_jobs.limits(hosted=True).slow == 0.0


def test_limits_fall_back_without_fork(monkeypatch):
    monkeypatch.setattr(app_jobs, "_CAN_FORK", False)
    assert not Limits(timeout=5).isolated
    assert run_job(_sleep, 0.01, limits=Limits(timeout=5)) == 0.01  # in-process, no hard limit


def test_local_jobs_run_in_process():
    assert not Limits().isolated
    with pytest.raises(KeyError):  # nothing wrapped: the app's catch-all sees the real error
        run_job(_boom, limits=Limits())


@needs_fork
def test_timeout_kills_the_job():
    t0 = time.perf_counter()
    with pytest.raises(JobLimit) as err:
        run_job(_sleep, 10, limits=Limits(timeout=1.0))
    assert err.value.kind == "time" and time.perf_counter() - t0 < 5


@needs_fork
@linux
def test_memory_cap_stops_the_job():
    with pytest.raises(JobLimit) as err:
        run_job(_hog, limits=Limits(timeout=30, memory_mb=2048))
    assert err.value.kind == "memory"


@needs_fork
@linux
def test_a_job_has_its_whole_memory_allowance():
    """The cap is the job's own size plus memory_mb, and the job server runs BLAS on one thread (per-thread buffers
    of a many-core machine would otherwise eat the allowance); this process keeps its own settings."""
    before = (os.environ.get("OPENBLAS_NUM_THREADS"), os.environ.get("OMP_NUM_THREADS"))
    room, openblas, omp = run_job(_headroom, limits=Limits(timeout=60, memory_mb=1024))
    assert room >= 0.9 * 1024 and (openblas, omp) == ("1", "1")
    assert (os.environ.get("OPENBLAS_NUM_THREADS"), os.environ.get("OMP_NUM_THREADS")) == before


@needs_fork
def test_a_job_killed_from_outside_is_a_memory_limit():
    with pytest.raises(JobLimit) as err:
        run_job(_killed, limits=Limits(timeout=60))
    assert err.value.kind == "memory"


@needs_fork
def test_errors_cross_the_process_boundary():
    with pytest.raises(ValueError, match="needs an alignment"):
        run_job(app_jobs.draw_job, DATA / "1UBQ.cif", _blob(color_by="conservation"), limits=Limits(timeout=60))
    with pytest.raises(app_jobs.JobFailed, match="KeyError"):
        run_job(_boom, limits=Limits(timeout=60))
    with pytest.raises(OSError, match="the upload is gone"):
        run_job(_gone, limits=Limits(timeout=60))


@needs_fork
def test_a_busy_app_says_so(monkeypatch):
    """Every job slot taken (held here, so nothing depends on timing): the next job waits BUSY_WAIT, then gives up."""
    monkeypatch.setattr(app_jobs, "BUSY_WAIT", 0.2)
    slots = app_jobs._semaphore(1)
    assert slots.acquire(timeout=60)
    try:
        with pytest.raises(JobLimit) as err:
            run_job(_sleep, 0, limits=Limits(timeout=60, jobs=1))
        assert err.value.kind == "busy"
    finally:
        slots.release()
    assert run_job(_sleep, 0, limits=Limits(timeout=60, jobs=1)) == 0  # the slot is free again


def test_failures_are_remembered_by_layout_not_by_colour():
    app_jobs.clear()
    err = JobLimit("time", "slow")
    app_jobs.remember_failure("/x/8TOC.cif", {"theme": "publication", "palette": "okabe-ito"}, err)
    assert app_jobs.known_failure("/x/8TOC.cif", {"theme": "trace", "palette": "tol-bright"}) is err
    assert app_jobs.known_failure("/x/8TOC.cif", {"chains": "A"}) is None


def test_a_busy_moment_is_never_remembered():
    app_jobs.clear()
    app_jobs.remember_failure("/x/1UBQ.cif", {}, JobLimit("busy", "busy"))
    assert app_jobs.known_failure("/x/1UBQ.cif", {}) is None and not app_jobs.FAILED


def test_sequence_and_3d_failures_are_kept_apart_from_the_figure():
    app_jobs.clear()
    slow = JobLimit("time", "slow")
    app_jobs.remember_failure("/x/8TOC.cif", {}, slow, job="page")
    assert app_jobs.known_failure("/x/8TOC.cif", {}) is None
    assert app_jobs.known_failure("/x/8TOC.cif", {}, job="sequence") is None
    assert app_jobs.known_failure("/x/8TOC.cif", {}, job="page") is slow


def test_failure_memory_keeps_the_latest():
    app_jobs.clear()
    for k in range(300):
        app_jobs.remember_failure(f"/x/{k}.cif", {}, JobLimit("time", "slow"))
    assert len(app_jobs.FAILED) == 256
    assert app_jobs.known_failure("/x/0.cif", {}) is None and app_jobs.known_failure("/x/299.cif", {}) is not None


def test_limit_messages_say_what_to_do():
    lim = app_jobs.HOSTED
    text = app_jobs.size_message("8TOC", {"residues": 28299, "chains": 181}, lim)
    assert text.startswith("8TOC is too large for this app: 28,299 residues in 181 chain(s)")
    assert "foldmap plot 8TOC -o 8TOC.svg" in text
    assert app_jobs.limit_message("8C3A", JobLimit("time", ""), lim).startswith(
        "8C3A took longer than 90 s to draw, the limit of this app."
    )
    assert app_jobs.limit_message("8C3A", JobLimit("memory", ""), lim) == (
        "8C3A needs more memory than this app allows (2048 MB). Draw part of it (Chains & assembly), or on your "
        "own computer: pip install foldmap, then foldmap plot 8C3A -o 8C3A.svg."
    )
    assert app_jobs.limit_message("8C3A", JobLimit("busy", ""), lim) == (
        "The app is busy drawing other figures; try again in a minute."
    )


def test_one_handler_tells_every_failure_apart():
    from foldmap.app import trouble

    lim = app_jobs.HOSTED
    explain = app_jobs.explain
    assert explain(JobLimit("time", ""), "8C3A", lim) == (
        app_jobs.limit_message("8C3A", JobLimit("time", ""), lim),
        "Chains (Chains & assembly)",
    )
    big = JobLimit("size", "8TOC is too large for this app")
    assert explain(big, "8TOC", lim) == ("8TOC is too large for this app", None)
    assert explain(ValueError("--move needs ref=dx,dy, got 'α2'"), "2LZM", lim) == (
        "--move needs ref=dx,dy, got 'α2'",
        "Move (Edit layout)",
    )
    crash = "foldmap could not draw this setting (IndexError: boom). Please report it at https://github.com/dzyla/foldmap/issues."
    # a crash in a job process (JobFailed is a RuntimeError, told apart first) and one in this process read alike
    assert explain(app_jobs.JobFailed("IndexError: boom"), "2LZM", lim, "domains") == (crash, "Domains (Edit layout)")
    assert explain(IndexError("boom"), "2LZM", lim, "domains") == (crash, "Domains (Edit layout)")
    gone = RuntimeError("gemmi cannot read /srv/uploads/x.cif")
    assert explain(gone, "x", lim, "domains") == (trouble(gone), None)  # no server path, no setting to blame
    assert explain(FileNotFoundError("/srv/uploads/x.cif"), "x", lim) == (trouble(OSError()), None)


def test_download_files_come_from_the_render_when_it_brought_them(monkeypatch):
    app_jobs.clear()
    lim = Limits()
    out = app_jobs.draw_job(DATA / "1UBQ.cif", _blob(), 0.0)
    eager = tuple(sorted(out["files"].items()))
    monkeypatch.setattr(app_jobs, "run_job", lambda *a, **k: pytest.fail("drew again"))
    assert app_jobs.figure_file(str(DATA / "1UBQ.cif"), _blob(), "png", lim, eager) == out["files"]["png"]


def test_download_callables_draw_on_demand():
    app_jobs.clear()
    lim, path = Limits(), str(DATA / "1UBQ.cif")
    assert app_jobs.figure_file(path, _blob(), "pdf", lim)[:4] == b"%PDF"
    assert app_jobs.figure_file(path, _blob(), "svg", lim, ()).startswith(b"<svg")  # no eager files: drawn
    view = json.dumps({"columns": 30, "ss_colour": "figure"}, sort_keys=True)
    assert app_jobs.sequence_file(path, _blob(), view, "svg", lim).startswith(b"<svg")
    assert app_jobs.sequence_file(path, _blob(), view, "png", lim)[:4] == b"\x89PNG"
    page = app_jobs.page_file(path, _blob(), lim)
    assert isinstance(page, bytes) and b"<html" in page[:200].lower()
    assert app_jobs.figure_file.cache_info().currsize == 2
    app_jobs.clear()
    assert app_jobs.figure_file.cache_info().currsize == 0 and not app_jobs._DRAWN


def test_sequence_job_is_the_sequence_download():
    from foldmap.app import sequence_files

    view = {"columns": 30, "ss_colour": "figure"}
    old = sequence_files(DATA / "1UBQ.cif", {"theme": "publication"}, view)
    got = app_jobs.sequence_job(DATA / "1UBQ.cif", _blob(), json.dumps(view), "svg")
    assert got == old["svg"]


def test_contact_problems_read_as_the_command_line_prints_them():
    from types import SimpleNamespace as NS

    from foldmap.cli import contact_problems

    report = NS(
        beads=[
            NS(problem="", nucleotide="1:DA5", residue="3:ARG10", pair=""),
            NS(problem="false", nucleotide="1:DT6", residue="3:LYS12", pair="3:LYS12-1:DT6"),
            NS(problem="unmatched", nucleotide="2:DG7", residue="4:SER3", pair=""),
        ],
        missing=[NS(nucleotide="2:DC8", pair="4:ASN5-2:DC8")],
    )
    assert contact_problems(report) == [
        "false: 1:DT6 <- 3:LYS12 (3:LYS12-1:DT6)",
        "unmatched: 2:DG7 <- 4:SER3",
        "missing: 2:DC8 <- 4:ASN5-2:DC8",
    ]


def test_a_memory_cap_always_comes_with_a_timeout(monkeypatch):
    for name in ENV:
        monkeypatch.delenv(name, raising=False)
    assert Limits(memory_mb=512).timeout == app_jobs.ISOLATED_TIMEOUT == 90
    monkeypatch.setenv("FOLDMAP_APP_TIMEOUT", "0")
    assert app_jobs.limits(hosted=True).timeout == 90
    assert Limits(timeout=0).timeout == 0  # no isolation asked for: runs here, untimed


@needs_fork
def test_a_job_dying_otherwise_is_a_failure_not_a_memory_limit():
    with pytest.raises(app_jobs.JobFailed, match="exit code 3"):
        run_job(_exits, limits=Limits(timeout=60))


@needs_fork
def test_a_result_too_large_to_hand_back_is_a_failure(monkeypatch):
    monkeypatch.setattr(app_jobs, "RESULT_MAX", 1000)
    with pytest.raises(app_jobs.JobFailed, match="larger than"):
        run_job(_large, limits=Limits(timeout=60))
    monkeypatch.setattr(app_jobs, "RESULT_MAX", 10_000)
    assert run_job(_large, limits=Limits(timeout=60)) == "x" * 5000


def test_no_isolation_is_said_once(monkeypatch, caplog):
    monkeypatch.setattr(app_jobs, "_CAN_FORK", False)
    monkeypatch.setattr(app_jobs, "_WARNED", set())
    with caplog.at_level("WARNING", logger="foldmap.app_jobs"):
        for _ in range(3):
            assert run_job(_sleep, 0, limits=Limits(timeout=5)) == 0
    assert len(caplog.records) == 1 and "no time or memory limit" in caplog.text


@needs_fork
def test_no_memory_cap_is_said_once(monkeypatch, caplog):
    monkeypatch.setattr(app_jobs, "_set_cap", _refuse_cap)
    monkeypatch.setattr(app_jobs, "_WARNED", set())
    with caplog.at_level("WARNING", logger="foldmap.app_jobs"):
        for _ in range(3):
            assert run_job(_sleep, 0, limits=Limits(timeout=60, memory_mb=512)) == 0
    assert len(caplog.records) == 1 and "no memory cap" in caplog.text


def test_a_sequence_failure_is_kept_with_its_view():
    app_jobs.clear()
    slow = JobLimit("time", "slow")
    app_jobs.remember_failure("/x/8TOC.cif", {}, slow, job="sequence", view_json='{"columns": 60}')
    assert app_jobs.known_failure("/x/8TOC.cif", {}, job="sequence", view_json='{"columns": 60}') is slow
    assert app_jobs.known_failure("/x/8TOC.cif", {}, job="sequence", view_json='{"columns": 30}') is None


def test_limits_point_at_chains_and_assembly():
    from foldmap.app_edits import control_for

    assert control_for("", "chains") == "Chains (Chains & assembly)"
    lim = app_jobs.HOSTED
    for kind in ("time", "memory"):
        assert app_jobs.explain(JobLimit(kind, ""), "8C3A", lim)[1] == "Chains (Chains & assembly)"
    assert app_jobs.explain(JobLimit("busy", ""), "8C3A", lim)[1] is None


def test_contacts_key_names_only_what_changes_the_contacts():
    from foldmap.app_jobs import contacts_key

    base = {"chains": "", "na_extent": "auto", "palette": "okabe-ito", "theme": "publication"}
    look = {"palette": "tol-bright", "theme": "trace", "background": "#000000", "rotate": 30.0}
    assert contacts_key(base) == contacts_key({**base, **look})
    assert contacts_key(base) != contacts_key({**base, "chains": "A"})
    assert contacts_key(base) != contacts_key({**base, "na_extent": "all"})
