"""The interactive page in a real browser (Chromium, driven by Playwright): what the stub harness of
tests/test_na_layout.py cannot see, such as where a click lands. Validation only: Playwright and Chromium are not
dependencies of foldmap. The tests skip unless node, the Playwright package (FOLDMAP_PLAYWRIGHT, a directory whose
node_modules holds it; ~/foldmap-remote/pw by default) and a Chromium build (FOLDMAP_CHROMIUM, else the newest
under ~/.cache/ms-playwright) are present."""

import json
import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest

from foldmap.interactive import HIT_PX, build_page
from foldmap.style import Style

pytestmark = pytest.mark.browser

DATA = Path(__file__).parent / "data"


def _playwright() -> Path | None:
    root = Path(os.environ.get("FOLDMAP_PLAYWRIGHT", Path.home() / "foldmap-remote" / "pw"))
    module = root / "node_modules" / "playwright"
    return module if (module / "package.json").is_file() else None


def _chromium() -> Path | None:
    if os.environ.get("FOLDMAP_CHROMIUM"):
        path = Path(os.environ["FOLDMAP_CHROMIUM"])
        return path if path.is_file() else None
    builds = Path.home().glob(".cache/ms-playwright/chromium-*/chrome-linux*/chrome")
    found = sorted(builds, key=lambda p: int(re.search(r"chromium-(\d+)", str(p)).group(1)))
    return found[-1] if found else None


node, playwright, chromium = shutil.which("node"), _playwright(), _chromium()
needs_browser = pytest.mark.skipif(
    node is None or playwright is None or chromium is None,
    reason="needs node, Playwright (FOLDMAP_PLAYWRIGHT) and Chromium (FOLDMAP_CHROMIUM); none is downloaded here",
)

# Opens the page with the 3D viewer stubbed (the CDN is blocked: no network) and clicks in the figure. Each check
# finds its own point from the drawing's geometry, PX screen pixels off a line's centre, clicks there and reads the
# selection; a control click at the same point with the click target switched off shows what the target changed.
_DRIVER = r"""
const [pagePath, chrome, playwright, px, wanted] = process.argv.slice(2);
const { chromium } = require(playwright);
const PX = Number(px);

const HELPERS = () => {
  window.__t = {
    data: JSON.parse(document.getElementById('topo-data').textContent),
    ink(id) {  // the path drawn under a group id
      const g = document.getElementById(id);
      return g && (g.tagName === 'path' ? g : g.querySelector('path'));
    },
    screen(el, s) {  // a point along a path, in client pixels
      const p = el.getPointAtLength(s), m = el.getScreenCTM();
      return { x: m.a * p.x + m.c * p.y + m.e, y: m.b * p.x + m.d * p.y + m.f };
    },
    offsets(el, px, n) {  // points px client pixels either side of a path's centre line, at n places along it
      const L = el.getTotalLength(), out = [];
      for (let i = 1; i < n; i++) {
        const s = L * i / n, a = this.screen(el, Math.max(0, s - 0.5)), b = this.screen(el, Math.min(L, s + 0.5));
        const c = this.screen(el, s), dx = b.x - a.x, dy = b.y - a.y, h = Math.hypot(dx, dy);
        if (h < 0.2) continue;
        for (const sign of [1, -1]) out.push({ x: c.x - sign * dy / h * px, y: c.y + sign * dx / h * px });
      }
      return out;
    },
    around(el, px, n) {  // points px client pixels from a path (no part of it nearer), in 16 directions round n places
      const L = el.getTotalLength(), on = [], out = [];
      for (let i = 0; i <= n; i++) on.push(this.screen(el, L * i / n));
      on.forEach(c => { for (let k = 0; k < 16; k++) {
        const p = { x: c.x + px * Math.cos(k * Math.PI / 8), y: c.y + px * Math.sin(k * Math.PI / 8) };
        if (on.every(q => Math.hypot(q.x - p.x, q.y - p.y) > px - 0.15)) out.push(p);
      } });
      return out;
    },
    widthPx(el) {  // a path's drawn stroke width in client pixels
      return parseFloat(getComputedStyle(el).strokeWidth) * Math.hypot(el.getScreenCTM().a, el.getScreenCTM().b);
    },
    listed(hitId) {  // the selection key a nucleic-acid target stands for, when the page lists that run or stem
      const m = /^na:(run|stem):(\d+):hit$/.exec(hitId), key = m && m[1] + ':' + m[2];
      return key && this.data.links.nucleics.some(x => x.key === key) ? 'nucleic:' + key : null;
    },
  };
};

(async () => {
  const browser = await chromium.launch({ executablePath: chrome });
  const page = await browser.newPage({ viewport: { width: 1400, height: 900 } });
  await page.route(/^https?:/, r => r.abort());
  await page.addInitScript(() => {  // the 3D panel's library, stubbed: every call a no-op
    globalThis.$3Dmol = { createViewer: () => new Proxy({}, { get: () => () => {} }) };
  });
  const errors = [];
  page.on('pageerror', e => errors.push(String(e)));
  await page.goto("file://" + pagePath);
  await page.waitForFunction(() => typeof window.foldmapSelection === 'function');
  await page.evaluate(HELPERS);

  async function clickAt(p) {
    await page.mouse.click(p.x, p.y);
    await page.waitForTimeout(20);
    const sel = await page.evaluate(() => window.foldmapSelection());
    await page.keyboard.press('Escape');
    return sel;
  }
  async function control(hitId, p) {  // the same click with that target switched off
    await page.evaluate(id => { document.getElementById(id).style.pointerEvents = 'none'; }, hitId);
    const sel = await clickAt(p);
    await page.evaluate(id => { document.getElementById(id).style.pointerEvents = ''; }, hitId);
    return sel;
  }
  // a click PX off a run's or stem's line, where its target is the topmost thing under the pointer (strict) or,
  // for a backbone drawn wider than the target, its own drawn line; at the fitted view, or zoomed in on that line
  // (four wheel steps, 1.75x: the target keeps its width on screen)
  async function offLine(zoomLabel, strict) {
    const hits = await page.evaluate(() => [...document.querySelectorAll('path[id^="na:"][id$=":hit"]')]
                                             .map(h => h.id).filter(id => __t.listed(id)));
    const res = [];
    for (const id of hits) {
      await page.click('#zfit');
      if (zoomLabel === 'in') {
        const m = await page.evaluate(id => { const h = document.getElementById(id);
                                              return __t.screen(h, h.getTotalLength() / 2); }, id);
        await page.mouse.move(m.x, m.y);
        for (let k = 0; k < 4; k++) await page.mouse.wheel(0, -100);
        await page.waitForTimeout(50);
      }
      const f = await page.evaluate(([id, px, strict]) => {
        const hit = document.getElementById(id), base = id.slice(0, -4), run = base.startsWith('na:run:');
        const ink = run ? __t.ink(base) : __t.ink(base + ':rails');
        const own = n => {  // the run's or stem's own drawn line
          const g = n && n.closest('g[id]');
          return !!g && (g.id === base || g.id.startsWith(base + ':'));
        };
        const pts = [];
        for (const p of __t.offsets(run ? ink : hit, px, 80)) {
          const top = document.elementFromPoint(p.x, p.y);
          if (top === hit || (!strict && own(top))) pts.push({ ...p, on: top === hit ? 'target' : 'ink' });
        }
        if (!pts.length) return null;
        return { want: __t.listed(id), p: pts[Math.floor(pts.length / 2)], ink: ink && __t.widthPx(ink) };
      }, [id, PX, strict]);
      if (!f) continue;
      res.push({ zoom: zoomLabel, hit: id, ...f, selection: await clickAt(f.p), without: await control(id, f.p) });
    }
    await page.click('#zfit');
    return res;
  }

  // a click PX off a protein loop's line where a nucleic-acid target also lies: the loop wins
  async function offLoopNearNA(zoomLabel) {
    const f = await page.evaluate(px => {
      const na = [...document.querySelectorAll('path[id^="na:"][id$=":hit"]')].filter(n => __t.listed(n.id));
      for (const hit of document.querySelectorAll('path[id^="loop:"][id$=":hit"]')) {
        const ink = __t.ink(hit.id.slice(0, -4));
        for (const p of __t.around(ink, px, Math.ceil(ink.getTotalLength()))) {
          if (document.elementFromPoint(p.x, p.y) !== hit) continue;
          const under = document.elementsFromPoint(p.x, p.y).find(n => na.includes(n));
          if (!under) continue;
          return { hit: hit.id, want: hit.id.slice(0, -4), na: __t.listed(under.id), p, ink: __t.widthPx(ink) };
        }
      }
      return null;
    }, PX);
    if (!f) return null;
    // just off the same loop's line, over its white casing (and nothing else)
    const close = await page.evaluate(id => {
      const hit = document.getElementById(id), ink = __t.ink(id.slice(0, -4));
      const casing = __t.ink(id.slice(0, -4).replace('loop:', 'loopcase:'));
      const px = (__t.widthPx(ink) + __t.widthPx(casing)) / 4;  // halfway between the line's edge and the casing's
      const pts = __t.around(ink, px, Math.ceil(ink.getTotalLength())).filter(p => {
        const local = new DOMPoint(p.x, p.y).matrixTransform(casing.getScreenCTM().inverse());
        return document.elementFromPoint(p.x, p.y) === hit && casing.isPointInStroke(local);
      });
      return pts.length ? { ...pts[Math.floor(pts.length / 2)], px } : null;
    }, f.hit);
    return { zoom: zoomLabel, ...f, selection: await clickAt(f.p), without: await control(f.hit, f.p),
             close: close && { p: close, selection: await clickAt(close) } };
  }

  // a hover and a click PX off an interaction line (ia:<kind>:...), where its target is the topmost thing there
  async function offInteraction() {
    await page.click('#zfit');
    const found = await page.evaluate(px => {
      const out = [];
      for (const hit of document.querySelectorAll('path[id^="ia:"][id$=":hit"]')) {
        const x = __t.data.interactions.find(x => x.key === hit.id.slice(0, -4));
        const pts = x ? __t.offsets(hit, px, 40).filter(p => document.elementFromPoint(p.x, p.y) === hit) : [];
        if (pts.length) out.push({ hit: hit.id, want: x.key, note: x.text, p: pts[Math.floor(pts.length / 2)] });
      }
      return out.slice(0, 3);
    }, PX);
    const res = [];
    for (const f of found) {
      await page.mouse.move(f.p.x, f.p.y);
      await page.waitForTimeout(20);
      const info = await page.evaluate(() => document.getElementById('info').textContent);
      res.push({ ...f, info, selection: await clickAt(f.p), without: await control(f.hit, f.p) });
    }
    return res;
  }

  const out = { errors, lines: [], loop: null, interactions: [] };
  if (wanted.includes('interactions')) out.interactions = await offInteraction();
  if (wanted.includes('lines')) {
    const strict = !wanted.includes('ink');
    out.lines.push(...await offLine('fit', strict));
    out.lines.push(...await offLine('in', strict));
  }
  if (wanted.includes('loop')) out.loop = await offLoopNearNA('fit');
  console.log(JSON.stringify(out));
  await browser.close();
})().catch(e => { console.error(e); process.exit(1); });
"""


def _drive(tmp_path, page: str, wanted: str) -> dict:
    (tmp_path / "page.html").write_text(page, encoding="utf-8")
    (tmp_path / "driver.js").write_text(_DRIVER, encoding="utf-8")
    args = [node, str(tmp_path / "driver.js"), str(tmp_path / "page.html"), str(chromium), str(playwright)]
    res = subprocess.run([*args, "3", wanted], capture_output=True, text=True, timeout=300)
    assert res.returncode == 0, res.stderr
    out = json.loads(res.stdout)
    assert out["errors"] == []
    return out


def _check_lines(out):
    lines = out["lines"]
    kinds = {(r["zoom"], r["hit"].split(":")[1]) for r in lines}
    assert {("fit", "run"), ("in", "run"), ("fit", "stem"), ("in", "stem")} <= kinds, lines
    for r in lines:
        assert r["selection"] == [r["want"]], r  # 3 px off its line: the run or stem is selected ...
        if r["p"]["on"] == "target":  # ... off its drawn line, where without the target the click selects nothing
            assert r["without"] == [], r
        else:  # ... on its own drawn line, 6 px wide or more there, which is a target in its own right
            assert r["ink"] > 2 * 3 - 0.5 and r["without"] == [r["want"]], r
    return lines


@needs_browser
def test_a_click_three_pixels_off_a_run_or_stem_selects_it(tmp_path):
    # Q-N3 / 13b-5: on the page a run or stem is selected within a few pixels of its line, at any zoom
    assert HIT_PX >= 2 * 3 + 2
    lines = _check_lines(_drive(tmp_path, build_page(DATA / "1LMB.cif", viewer="3dmol"), "lines"))
    assert all(r["p"]["on"] == "target" for r in lines)


@needs_browser
def test_a_thick_backbone_keeps_its_click_target(tmp_path):
    # na_line (46a5291) thickens the drawn backbone: the target is still at least HIT_PX wide on screen, and the
    # thicker line is its own target where it is wider
    page = build_page(DATA / "1LMB.cif", viewer="3dmol", look=Style(na_line=2.0))
    lines = _check_lines(_drive(tmp_path, page, "lines ink"))
    assert any(r["p"]["on"] == "ink" for r in lines if r["zoom"] == "in")  # zoomed in, the run is wider than 8 px


@needs_browser
def test_a_click_three_pixels_off_a_loop_beside_a_strand_still_selects_the_loop(tmp_path):
    # the strand's wider target must not take the clicks of a protein loop running beside it: 3 px off the loop,
    # inside the strand's target, the loop is selected; with the loop's own target off, the strand would have been
    # 1IHF: a loop of chain A passes within about 1 px of a single strand at the fitted view
    out = _drive(tmp_path, build_page(DATA / "na" / "1IHF.cif.gz", viewer="3dmol"), "loop")
    loop = out["loop"]
    assert loop is not None, "no protein loop runs within reach of a strand's target"
    assert loop["selection"] == [loop["want"]], loop
    assert loop["ink"] < 2 * 3, loop  # 3 px off the loop's centre is off its drawn line
    assert loop["without"] == [loop["na"]], loop
    # just off its line, over its white casing, the loop is selected too (the casing lets clicks through to the target)
    assert loop["close"] is not None and loop["close"]["selection"] == [loop["want"]], loop


@needs_browser
def test_a_hover_or_click_three_pixels_off_an_interaction_line_shows_and_selects_it(tmp_path):
    # plan 13c-5: an interaction line is a thin dash; its click target makes it a few pixels wide on screen, so a
    # hover shows its note (residues, atoms, distances) and a click selects it, where without the target nothing is
    page = build_page(DATA / "1LMB.cif", viewer="3dmol", look=Style(interactions="salt,cation-pi"))
    found = _drive(tmp_path, page, "interactions")["interactions"]
    assert len(found) == 3, found
    for r in found:
        assert r["info"] == r["note"] and r["selection"] == [r["want"]], r
        assert r["want"] not in r["without"], r
