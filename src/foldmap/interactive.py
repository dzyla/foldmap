"""One self-contained HTML page linking three views of a structure: the topology figure, the CA distance
(contact) map and the 3D model. Hovering an element, a map cell or an atom highlights the same residues in
the other two views."""

from __future__ import annotations

import html
import io
import json
from pathlib import Path

import gemmi
import numpy as np

from .io import read_model
from .render import draw, element_colours, save_svg
from .route import route_loops
from .style import Style

THREEDMOL = "https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.1.0/3Dmol-min.js"
_LOOP_GREY = "#9aa3ad"

# DOM-free helpers: unit-tested under node (tests/test_interactive.py)
_LIB = r"""
(function (g) {
  const TopoLib = {
    binning(n, max) {
      const size = Math.max(1, Math.ceil(n / max));
      return { size, count: Math.ceil(n / size) };
    },
    distance(ca, i, j) {
      const a = ca[i], b = ca[j];
      return Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
    },
    matrix(ca, max) {  // closest CA-CA distance within each pair of residue bins
      const { size, count } = this.binning(ca.length, max);
      const values = new Float32Array(count * count);
      for (let I = 0; I < count; I++) {
        for (let J = I; J < count; J++) {
          let m = Infinity;
          for (let i = I * size; i < Math.min(ca.length, (I + 1) * size); i++) {
            for (let j = J * size; j < Math.min(ca.length, (J + 1) * size); j++) {
              const d = this.distance(ca, i, j);
              if (d < m) m = d;
            }
          }
          values[I * count + J] = m;
          values[J * count + I] = m;
        }
      }
      return { size, count, values };
    },
    index(data) {
      const byKey = {};
      data.residues.forEach((r, k) => { byKey[r.c + ':' + r.n + (r.i || '')] = k; });
      const byElement = {};
      data.elements.forEach(e => { byElement[e.id] = e; });
      return { byKey, byElement };
    },
    elementAt(data, k) { return data.residues[k] ? data.residues[k].e : null; },
    span(data, id) {
      const e = data.elements.find(x => x.id === id);
      return e ? [e.start, e.end] : null;
    },
    loopSpan(data, a, b) {  // residues strictly between two elements
      const ea = data.elements.find(x => x.id === a), eb = data.elements.find(x => x.id === b);
      return ea && eb ? [ea.end + 1, eb.start - 1] : null;
    },
    ramp(t) {  // close (dark blue) to far (pale): distances drawn on the contact map
      const stops = [[16, 42, 92], [33, 102, 172], [67, 162, 202], [168, 221, 181], [247, 252, 240]];
      const x = Math.min(1, Math.max(0, t)) * (stops.length - 1), k = Math.min(stops.length - 2, Math.floor(x)), f = x - k;
      return stops[k].map((v, i) => Math.round(v + f * (stops[k + 1][i] - v)));
    },
  };
  g.TopoLib = TopoLib;
})(globalThis);
"""

_APP = r"""
(function () {
  const data = JSON.parse(document.getElementById('topo-data').textContent);
  const lib = globalThis.TopoLib, idx = lib.index(data);
  const info = document.getElementById('info');
  const svg = document.querySelector('#topology svg');
  const MAX = 700, CUT = 8;

  // ---- topology: element groups by id
  const parts = {};  // element id -> [svg nodes]
  svg.querySelectorAll('[id]').forEach(node => {
    const m = node.id.match(/^(strand|helix|helix-back|eta|label|ghost|ghostlabel|resnum):(.+?)(:[NC])?$/);
    if (m) (parts[m[2]] = parts[m[2]] || []).push(node);
  });
  svg.querySelectorAll('[id^="loop:"], [id^="loop-arrow:"]').forEach(node => {
    const [a, b] = node.id.split(':').slice(1).join(':').split('>');
    const key = 'loop ' + a + '>' + b;
    (parts[key] = parts[key] || []).push(node);
  });
  function svgOn(keys) {
    svg.classList.toggle('focus', keys.length > 0);
    svg.querySelectorAll('.on').forEach(n => n.classList.remove('on'));
    keys.forEach(k => (parts[k] || []).forEach(n => n.classList.add('on')));
  }

  // ---- contact map
  const base = document.getElementById('map'), over = document.getElementById('map-over');
  const M = lib.matrix(data.ca, MAX);
  [base, over].forEach(c => { c.width = M.count; c.height = M.count; });
  const ctx = base.getContext('2d'), img = ctx.createImageData(M.count, M.count);
  for (let p = 0; p < M.values.length; p++) {
    const [r, g, b] = lib.ramp(M.values[p] / 30);
    img.data.set([r, g, b, 255], p * 4);
  }
  ctx.putImageData(img, 0, 0);
  const strip = document.getElementById('strip');
  strip.width = M.count; strip.height = 1;
  const sctx = strip.getContext('2d');
  data.elements.forEach(e => {  // the sequence strip under the map: each element in its figure colour
    sctx.fillStyle = e.colour;
    sctx.fillRect(Math.floor(e.start / M.size), 0, Math.max(1, Math.ceil((e.end - e.start + 1) / M.size)), 1);
  });
  const octx = over.getContext('2d');
  function mapOn(ranges, pair) {
    octx.clearRect(0, 0, M.count, M.count);
    octx.fillStyle = 'rgba(255, 106, 0, 0.28)';
    ranges.forEach(([a, b]) => {
      const s = Math.floor(a / M.size), w = Math.max(1, Math.floor(b / M.size) - s + 1);
      octx.fillRect(s, 0, w, M.count);
      octx.fillRect(0, s, M.count, w);
    });
    if (pair) {
      octx.fillStyle = 'rgba(255, 106, 0, 0.95)';
      const [i, j] = pair.map(k => Math.floor(k / M.size));
      octx.fillRect(j, 0, 1, M.count);
      octx.fillRect(0, i, M.count, 1);
    }
  }

  // ---- 3D
  const viewer = $3Dmol.createViewer(document.getElementById('mol'), { backgroundColor: 'white' });
  viewer.addModel(data.model, 'cif');
  function residueOf(atom) { return idx.byKey[atom.chain + ':' + atom.resi + (atom.icode || '').trim()]; }
  function colourOf(atom) {
    const k = residueOf(atom);
    if (k === undefined) return '#c7ccd2';
    const e = data.residues[k].e;
    return e ? idx.byElement[e].colour : '#9aa3ad';
  }
  function baseStyle() {
    viewer.setStyle({}, { cartoon: { colorfunc: colourOf } });
    viewer.setStyle({ hetflag: true }, { stick: { radius: 0.15, colorscheme: 'grayCarbon' } });
  }
  baseStyle();
  viewer.zoomTo();
  viewer.render();
  function molOn(residues) {
    baseStyle();
    if (residues.length) {
      const sel = { or: residues.map(k => ({ chain: data.residues[k].c, resi: data.residues[k].n })) };
      viewer.setStyle(sel, { cartoon: { color: '#ff6a00' }, stick: { radius: 0.2, color: '#ff6a00' } });
    }
    viewer.render();
  }

  // ---- linking
  function range(a, b) { const out = []; for (let k = a; k <= b; k++) out.push(k); return out; }
  function name(k) { const r = data.residues[k]; return r.c + ' ' + r.r + r.n + (r.i || ''); }
  function show(state) {
    svgOn(state.keys || []);
    mapOn(state.ranges || [], state.pair || null);
    molOn(state.residues || []);
    info.textContent = state.text || info.dataset.idle;
  }
  function element(id) {
    const e = idx.byElement[id];
    show({ keys: [id], ranges: [[e.start, e.end]], residues: range(e.start, e.end),
           text: e.label + ' — ' + e.chain + ' ' + e.first + '–' + e.last + ' (' + (e.end - e.start + 1) + ' residues)' });
  }
  function clear() { show({}); }
  info.dataset.idle = info.textContent;

  svg.addEventListener('mouseover', ev => {
    for (let n = ev.target; n && n !== svg; n = n.parentNode) {
      if (!n.id) continue;
      const m = n.id.match(/^(strand|helix|helix-back|eta|label|ghost):(.+)$/);
      if (m && idx.byElement[m[2]]) return element(m[2]);
      if (n.id.startsWith('loop:')) {
        const [a, b] = n.id.slice(5).split('>'), s = lib.loopSpan(data, a, b);
        if (s && s[1] >= s[0]) return show({ keys: ['loop ' + a + '>' + b], ranges: [s], residues: range(s[0], s[1]),
                                             text: 'loop ' + name(s[0]) + ' – ' + name(s[1]) });
      }
    }
  });
  svg.addEventListener('mouseleave', clear);
  over.addEventListener('mousemove', ev => {
    const box = over.getBoundingClientRect();
    const J = Math.floor((ev.clientX - box.left) / box.width * M.count), I = Math.floor((ev.clientY - box.top) / box.height * M.count);
    const i = Math.min(data.ca.length - 1, I * M.size), j = Math.min(data.ca.length - 1, J * M.size);
    const keys = [lib.elementAt(data, i), lib.elementAt(data, j)].filter(Boolean);
    show({ keys, pair: [i, j], residues: [i, j],
           text: name(i) + ' · ' + name(j) + ' — ' + lib.distance(data.ca, i, j).toFixed(1) + ' Å' +
                 (lib.distance(data.ca, i, j) <= CUT ? ' (contact)' : '') });
  });
  over.addEventListener('mouseleave', clear);
  viewer.setHoverable({}, true, atom => {
    const k = residueOf(atom);
    if (k === undefined) return;
    const e = data.residues[k].e;
    if (e) element(e); else show({ residues: [k], ranges: [[k, k]], text: name(k) + ' (loop)' });
  }, clear);
})();
"""

_CSS = """
:root { --bg: #f4f6f8; --panel: #ffffff; --ink: #1d2733; --muted: #5d6877; --line: #d8dee5; --accent: #ff6a00;
        --sans: "IBM Plex Sans", "Helvetica Neue", Arial, sans-serif; --mono: "IBM Plex Mono", Menlo, monospace; }
@media (prefers-color-scheme: dark) { :root { --bg: #11161c; --panel: #1a2129; --ink: #e5eaf0; --muted: #98a3b1;
        --line: #2b343f; color-scheme: dark } }
* { box-sizing: border-box }
body { margin: 0; background: var(--bg); color: var(--ink); font: 14px/1.45 var(--sans); padding: 16px }
header { display: flex; flex-wrap: wrap; align-items: baseline; gap: 4px 16px; margin-bottom: 12px }
h1 { font-size: 18px; font-weight: 600; margin: 0 }
header p { margin: 0; color: var(--muted) }
.grid { display: grid; grid-template-columns: minmax(0, 3fr) minmax(0, 2fr); gap: 12px }
@media (max-width: 900px) { .grid { grid-template-columns: minmax(0, 1fr) } }
.panel { background: var(--panel); border: 1px solid var(--line); border-radius: 6px; padding: 10px; min-width: 0 }
.panel h2 { font: 600 12px var(--mono); text-transform: uppercase; letter-spacing: .06em; color: var(--muted); margin: 0 0 8px }
#topology { grid-row: span 2; overflow: auto }
#topology svg { width: 100%; height: auto; display: block; background: #fff; border-radius: 4px }
#topology svg [id^="strand:"], #topology svg [id^="helix:"], #topology svg [id^="loop:"] { cursor: pointer }
#topology svg.focus [id^="strand:"]:not(.on), #topology svg.focus [id^="helix:"]:not(.on),
#topology svg.focus [id^="helix-back:"]:not(.on), #topology svg.focus [id^="eta:"]:not(.on), #topology svg.focus [id^="loop:"]:not(.on),
#topology svg.focus [id^="loop-arrow:"]:not(.on), #topology svg.focus [id^="ghost:"]:not(.on) { opacity: .22; transition: opacity .12s }
#mol { position: relative; width: 100%; height: 420px; border-radius: 4px; overflow: hidden }
.map { position: relative; width: 100%; max-width: 520px; aspect-ratio: 1 }
.map canvas { position: absolute; inset: 0; width: 100%; height: 100%; image-rendering: pixelated }
#map-over { cursor: crosshair }
#strip { display: block; width: 100%; max-width: 520px; height: 10px; image-rendering: pixelated; margin-top: 4px }
#info { font: 13px var(--mono); padding: 8px 10px; background: var(--panel); border: 1px solid var(--line);
        border-radius: 6px; margin-top: 12px; min-height: 2.6em }
.legend { color: var(--muted); font-size: 12px; margin-top: 6px }
"""


def _residues(bb, sses) -> list[dict]:
    owner = {}
    for s in sses:
        for k in range(s.start, s.end + 1):
            owner[k] = s.id
    return [{"c": l.chain, "n": l.seq, "i": l.icode, "r": l.name.title(), "e": owner.get(k)}
            for k, l in enumerate(bb.labels)]


def _model_cif(path, assembly: str) -> str:
    model, _ = read_model(path, assembly)
    st = gemmi.Structure()
    st.add_model(model)
    st.remove_waters()
    st.setup_entities()
    return st.make_mmcif_document().as_string()


def build_page(path, mode: str = "projected", look: Style | None = None, title: str | None = None,
               **layout_options) -> str:
    """The HTML page as a string. layout_options go to cli.make_layout (symmetry, assembly, swap, ...)."""
    from .cli import make_layout

    look = look or Style()
    layout, sses, bb = make_layout(path, mode, look=look, **layout_options)
    loops = route_loops(layout, sses, bb)
    svg = save_svg(draw(layout, loops, sses, title, look=look))
    colours, _ = element_colours(layout, sses, look)
    data = {
        "title": title or Path(path).stem,
        "residues": _residues(bb, sses),
        "ca": np.round(bb.ca, 2).tolist(),
        "elements": [{"id": s.id, "kind": s.kind, "chain": s.chain, "label": layout.placed[s.id].label,
                      "first": s.first.seq, "last": s.last.seq, "start": s.start, "end": s.end,
                      "colour": colours[s.id].lower()} for s in sses if s.id in layout.placed],
        "model": _model_cif(path, layout_options.get("assembly", "auto")),
    }
    blob = json.dumps(data, separators=(",", ":")).replace("</", "<\\/")
    name = html.escape(data["title"])
    n_res, n_el = len(data["residues"]), len(data["elements"])
    chains = len({r["c"] for r in data["residues"]})
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>{name} · topology explorer</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;600&family=IBM+Plex+Mono&display=swap">
<style>{_CSS}</style>
<script src="{THREEDMOL}"></script>
</head><body>
<header><h1>{name}</h1><p>{chains} chain(s) · {n_res} residues · {n_el} helices and strands. Hover anything: the same residues light up in all three views.</p></header>
<div class="grid">
  <section class="panel" id="topology"><h2>Topology</h2>{svg}</section>
  <section class="panel"><h2>3D model</h2><div id="mol"></div></section>
  <section class="panel"><h2>Contact map · CA–CA distance</h2>
    <div class="map"><canvas id="map"></canvas><canvas id="map-over"></canvas></div>
    <canvas id="strip"></canvas>
    <p class="legend">Dark: close (contacts under 8 Å), pale: 30 Å or more. The strip shows each element in its figure colour.</p></section>
</div>
<div id="info">Hover a helix, strand or loop in the topology, a cell of the contact map, or an atom in 3D.</div>
<script type="application/json" id="topo-data">{blob}</script>
<script id="topo-lib">{_LIB}</script>
<script id="topo-app">{_APP}</script>
</body></html>
"""


def write_page(path, out, **kw) -> Path:
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(build_page(path, **kw))
    return out
