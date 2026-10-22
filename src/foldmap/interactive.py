# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""One HTML page linking three views of a structure: the topology figure, the CA distance (contact) map and
the 3D model. Hovering an element, a map cell or an atom highlights the same residues in the other two views.

The page carries its own data, figure and code; only the 3D viewer (Mol* or 3Dmol.js) loads from a CDN, pinned
to one version and checked by Subresource Integrity, so the 3D panel needs a network connection."""

from __future__ import annotations

import html
import json
import re
from pathlib import Path

import gemmi
import numpy as np

from . import CREDIT
from .io import read_model
from .render import INTERACTION_ONE, draw, element_colours, save_svg
from .route import route_loops
from .style import Style

_NOTICE_URL = "https://github.com/dzyla/foldmap/blob/main/NOTICE"
_CREDIT_HTML = CREDIT.replace("see NOTICE", f'see <a href="{_NOTICE_URL}">NOTICE</a>')

THREEDMOL = "https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.1.0/3Dmol-min.js"
MOLSTAR_JS = "https://cdn.jsdelivr.net/npm/pdbe-molstar@3.12.0/build/pdbe-molstar-plugin.js"  # Mol* (PDBe build)
MOLSTAR_CSS = "https://cdn.jsdelivr.net/npm/pdbe-molstar@3.12.0/build/pdbe-molstar-light.css"
# Subresource Integrity: the browser refuses a CDN file whose SHA-384 differs. Recompute when a URL changes:
#   curl -fsSL -o f URL; openssl dgst -sha384 -binary f | openssl base64 -A  (hash a saved file)
SRI = {
    THREEDMOL: "sha384-TxfnTy2TY3ms7kc3hVNjePx2ekn5oouyllBzqjsaZMmtran9oNs+yTKd0TY+x0eY",
    MOLSTAR_JS: "sha384-Whyx5/jg1RTZAEdrf9uj7bdQcSHC2EVtMGgnXGVdgCOGxLbhf9XC4k8brO4hV6yp",
    MOLSTAR_CSS: "sha384-kARlbziDLc0Hj8qk+KXR1P27YYVjD71TknAF2jm8h4pCtySxFh7ivRH8IXXJNC/B",
}
VIEWERS = ("molstar", "3dmol")
_LOOP_GREY = "#9aa3ad"
# The page's click targets (Q-N3): along each nucleic-acid run and stem, and each protein loop, an invisible stroke
# HIT_PX screen pixels wide at any zoom (vector-effect: non-scaling-stroke), so a click up to HIT_PX / 2 px from the
# drawn line's centre selects it; a backbone drawn wider than that (na_line, zoomed in) is its own target. 8, not 6:
# a 6 px stroke ends exactly 3 px from the centre, where a click 3 px off the line would fall on its edge.
HIT_PX = 8
# the stroked parts of a run (na:run:K) and of a stem (na:stem:K:rails, :rung(s), :back, :front, :block, :N for a
# second colour) or a bent duplex drawn as a curve (na:curve:K:..., Style.na_bend curved), each a group whose first
# child is the path matplotlib drew; a stem's texts (:first) are not ink
_NA_INK = re.compile(
    r'<g id="na:(?:(run):(\d+)|(stem|curve):(\d+):(?:rails|rungs?|back|front|block)(?::\d+)?)">\s*'
    r'<path d="([^"]*)"(?: clip-path="([^"]*)")?'
)
_LOOP_INK = re.compile(r'<g id="loop:([^"]+)">\s*<path d="([^"]*)"(?: clip-path="([^"]*)")?')  # a loop's whole line
# an interaction line (hbond:<a>-<b>, ia:<kind>:<a>-<b>): its dashed path, in pieces where it passes a letter
_LINE_INK = re.compile(r'<g id="((?:hbond|ia:[a-z-]+):[^"]+)">\s*<path d="([^"]*)"(?: clip-path="([^"]*)")?')
# the first group drawn at loop level or above (loops and their casings, nucleic-acid ink): the targets go just
# under it, above the panels, plates and membrane band (which would hide them) and under every line a reader clicks
_LOOP_LEVEL = re.compile(r'<g id="(?:loopcase:|loop:|loop-seg:|loop-arrow:|na:)')

# DOM-free helpers: unit-tested under node (tests/test_interactive.py)
_LIB = r"""
(function (g) {
  const TopoLib = {
    itemKey(item) {  // one string per selectable thing
      switch (item.type) {
        case 'element': return 'element:' + item.id;
        case 'loop': return 'loop:' + item.a + '>' + item.b;
        case 'chain': return 'chain:' + item.chain;
        case 'pair': return 'pair:' + item.i + '-' + item.j;
        case 'cell': return 'cell:' + item.i + '-' + item.j;
        case 'disulfide': case 'glycan': case 'ligand': case 'feature': case 'nucleic': return item.type + ':' + item.key;
        case 'interaction': return item.key;  // its gid: hbond:A:44-1:4, ia:salt:A:19-1:7
        case 'bead': return 'bead:' + (item.mark || item.chain + ':' + item.idx);
        default: return 'residue:' + item.k;
      }
    },
    nucleicKey(id) {  // na:stem:3, na:stem:3:rails, na:stem:3:hit, na:label:3 -> 'stem:3'; na:run:5, na:run:5:hit
      // (the page's wide click target) -> 'run:5'; na:curve:0, na:curve:0:rails, na:curve:0:label (a bent duplex drawn
      // as one curve) -> 'curve:0'; anything else -> null
      const m = /^na:(stem|label|run|curve):(\d+)(?::|$)/.exec(id || '');
      return m ? (m[1] === 'run' || m[1] === 'curve' ? m[1] + ':' : 'stem:') + m[2] : null;
    },
    toggle(pinned, item, additive) {  // click: pin it (or clear if it is the only pin); shift-click: add/remove
      const key = TopoLib.itemKey(item), has = pinned.some(p => TopoLib.itemKey(p) === key);
      if (additive) return has ? pinned.filter(p => TopoLib.itemKey(p) !== key) : pinned.concat([item]);
      return has && pinned.length === 1 ? [] : [item];
    },
    chainResidues(data, chain) {
      const out = [];
      data.residues.forEach((r, k) => { if (r.c === chain) out.push(k); });
      return out;
    },
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
    nearestElement(data, k) {  // a residue's element or, in a loop, the closest one along its own chain (fewest
      // residues away, a tie to the N-terminal one); null for a residue whose chain has no helix or strand
      const r = data.residues[k];
      if (!r) return null;
      if (r.e) return r.e;
      let best = null, gap = Infinity;
      data.elements.forEach(e => {
        if (e.chain !== r.c) return;
        const d = Math.max(e.start - k, k - e.end);
        if (d < gap || (d === gap && e.start < best.start)) { best = e; gap = d; }
      });
      return best ? best.id : null;
    },
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
  const svg = document.querySelector('#stage svg');
  const stage = document.getElementById('stage');

  // ---- topology: fit to the panel, wheel to zoom, drag to pan
  const full = (svg.getAttribute('viewBox') || '').split(/\s+/).map(Number);
  if (full.length !== 4 || full.some(isNaN)) {
    const w = parseFloat(svg.getAttribute('width')), h = parseFloat(svg.getAttribute('height'));
    full.splice(0, 4, 0, 0, w, h);
  }
  svg.setAttribute('viewBox', full.join(' '));
  svg.removeAttribute('width'); svg.removeAttribute('height');
  svg.setAttribute('preserveAspectRatio', 'xMidYMid meet');
  let view = full.slice();
  function setView(v) { view = v; svg.setAttribute('viewBox', v.join(' ')); }
  function zoom(f, cx, cy) {
    const [x, y, w, h] = view;
    const px = cx === undefined ? x + w / 2 : cx, py = cy === undefined ? y + h / 2 : cy;
    const nw = Math.min(full[2] * 4, Math.max(full[2] / 12, w * f)), k = nw / w;
    setView([px - (px - x) * k, py - (py - y) * k, w * k, h * k]);
  }
  function toSvg(ev) {
    const pt = svg.createSVGPoint(); pt.x = ev.clientX; pt.y = ev.clientY;
    return pt.matrixTransform(svg.getScreenCTM().inverse());
  }
  stage.addEventListener('wheel', ev => { ev.preventDefault(); const p = toSvg(ev); zoom(ev.deltaY > 0 ? 1.15 : 1 / 1.15, p.x, p.y); },
                         { passive: false });
  let drag = null;
  stage.addEventListener('pointerdown', ev => { drag = { x: ev.clientX, y: ev.clientY, v: view.slice(), moved: false }; });
  window.addEventListener('pointermove', ev => {
    if (!drag || drag.up) return;  // released: a quick move right after a click must not pan
    const dx = ev.clientX - drag.x, dy = ev.clientY - drag.y;
    if (Math.abs(dx) + Math.abs(dy) > 3) { drag.moved = true; stage.classList.add('dragging'); }
    if (!drag.moved) return;
    const s = drag.v[2] / stage.clientWidth > drag.v[3] / stage.clientHeight ? drag.v[2] / stage.clientWidth
                                                                             : drag.v[3] / stage.clientHeight;
    setView([drag.v[0] - dx * s, drag.v[1] - dy * s, drag.v[2], drag.v[3]]);
  });
  window.addEventListener('pointerup', () => {
    if (drag) { drag.up = true; setTimeout(() => { drag = null; }, 0); }
    stage.classList.remove('dragging');
  });
  document.getElementById('zin').onclick = () => zoom(1 / 1.3);
  document.getElementById('zout').onclick = () => zoom(1.3);
  document.getElementById('zfit').onclick = () => setView(full.slice());
  stage.addEventListener('dblclick', () => setView(full.slice()));
  const MAX = 700, CUT = 8;

  // ---- topology: element groups by id
  const parts = {};  // element id -> [svg nodes]
  svg.querySelectorAll('[id]').forEach(node => {
    const m = node.id.match(/^(strand|helix|helix-back|eta|label|ghost|ghostlabel|resnum):(.+?)(:[NC])?$/);
    if (m) (parts[m[2]] = parts[m[2]] || []).push(node);
  });
  svg.querySelectorAll('[id^="loop:"], [id^="loop-arrow:"]').forEach(node => {  // a loop's click target (:hit) too
    const [a, b] = node.id.replace(/:hit$/, '').split(':').slice(1).join(':').split('>');
    const key = 'loop ' + a + '>' + b;
    (parts[key] = parts[key] || []).push(node);
  });
  svg.querySelectorAll('[id^="disulfide:"], [id^="ss-dot:"], [id^="disulfide-end:"], [id^="glycan:"], [id^="ligand"], [id^="feature"], [id^="na:"]').forEach(node => {
    let m = node.id.match(/^(disulfide|ss-dot|disulfide-end):(.+)$/), key = null;
    const na = lib.nucleicKey(node.id);
    if (na) key = 'nucleic ' + na;
    if (m) key = 'disulfide ' + m[2];
    else if ((m = node.id.match(/^glycan:(\d+):/))) key = 'glycan ' + m[1];
    else if ((m = node.id.match(/^feature(?:-ring)?:(\d+)$/))) key = 'feature ' + m[1];
    else if ((m = node.id.match(/^ligand(?:-label|-tether)?:(.+?)(?::\d+)?$/))) key = 'ligand ' + m[1];
    if (key) (parts[key] = parts[key] || []).push(node);
  });
  function svgOn(keys) {
    svg.classList.toggle('focus', keys.length > 0);
    svg.querySelectorAll('.on').forEach(n => n.classList.remove('on'));
    keys.forEach(k => (parts[k] || []).forEach(n => n.classList.add('on')));
  }

  // ---- contact map
  const base = document.getElementById('map'), over = document.getElementById('map-over');
  const hasMap = data.ca.length > 0;  // a nucleic-acid-only file has no protein residues: no map
  if (!hasMap) document.getElementById('mappanel').hidden = true;
  const M = hasMap ? lib.matrix(data.ca, MAX) : { count: 0, size: 1, values: [] };
  [base, over].forEach(c => { c.width = M.count; c.height = M.count; });
  const ctx = base.getContext('2d'), img = hasMap ? ctx.createImageData(M.count, M.count) : null;
  for (let p = 0; p < M.values.length; p++) {
    const [r, g, b] = lib.ramp(M.values[p] / 30);
    img.data.set([r, g, b, 255], p * 4);
  }
  if (hasMap) ctx.putImageData(img, 0, 0);
  const strip = document.getElementById('strip');
  strip.width = M.count; strip.height = 1;
  const sctx = strip.getContext('2d');
  data.elements.forEach(e => {  // the sequence strip under the map: each element in its figure colour
    sctx.fillStyle = e.colour;
    sctx.fillRect(Math.floor(e.start / M.size), 0, Math.max(1, Math.ceil((e.end - e.start + 1) / M.size)), 1);
  });
  const octx = over.getContext('2d');
  function mapOn(ranges, pair) {
    if (!hasMap) return;
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

  // ---- 3D: one interface over Mol* (default) or 3Dmol: base colours, highlight, focus, reset, hover callbacks;
  // the camera moves when the reader moves it (or presses reset) and when a click in the figure or on the map selects
  // something (focus: a short glide onto it, the view's orientation kept), never on a hover, a click in 3D or a cleared
  // selection
  const molDiv = document.getElementById('mol');
  const LOOP = '#9aa3ad', PALE = '#dde2e7', LOOP_ON = '#4a4f57', HET = '#e69f00';
  const FOCUS_MS = 250;  // the glide onto a selection: Mol*'s own default, short enough to keep the reader's place
  function colourOfResidue(k) {
    const e = data.residues[k].e;
    return e ? idx.byElement[e].colour : LOOP;
  }
  function runs(residues) {  // contiguous stretches of one chain and one colour, for range-based selections
    const ks = residues.slice().sort((a, b) => a - b), out = [];
    ks.forEach(k => {
      const r = data.residues[k], colour = colourOfResidue(k), last = out[out.length - 1];
      if (last && last.chain === r.c && last.colour === colour && last.endIndex === k - 1) {
        last.end = r.n; last.endIndex = k;
      } else {
        out.push({ chain: r.c, begin: r.n, end: r.n, endIndex: k, colour });
      }
    });
    return out;
  }
  const hooks = {};

  function molstarViewer() {
    const plugin = new PDBeMolstarPlugin();
    const url = URL.createObjectURL(new Blob([data.model], { type: 'text/plain' }));
    const query = run => ({ auth_asym_id: run.chain, beg_auth_seq_id: run.begin, end_auth_seq_id: run.end });
    const baseData = data.elements.map(e => ({ auth_asym_id: e.chain, beg_auth_seq_id: e.first,
                                               end_auth_seq_id: e.last, color: e.colour }));
    let ready = false, busy = false, pending = null, focusCamera = null;
    async function apply(params) {  // Mol* calls are asynchronous: keep only the latest request while one runs
      if (!ready || busy) { pending = params; return; }
      busy = true;
      try { await plugin.visual.select(params); } catch (err) { console.warn(err); }
      busy = false;
      if (pending) { const next = pending; pending = null; apply(next); }
    }
    const base = () => apply({ data: baseData, nonSelectedColor: LOOP });
    // pdbe-molstar 3.12.0's camera-focus-loci behaviour moves the camera on a click: onto the element clicked
    // (clickCenterFocus) and, on nothing, back to the whole model (clickResetCameraOnEmpty, PluginCommands.Camera.
    // Reset). Its bindings, emptied here (render's focusBindings), match no click: the camera moves only when the
    // reader drags, scrolls or presses a reset control, which call the camera manager directly
    const never = { triggers: [], action: '', description: '' };  // Mol*'s Binding.Empty
    const focusBindings = { clickCenterFocus: never, clickCenterFocusSelectMode: never,
                            clickResetCameraOnEmpty: never, clickResetCameraOnEmptySelectMode: never };
    plugin.render(molDiv, {
      focusBindings,
      customData: { url, format: 'mmcif', binary: false }, bgColor: '#ffffff', hideControls: true,
      sequencePanel: false, leftPanel: false, pdbeLink: false, expanded: false, landscape: false,
      loadingOverlay: false, visualStyle: 'cartoon', hideStructure: ['water'],
      hideCanvasControls: ['expand', 'selection', 'animation', 'controlToggle', 'controlInfo'],
    });
    plugin.events.loadComplete.subscribe(ok => {
      if (!ok) return;
      ready = true;
      const cam = plugin.plugin && plugin.plugin.managers && plugin.plugin.managers.camera;
      // the camera manager's own focus, kept for a click in the figure; Mol*'s callers get none (no zoom on a click in 3D)
      if (cam && cam.focusLoci && !focusCamera) focusCamera = cam.focusLoci.bind(cam);
      if (cam) ['focusLoci', 'focusSpheres'].forEach(f => { if (cam[f]) cam[f] = () => {}; });
      base();
    });
    function key(d) { return idx.byKey[d.auth_asym_id + ':' + d.auth_seq_id + ((d.ins_code || '').trim())]; }
    molDiv.addEventListener('PDB.molstar.mouseover', ev => {
      const k = key(ev.eventData || {});
      if (k !== undefined && hooks.hover) hooks.hover(k);
    });
    molDiv.addEventListener('PDB.molstar.mouseout', () => hooks.leave && hooks.leave());
    molDiv.addEventListener('PDB.molstar.click', ev => {
      const k = key(ev.eventData || {});
      if (k !== undefined && hooks.click) hooks.click(k);
    });
    new ResizeObserver(() => { const c = plugin.plugin && plugin.plugin.canvas3d; if (c) c.handleResize(); })
      .observe(molDiv);
    return {
      name: 'molstar',
      on(residues, het = []) {
        if (!residues.length && !het.length) return base();
        const parts = runs(residues).map(r => Object.assign(query(r), { color: r.colour === LOOP ? LOOP_ON : r.colour,
                                                                         sideChain: true }));
        het.forEach(([chain, seq]) => parts.push({ auth_asym_id: chain, auth_seq_id: seq, color: HET,
                                                   representation: 'ball-and-stick', representationColor: HET }));
        apply({ data: parts, nonSelectedColor: PALE });
      },
      focus(residues, het = []) {  // centre and frame these residues: no turn, so the reader's orientation stays
        if (!ready || !focusCamera || (!residues.length && !het.length)) return;
        const params = runs(residues).map(query).concat(het.map(([chain, seq]) => ({ auth_asym_id: chain,
                                                                                         auth_seq_id: seq })));
        try { focusCamera(plugin.getLociForParams(params), { durationMs: FOCUS_MS }); } catch (err) { console.warn(err); }
      },
      reset() { if (ready) plugin.visual.reset({ camera: true }); },
      plugin,
    };
  }

  function threeDmolViewer() {
    const viewer = $3Dmol.createViewer(molDiv, { backgroundColor: 'white', antialias: true });
    viewer.addModel(data.model, 'cif');
    const residueOf = atom => idx.byKey[atom.chain + ':' + atom.resi + (atom.icode || '').trim()];
    const colourOf = atom => { const k = residueOf(atom); return k === undefined ? '#c7ccd2' : colourOfResidue(k); };
    function base() {
      viewer.setStyle({}, { cartoon: { colorfunc: colourOf } });
      viewer.setStyle({ hetflag: true }, { stick: { radius: 0.15, colorscheme: 'grayCarbon' } });
    }
    base(); viewer.zoomTo(); viewer.render();
    new ResizeObserver(() => { viewer.resize(); viewer.render(); }).observe(molDiv);
    setTimeout(() => { viewer.resize(); viewer.zoomTo(); viewer.render(); }, 50);
    viewer.setHoverable({}, true, atom => {
      const k = residueOf(atom);
      if (k !== undefined && hooks.hover) hooks.hover(k);
    }, () => hooks.leave && hooks.leave());
    viewer.setClickable({}, true, atom => {
      const k = residueOf(atom);
      if (k !== undefined && hooks.click) hooks.click(k);
    });
    const sel = (residues, het = []) => ({ or: residues.map(k => ({ chain: data.residues[k].c, resi: data.residues[k].n }))
                                         .concat(het.map(([chain, seq]) => ({ chain, resi: seq }))) });
    return {
      name: '3dmol',
      on(residues, het = []) {  // like the topology: the selection in its figure colours, everything else pale
        if (!residues.length && !het.length) { base(); viewer.render(); return; }
        viewer.setStyle({}, { cartoon: { color: PALE, opacity: 0.85 } });
        if (residues.length) {
          viewer.setStyle(sel(residues), { cartoon: { colorfunc: colourOf }, stick: { radius: 0.18, colorscheme: 'grayCarbon' } });
        }
        het.forEach(([chain, seq]) => viewer.setStyle({ chain, resi: seq }, { stick: { radius: 0.25, color: HET },
                                                                              sphere: { scale: 0.3, color: HET } }));
        viewer.render();
      },
      focus(residues, het = []) {  // centre and frame these residues (3Dmol's zoomTo keeps the rotation)
        if (!residues.length && !het.length) return;
        viewer.zoomTo(sel(residues, het), FOCUS_MS); viewer.render();
      },
      reset() { viewer.zoomTo({}, 400); viewer.render(); },
      viewer,
    };
  }

  const mol = data.viewer === '3dmol' ? threeDmolViewer() : molstarViewer();
  window.foldmapViewer = mol;
  document.getElementById('molreset').onclick = () => mol.reset();
  const molOn = (residues, het) => mol.on(residues, het || []);

  // ---- linking: hover previews, click pins (shift-click adds), chains via legend or termini
  function range(a, b) { const out = []; for (let k = a; k <= b; k++) out.push(k); return out; }
  function name(k) { const r = data.residues[k]; return r.c + ' ' + r.r + r.n + (r.i || ''); }
  function runsOf(ks) {  // index runs as [start, end] pairs, for the contact map
    const s = ks.slice().sort((a, b) => a - b), out = [];
    s.forEach(k => { const last = out[out.length - 1]; if (last && last[1] === k - 1) last[1] = k; else out.push([k, k]); });
    return out;
  }
  const chainParts = {};  // chain -> svg nodes of its legend entry and termini (lit when the chain is selected)
  svg.querySelectorAll('[id^="legend-chain"], [id^="terminus:"], [id^="stub:"]').forEach(n => {
    const c = n.id.split(':').pop();
    (chainParts[c] = chainParts[c] || []).push(n);
  });
  function stateOf(item) {
    if (item.type === 'element') {
      const e = idx.byElement[item.id];
      return { keys: [item.id], residues: range(e.start, e.end), label: e.label + ' (' + e.chain + ')',
               text: e.label + ' — ' + e.chain + ' ' + e.first + '–' + e.last + ' (' + (e.end - e.start + 1) + ' residues)' };
    }
    if (item.type === 'loop') {
      const s = lib.loopSpan(data, item.a, item.b);
      return { keys: ['loop ' + item.a + '>' + item.b], residues: range(s[0], s[1]), label: 'loop ' + name(s[0]),
               text: 'loop ' + name(s[0]) + ' – ' + name(s[1]) };
    }
    if (item.type === 'chain') {
      const ks = lib.chainResidues(data, item.chain);
      const keys = Object.keys(parts).filter(k => k.replace(/^loop /, '').startsWith(item.chain + ':'));
      const n = data.elements.filter(e => e.chain === item.chain).length;
      return { keys, chain: item.chain, residues: ks, label: 'chain ' + item.chain,
               text: 'chain ' + item.chain + ' — ' + ks.length + ' residues, ' + n + ' helices and strands' };
    }
    if (item.type === 'pair') {
      const d = lib.distance(data.ca, item.i, item.j);
      return { keys: [lib.elementAt(data, item.i), lib.elementAt(data, item.j)].filter(Boolean), pair: [item.i, item.j],
               residues: [item.i, item.j], label: name(item.i) + '·' + name(item.j),
               text: name(item.i) + ' · ' + name(item.j) + ' — ' + d.toFixed(1) + ' Å' + (d <= CUT ? ' (contact)' : '') };
    }
    if (item.type === 'cell') {  // a click on the map: the whole element of each residue, or the nearest one along its chain
      const d = lib.distance(data.ca, item.i, item.j), keys = [], res = new Set(), names = [];
      [item.i, item.j].forEach(k => {
        const id = lib.nearestElement(data, k);
        if (!id) {  // a chain with no helix or strand: the residue itself
          if (!res.has(k)) names.push(name(k) + ' (loop)');
          res.add(k);
          return;
        }
        if (keys.includes(id)) return;
        const e = idx.byElement[id];
        keys.push(id);
        range(e.start, e.end).forEach(x => res.add(x));
        names.push(e.label + ' (' + e.chain + ')');
      });
      return { keys, residues: Array.from(res), pair: [item.i, item.j], label: names.join('·'),
               text: names.join(' · ') + ' — from ' + name(item.i) + ' · ' + name(item.j) + ', ' + d.toFixed(1) + ' Å'
                     + (d <= CUT ? ' (contact)' : '') };
    }
    if (item.type === 'bead') {  // a nucleotide the protein touches: who touches it, and how closely
      const x = item.mark ? data.beads.marks[item.mark] : data.beads.chains[item.chain][item.idx];
      return { keys: [], residues: [x.k], label: x.text, text: x.text };
    }
    if (item.type === 'interaction') {  // an H-bond or other interaction line: its atoms, distances and angles
      const x = data.interactions.find(h => h.key === item.key);
      return { keys: [], residues: x.residues, het: x.het, label: x.label, text: x.text };
    }
    if (['disulfide', 'glycan', 'ligand', 'feature', 'nucleic'].includes(item.type)) {
      const list = data.links[item.type === 'disulfide' ? 'disulfides' : item.type + 's'];
      const x = list.find(d => d.key === item.key);
      const keys = [item.type + ' ' + item.key].concat(x.residues.map(k => data.residues[k].e).filter(Boolean));
      return { keys: ['ligand', 'nucleic'].includes(item.type) ? [keys[0]] : keys, residues: x.residues, het: x.het || [],
               label: item.type === 'ligand' ? x.text.split(',')[0] : item.type === 'nucleic' ? x.text.split(' · ')[0] : item.type,
               text: x.text };
    }
    return { keys: [], residues: [item.k], label: name(item.k), text: name(item.k) + ' (loop)' };
  }
  let pinned = [], hovered = null, shiftDown = false;
  function gather(states) {  // the union of some selections: figure keys, residues, het groups, a map pair
    const keys = new Set(), res = new Set();
    let pair = null;
    const het = [];
    states.forEach(st => { st.keys.forEach(k => keys.add(k)); st.residues.forEach(k => res.add(k)); if (st.pair) pair = st.pair;
                           (st.het || []).forEach(h => het.push(h)); });
    return { keys: Array.from(keys), res: Array.from(res), het, pair };
  }
  function render() {
    const items = hovered ? pinned.concat([hovered]) : pinned, states = items.map(stateOf);
    const { keys, res, het, pair } = gather(states);
    svgOn(keys);
    states.forEach(st => { if (st.chain) (chainParts[st.chain] || []).forEach(n => n.classList.add('on')); });
    mapOn(runsOf(res), pair);
    molOn(res, het);
    let text = info.dataset.idle;
    if (states.length === 1) text = states[0].text;
    else if (states.length > 1) text = states.length + ' selected: ' + states.map(st => st.label).join(', ') + ' — ' + res.length + ' residues';
    if (pinned.length) text += '   ·   click again or Esc to clear, shift-click to add';
    info.textContent = text;
    info.classList.toggle('pinned', pinned.length > 0);
  }
  function hover(item) { hovered = item; render(); }
  function unhover() { if (hovered) { hovered = null; render(); } }
  function pick(item, additive, focus) {
    pinned = lib.toggle(pinned, item, additive);
    hovered = null;
    render();
    // a click in the figure or on the map that selects (or adds) something takes the 3D camera to the whole
    // selection; a click that clears or removes moves nothing (no zoom-out), and nor does a click in the 3D view
    const key = lib.itemKey(item);
    if (focus && pinned.some(p => lib.itemKey(p) === key)) {
      const { res, het } = gather(pinned.map(stateOf));
      mol.focus(res, het);
    }
  }
  function clearPins() { pinned = []; hovered = null; render(); }
  window.foldmapSelection = () => pinned.map(lib.itemKey);  // for tests and scripting
  info.dataset.idle = info.textContent;
  window.addEventListener('keydown', ev => { if (ev.key === 'Shift') shiftDown = true; if (ev.key === 'Escape') clearPins(); });
  window.addEventListener('keyup', ev => { if (ev.key === 'Shift') shiftDown = false; });

  function itemAt(target) {  // the selectable thing under the pointer in the topology, if any
    for (let n = target; n && n !== svg; n = n.parentNode) {
      if (!n.id) continue;
      const m = n.id.match(/^(strand|helix|helix-back|eta|label|ghost|resnum):(.+?)(:[NC])?$/);
      if (m && idx.byElement[m[2]]) return { type: 'element', id: m[2] };
      if (n.id.startsWith('loop:') || n.id.startsWith('loop-seg:')) {
        const body = n.id.startsWith('loop:') ? n.id.slice(5).replace(/:hit$/, '') : n.id.slice(9).replace(/:\d+$/, '');
        const [a, b] = body.split('>'), s = lib.loopSpan(data, a, b);
        if (s && s[1] >= s[0]) return { type: 'loop', a, b };
      }
      const b = n.id.match(/^na:contact:(.+)$/), t = n.id.match(/^na:nt:(.+)$/);
      if (b && data.beads.chains[b[1]]) {
        let u = target;
        while (u && u !== n && u.tagName !== 'use') u = u.parentNode;
        const idx = u && u !== n ? Array.prototype.indexOf.call(n.querySelectorAll('use'), u) : -1;
        if (idx >= 0 && data.beads.chains[b[1]][idx]) return { type: 'bead', chain: b[1], idx };
      }
      if (t && data.beads.marks[t[1]]) return { type: 'bead', mark: t[1] };  // the bold letter of a touched nucleotide
      const c = n.id.match(/^(legend-chain|legend-chain-label|terminus:[NC]|stub:[NC]):(.+)$/);
      if (c) return { type: 'chain', chain: c[2] };
      let l = n.id.match(/^(?:disulfide|ss-dot|disulfide-end):(.+)$/);
      if (l) return { type: 'disulfide', key: l[1] };
      if ((l = n.id.match(/^glycan:(\d+):/))) return { type: 'glycan', key: l[1] };
      if ((l = n.id.match(/^((?:hbond|ia:[a-z-]+):.+?)(?::hit)?$/)) && data.interactions.some(h => h.key === l[1])) {
        return { type: 'interaction', key: l[1] };  // a line or its click target (:hit)
      }
      if ((l = n.id.match(/^feature(?:-ring)?:(\d+)$/))) return { type: 'feature', key: l[1] };
      if ((l = n.id.match(/^ligand(?:-label|-tether)?:(.+?)(?::\d+)?$/))) return { type: 'ligand', key: l[1] };
      const na = lib.nucleicKey(n.id);
      if (na && data.links.nucleics.some(x => x.key === na)) return { type: 'nucleic', key: na };
    }
    return null;
  }
  svg.addEventListener('mouseover', ev => { const item = itemAt(ev.target); if (item) hover(item); else unhover(); });
  svg.addEventListener('mouseleave', unhover);
  svg.addEventListener('click', ev => {
    if (drag && drag.moved) return;
    const item = itemAt(ev.target);
    if (item) pick(item, ev.shiftKey, true);
    else if (!ev.shiftKey) clearPins();  // a click on empty paper clears the selection
  });
  svg.querySelectorAll('[id^="legend-chain"], [id^="terminus:"]').forEach(n => { n.style.cursor = 'pointer'; });

  const mapBox = document.querySelector('.map'), stripEl = document.getElementById('strip');
  new ResizeObserver(() => { stripEl.style.width = mapBox.clientWidth + 'px'; }).observe(mapBox);
  function pairAt(ev) {
    const box = over.getBoundingClientRect();
    const J = Math.floor((ev.clientX - box.left) / box.width * M.count), I = Math.floor((ev.clientY - box.top) / box.height * M.count);
    return { type: 'pair', i: Math.min(data.ca.length - 1, I * M.size), j: Math.min(data.ca.length - 1, J * M.size) };
  }
  over.addEventListener('mousemove', ev => hover(pairAt(ev)));
  over.addEventListener('mouseleave', unhover);
  // a hover shows the residue pair; a click selects the whole elements the two residues belong to (or lie nearest)
  over.addEventListener('click', ev => { const p = pairAt(ev); pick({ type: 'cell', i: p.i, j: p.j }, ev.shiftKey, true); });

  const residueItem = k => (data.residues[k].e ? { type: 'element', id: data.residues[k].e } : { type: 'residue', k });
  hooks.hover = k => hover(residueItem(k));  // the 3D model
  hooks.leave = unhover;
  hooks.click = k => pick(residueItem(k), shiftDown);
})();
"""

_CSS = """
:root { --bg: #eef1f4; --panel: #ffffff; --ink: #1d2733; --muted: #5d6877; --line: #d5dce3; --accent: #ff6a00;
        --sans: "IBM Plex Sans", "Helvetica Neue", Arial, sans-serif; --mono: "IBM Plex Mono", Menlo, monospace;
        color-scheme: light }
* { box-sizing: border-box }
html, body { height: 100% }
body { margin: 0; background: var(--bg); color: var(--ink); font: 14px/1.45 var(--sans); display: flex;
       flex-direction: column; overflow: hidden }
header { display: flex; flex-wrap: wrap; align-items: baseline; gap: 2px 16px; padding: 10px 16px 6px }
h1 { font-size: 17px; font-weight: 600; margin: 0 }
header p { margin: 0; color: var(--muted); font-size: 13px }
main { flex: 1; min-height: 0; display: grid; grid-template-columns: minmax(0, 1fr) minmax(320px, 38%); gap: 10px;
       padding: 0 16px }
.side { display: flex; flex-direction: column; gap: 10px; min-height: 0 }
.panel { background: var(--panel); border: 1px solid var(--line); border-radius: 8px; padding: 8px 10px; min-width: 0;
         min-height: 0; display: flex; flex-direction: column; position: relative }
.panel h2 { font: 600 11px var(--mono); text-transform: uppercase; letter-spacing: .06em; color: var(--muted);
            margin: 0 0 6px; display: flex; align-items: center; gap: 8px }
.panel h2 .tools { margin-left: auto; display: flex; gap: 4px }
.tools button { font: 12px var(--mono); border: 1px solid var(--line); background: #f7f9fb; color: var(--ink);
                border-radius: 5px; padding: 1px 8px; cursor: pointer }
.tools button:hover { border-color: var(--accent) }
#topology { min-height: 0 }
#stage { flex: 1; min-height: 0; overflow: hidden; border-radius: 6px; background: #fff; cursor: grab;
         touch-action: none }
#stage.dragging { cursor: grabbing }
#stage svg { width: 100%; height: 100%; display: block }
#stage svg [id^="strand:"], #stage svg [id^="helix:"], #stage svg [id^="eta:"], #stage svg [id^="loop:"] { cursor: pointer }
#stage svg.focus [id^="strand:"]:not(.on), #stage svg.focus [id^="helix:"]:not(.on),
#stage svg.focus [id^="helix-back:"]:not(.on), #stage svg.focus [id^="eta:"]:not(.on),
#stage svg.focus [id^="loop:"]:not(.on), #stage svg.focus [id^="loop-arrow:"]:not(.on),
#stage svg.focus [id^="ghost:"]:not(.on),
#stage svg.focus [id^="disulfide:"]:not(.on), #stage svg.focus [id^="ss-dot:"]:not(.on),
#stage svg.focus [id^="disulfide-end:"]:not(.on),
#stage svg.focus [id^="glycan:"]:not(.on), #stage svg.focus [id^="ligand"]:not(.on),
#stage svg.focus [id^="feature"]:not(.on), #stage svg.focus [id^="na:stem:"]:not(.on),
#stage svg.focus [id^="na:label:"]:not(.on), #stage svg.focus [id^="na:run:"]:not(.on),
#stage svg.focus [id^="na:curve:"]:not(.on) { opacity: .2; transition: opacity .12s }
#stage svg [id^="disulfide:"], #stage svg [id^="ss-dot:"], #stage svg [id^="disulfide-end:"], #stage svg [id^="glycan:"], #stage svg [id^="ligand"],
#stage svg [id^="na:stem:"], #stage svg [id^="na:label:"], #stage svg [id^="na:run:"], #stage svg [id^="na:curve:"],
#stage svg [id^="hbond:"], #stage svg [id^="ia:"] { cursor: pointer }
#stage svg [id^="loopcase:"] { pointer-events: none }  /* a loop's white casing: clicks reach its click target */
#molpanel { flex: 1.15 }
#mol { flex: 1; min-height: 200px; position: relative; border-radius: 6px; overflow: hidden; background: #fff }
#mappanel { flex: 1 }
[hidden] { display: none !important }
.mapwrap { flex: 1; min-height: 0; display: flex; flex-direction: column; align-items: center }
.map { position: relative; height: 100%; max-width: 100%; aspect-ratio: 1 }
.map canvas { position: absolute; inset: 0; width: 100%; height: 100%; image-rendering: pixelated; border-radius: 3px }
#map-over { cursor: crosshair }
#strip { display: block; height: 8px; image-rendering: pixelated; margin-top: 4px; border-radius: 2px }
.legend { color: var(--muted); font-size: 11.5px; margin: 4px 0 0 }
#info.pinned { border-color: var(--accent) }
#info { font: 13px var(--mono); padding: 7px 12px; background: var(--panel); border: 1px solid var(--line);
        border-radius: 8px; margin: 10px 16px 0; min-height: 2.4em; white-space: nowrap; overflow: hidden;
        text-overflow: ellipsis }
#credit { padding: 6px 16px 10px; font: 11px var(--mono); color: var(--muted) }
@media (max-width: 900px) {
  body { overflow: auto; display: block }
  main { display: flex; flex-direction: column }
  #topology { height: 75vh } #molpanel { height: 60vh } #mappanel { height: 70vw }
}
"""


def _residues(bb, sses) -> list[dict]:
    owner = {}
    for s in sses:
        for k in range(s.start, s.end + 1):
            owner[k] = s.id
    return [
        {"c": l.chain, "n": l.seq, "i": l.icode, "r": l.name.title(), "e": owner.get(k)}
        for k, l in enumerate(bb.labels)
    ]


def _model_cif(path, assembly: str, axes=None, centre=None, model=None) -> str:
    """The model as mmCIF; with axes (rows: page-right, page-up, toward viewer) turned into page coordinates, so
    the 3D viewer's default camera sees it the way the figure shows it."""
    model = (model or read_model(path, assembly))[0].clone()  # moved below: a given model stays as read
    if axes is not None:
        rot = np.asarray(axes, float)
        shift = -rot @ np.asarray(centre, float)
        model.transform_pos_and_adp(gemmi.Transform(gemmi.Mat33(rot.tolist()), gemmi.Vec3(*shift)))
    st = gemmi.Structure()
    st.add_model(model)
    st.remove_waters()
    st.setup_entities()
    doc = st.make_mmcif_document()
    _type_components(doc.sole_block(), path)
    return doc.as_string()


def _type_components(block, path) -> None:
    """Fill _chem_comp.type (gemmi leaves it '.'): from the source file when it has them, else from gemmi's
    residue table. Viewers such as Mol* tell polymer from ligand by it."""
    known = {}
    try:
        src = gemmi.cif.read(str(path)).sole_block().find("_chem_comp.", ["id", "type"])
        known = {r[0]: gemmi.cif.as_string(r[1]) for r in src}
    except (RuntimeError, ValueError, IndexError):
        pass
    table = block.find("_chem_comp.", ["id", "type"])
    for row in table:
        name = row[0]
        kind = known.get(name) or _component_type(name)
        row[1] = gemmi.cif.quote(kind)


def _component_type(name: str) -> str:
    info = gemmi.find_tabulated_residue(name)
    if info is None:
        return "non-polymer"
    if info.is_amino_acid():
        return "peptide linking" if name == "GLY" else "L-peptide linking"
    if info.is_nucleic_acid():
        return "DNA linking" if name.startswith("D") else "RNA linking"
    return "non-polymer"


def _nucleics(layout) -> list[dict]:
    """Stems (a bent duplex drawn as a curve: one, curve:<n>) and single-stranded stretches as drawn: their nucleotides
    (for the 3D view) and the protein residues touching them."""
    na = layout.nucleic
    if na is None:
        return []
    toucher: dict[tuple[int, int], set[int]] = {}
    for r, nts in na.contacts.items():
        for x in nts:
            toucher.setdefault(x, set()).add(r)

    def span(first, last):  # "C -49 … 49", across chains "C -49 … E 49", one nucleotide "C -50"
        if (first.chain, first.tag) == (last.chain, last.tag):
            return f"{first.chain} {first.tag}"
        return f"{first.chain} {first.tag} \u2026 {'' if last.chain == first.chain else last.chain + ' '}{last.tag}"

    def entry(key, nts, text):
        labels = [na.strands[s].labels[i] for s, i in nts]
        residues = sorted({r for x in nts for r in toucher.get(x, ())})
        return {"key": key, "residues": residues, "het": [[lab.chain, lab.seq] for lab in labels], "text": text}

    out = []
    curves = getattr(layout, "na_curves", None) or []  # a bent duplex drawn as one curve: one entry, not its pieces
    hidden = {k for c in curves for k in c.pieces}
    for k, d in enumerate(layout.stems):
        if k in hidden:
            continue
        keys = sorted(d.axial)
        ends = span(na.strands[keys[0][0]].labels[keys[0][1]], na.strands[keys[-1][0]].labels[keys[-1][1]])
        out.append(entry(f"stem:{k}", keys, f"{d.name} \u00b7 {len(d.pairs)} bp \u00b7 {ends}"))
    for n, c in enumerate(curves):
        keys = sorted(c.at)
        ends = span(na.strands[keys[0][0]].labels[keys[0][1]], na.strands[keys[-1][0]].labels[keys[-1][1]])
        out.append(entry(f"curve:{n}", keys, f"{c.name} \u00b7 {len(c.pairs)} bp \u00b7 {ends}"))
    for n, p in enumerate(layout.na_paths):
        if p.nts:
            lab = na.strands[p.strand].labels
            text = f"{p.kind} {span(lab[p.nts[0]], lab[p.nts[-1]])} ({len(p.nts)} nt)"
            out.append(entry(f"run:{n}", [(p.strand, i) for i in p.nts], text))
    return out


def _links(layout, loops, sses, bb, look) -> dict:
    """Disulfides, glycans and ligands as drawn, each with its protein residues (Backbone indices) and, for
    glycans and ligands, their own residues (chain, number) for the 3D view."""
    from .render import ligand_marks

    links = layout.links
    out = {"disulfides": [], "glycans": [], "ligands": [], "features": [], "nucleics": _nucleics(layout)}
    if links is None and getattr(layout, "uniprot", None) is None:
        return out
    shown = lambda r: not layout.partial or layout.res_chain[r] not in layout.partial  # noqa: E731
    label = lambda r: f"{bb.labels[r].chain} {bb.labels[r].name.title()}{bb.labels[r].seq}"  # noqa: E731
    if look.disulfides and links is not None:
        for i, j in links.disulfides:
            if shown(i) and shown(j):
                out["disulfides"].append(
                    {"key": f"{i}-{j}", "residues": [i, j], "text": f"disulfide {label(i)} – {label(j)}"}
                )
    if look.glycans and links is not None:
        members = getattr(links, "glycan_members", []) or [[] for _ in links.glycans]
        for (r, sugars), het in zip(links.glycans, members):
            if shown(r):
                out["glycans"].append(
                    {
                        "key": str(r),
                        "residues": [r],
                        "het": [list(h) for h in het],
                        "text": f"glycan on {label(r)}: {'-'.join(sugars)}",
                    }
                )
    from .render import uniprot_marks

    for n, m in enumerate(uniprot_marks(layout, loops, sses, look)):
        f, k = m["feature"], m["residue"]
        out["features"].append(
            {"key": str(n), "residues": [k], "text": f"UniProt {f.type.lower()}: {f.label} ({label(k)})"}
        )
    for m in ligand_marks(layout, loops, sses, look):
        g = m["ligand"]
        key = f"{g.name}:{g.chain}{g.seq}"
        names = [label(r) for r in g.contacts] or [f"{c}:{n}{i}" for c, n, i in g.nucleotides]  # nucleotides (C33)
        held = ", ".join(names[:6]) + (" …" if len(names) > 6 else "")
        out["ligands"].append(
            {
                "key": key,
                "residues": list(g.contacts),
                "het": [list(h) for h in g.members or [(g.chain, g.seq)]],
                "text": f"{'metal' if g.metal else 'ligand'} {g.symbol}, held by {held}",
            }
        )
    return out


def _beads(fig, layout, bb, model) -> dict:
    """The contact beads as drawn: `chains`, per chain in drawing order, and `marks`, by nucleotide ("A:DG7",
    the id of its letter), each with the toucher (Backbone index) and the hover note "DG7 – contacted by
    A:Arg42 (2.9 Å)", the distance being the closest heavy-atom pair."""
    from .nacheck import verify_contacts

    drawn = getattr(fig, "foldmap_beads", None)
    out = {"chains": {}, "marks": {}}
    if drawn is None or not drawn.marks:
        return out
    pairs = {b.key: b.pair for b in verify_contacts(model, drawn, layout, bb).beads}
    for chain, items in drawn.marks.items():
        out["chains"][chain] = []
        for key, r, _ in items:
            lab, who = layout.nucleic.strands[key[0]].labels[key[1]], bb.labels[r]
            text = f"{lab.name}{lab.tag} – contacted by {who.chain}:{who.name.title()}{who.tag}"
            pair = pairs.get(key)
            x = {"k": r, "text": text + (f" ({pair.dist:.1f} Å)" if pair else "")}
            out["chains"][chain].append(x)
            out["marks"][f"{lab.chain}:{lab.tag}"] = x
    return out


def _interactions(fig) -> list:
    """The interaction lines as drawn (the H-bonds' hbond: lines and the other kinds' ia: lines), for the hover and
    the click: render.interaction_record's `key` (the line's gid), `kind`, `residues` (Backbone indices of the
    protein ends), `het` (its nucleotides, for the 3D view) and `text`, the note "H-bond A:Gln44 – 1:DA4: NE2 – N7
    3.3 Å; ..." or "Salt bridge A:Arg19 – 1:DT7: NH1 – OP1 3.1 Å" (atoms in the order of the residues, an angle
    where the kind has one), and `label`, what one line of its kind is called (H-bond, Salt bridge)."""
    overlay = getattr(fig, "foldmap_interactions", None)
    keep = ("key", "kind", "residues", "het", "text")
    return [
        {**{k: x[k] for k in keep}, "label": INTERACTION_ONE[x["kind"]]} for x in (overlay.records() if overlay else [])
    ]


def _hit_targets(svg: str) -> str:
    """The page's copy of the figure with click targets: an invisible path HIT_PX screen pixels wide along the drawn
    lines of each nucleic-acid run, stem and curve (na:run:K:hit, na:stem:K:hit, na:curve:K:hit; one per run, stem and
    curve, its rails, rungs and body together), of each protein loop (loop:<name>:hit) and of each interaction line
    (hbond:<a>-<b>:hit, ia:<kind>:<a>-<b>:hit; plan 13c-5). All sit under the loops and the nucleic-acid ink, so a click
    on anything drawn selects that thing; the loops' targets lie over the nucleic acid's, so a click a few pixels off a
    loop that runs beside a strand still selects the loop, and both lie over the interaction lines', which end on them.
    The static figure (save_svg) has none."""
    lines: dict[str, list[str]] = {}
    clips: dict[str, str | None] = {}

    def add(key, d, clip):
        lines.setdefault(key, []).append(d.strip())
        clips.setdefault(key, clip)

    for m in _LINE_INK.finditer(svg):
        add(f"{m[1]}:hit", m[2], m[3])
    for m in _NA_INK.finditer(svg):
        add(f"na:{m[1] or m[3]}:{m[2] or m[4]}:hit", m[5], m[6])
    for m in _LOOP_INK.finditer(svg):
        add(f"loop:{m[1]}:hit", m[2], m[3])
    under = _LOOP_LEVEL.search(svg) or _LINE_INK.search(svg)
    if not lines or under is None:
        return svg
    style = (
        f"fill: none; stroke: transparent; stroke-width: {HIT_PX}px; stroke-linecap: round; stroke-linejoin: round; "
        "vector-effect: non-scaling-stroke; pointer-events: stroke"
    )
    hits = "".join(
        f'<path id="{key}" d="{" ".join(d)}"'
        + (f' clip-path="{clips[key]}"' if clips[key] else "")
        + f' style="{style}"/>\n   '
        for key, d in lines.items()
    )
    return svg[: under.start()] + hits + svg[under.start() :]


def _viewer_head(viewer: str) -> str:
    def script(url: str) -> str:
        return f'<script src="{url}" integrity="{SRI[url]}" crossorigin="anonymous"></script>'

    if viewer == "3dmol":
        return script(THREEDMOL)
    css = f'<link rel="stylesheet" href="{MOLSTAR_CSS}" integrity="{SRI[MOLSTAR_CSS]}" crossorigin="anonymous">'
    return css + "\n" + script(MOLSTAR_JS)


def _script_json(data) -> str:
    """JSON safe inside a <script> element: no character of it can end the element or open a comment."""
    text = json.dumps(data, separators=(",", ":"))
    return text.replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")


def build_page(
    path,
    mode: str = "projected",
    look: Style | None = None,
    title: str | None = None,
    viewer: str = "molstar",
    provenance_out: list | None = None,
    **layout_options,
) -> str:
    """The HTML page as a string. viewer: 'molstar' (Mol*, PDBe build) or '3dmol' for the 3D panel.
    provenance_out: a list the figure's provenance items (how it was made) are added to, for the caller to say.
    layout_options go to cli.make_layout (symmetry, assembly, swap, ...)."""
    from .cli import make_layout

    if viewer not in VIEWERS:
        raise ValueError(f"unknown 3D viewer {viewer!r}; choose from {', '.join(VIEWERS)}")

    look = look or Style()
    model = read_model(path, layout_options.get("assembly", "auto"))  # once: layout, hover notes and 3D model
    layout, sses, bb = make_layout(path, mode, look=look, model=model, **layout_options)
    loops = route_loops(layout, sses, bb)
    fig = draw(layout, loops, sses, title, model=model, look=look)
    svg = _hit_targets(save_svg(fig))
    if provenance_out is not None:
        provenance_out.extend(dict(i) for i in getattr(fig, "foldmap_provenance", []))
    colours, _ = element_colours(layout, sses, look)
    data = {
        "title": title or Path(path).stem,
        "residues": _residues(bb, sses),
        "ca": np.round(bb.ca, 2).tolist(),
        "elements": [
            {
                "id": s.id,
                "kind": s.kind,
                "chain": s.chain,
                "label": layout.placed[s.id].label,
                "first": s.first.seq,
                "last": s.last.seq,
                "start": s.start,
                "end": s.end,
                "colour": colours[s.id].lower(),
            }
            for s in sses
            if s.id in layout.placed
        ],
        "model": _model_cif(path, "", layout.page_axes, layout.page_centre, model),
        "viewer": viewer,
        "links": _links(layout, loops, sses, bb, look),
        "beads": _beads(fig, layout, bb, model),
        "interactions": _interactions(fig),
    }
    blob = _script_json(data)
    name = html.escape(data["title"])
    n_res, n_el = len(data["residues"]), len(data["elements"])
    nucleic = layout.nucleic.strands if layout.nucleic is not None else []
    n_nt = sum(len(st.labels) for st in nucleic)
    chains = len({r["c"] for r in data["residues"]} | {lab.chain for st in nucleic for lab in st.labels})
    counts = " · ".join(
        [f"{chains} chain{'s' if chains != 1 else ''}"]
        + ([f"{n_res} residues"] if n_res or not n_nt else [])
        + ([f"{n_nt} nt"] if n_nt else [])
        + ([f"{n_el} helices and strands"] if n_el or not n_nt else [])
    )
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>{name} · topology explorer</title>
<style>{_CSS}</style>
{_viewer_head(viewer)}
</head><body>
<header><h1>{name}</h1><p>{counts} · hover anything: the
same residues light up in all three views · click an element to keep it lit</p></header>
<main>
  <section class="panel" id="topology"><h2>Topology <span class="tools"><button id="zin" title="zoom in">+</button>
    <button id="zout" title="zoom out">−</button><button id="zfit" title="fit the whole figure">fit</button></span></h2>
    <div id="stage">{svg}</div></section>
  <div class="side">
    <section class="panel" id="molpanel"><h2>3D model <span class="tools"><button id="molreset" title="show the whole
      model">reset view</button></span></h2><div id="mol"></div></section>
    <section class="panel" id="mappanel"><h2>Contact map · CA–CA distance</h2>
      <div class="mapwrap"><div class="map"><canvas id="map"></canvas><canvas id="map-over"></canvas></div>
      <canvas id="strip"></canvas></div>
      <p class="legend">Dark: in contact (under 8 Å) · pale: 30 Å or more · strip: elements in figure colours.</p></section>
  </div>
</main>
<div id="info">Hover a helix, strand or loop in the topology, a cell of the contact map, or an atom in 3D.</div>
<footer id="credit">{_CREDIT_HTML}</footer>
<script type="application/json" id="topo-data">{blob}</script>
<script id="topo-lib">{_LIB}</script>
<script id="topo-app">{_APP}</script>
</body></html>
"""


def write_page(path, out, **kw) -> Path:
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(build_page(path, **kw), encoding="utf-8")
    return out
