"""One self-contained HTML page linking three views of a structure: the topology figure, the CA distance
(contact) map and the 3D model. Hovering an element, a map cell or an atom highlights the same residues in
the other two views."""

from __future__ import annotations

import html
import json
from pathlib import Path

import gemmi
import numpy as np

from . import CREDIT
from .io import read_model
from .render import draw, element_colours, save_svg
from .route import route_loops
from .style import Style

THREEDMOL = "https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.1.0/3Dmol-min.js"
MOLSTAR_JS = "https://cdn.jsdelivr.net/npm/pdbe-molstar@3.12.0/build/pdbe-molstar-plugin.js"  # Mol* (PDBe build)
MOLSTAR_CSS = "https://cdn.jsdelivr.net/npm/pdbe-molstar@3.12.0/build/pdbe-molstar-light.css"
VIEWERS = ("molstar", "3dmol")
_LOOP_GREY = "#9aa3ad"

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
        case 'disulfide': case 'glycan': case 'ligand': return item.type + ':' + item.key;
        default: return 'residue:' + item.k;
      }
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
  svg.querySelectorAll('[id^="loop:"], [id^="loop-arrow:"]').forEach(node => {
    const [a, b] = node.id.split(':').slice(1).join(':').split('>');
    const key = 'loop ' + a + '>' + b;
    (parts[key] = parts[key] || []).push(node);
  });
  svg.querySelectorAll('[id^="disulfide:"], [id^="ss-dot:"], [id^="glycan:"], [id^="ligand"]').forEach(node => {
    let m = node.id.match(/^(disulfide|ss-dot):(.+)$/), key = null;
    if (m) key = 'disulfide ' + m[2];
    else if ((m = node.id.match(/^glycan:(\d+):/))) key = 'glycan ' + m[1];
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

  // ---- 3D: one interface over Mol* (default) or 3Dmol: base colours, highlight, fly-to, reset, hover callbacks
  const molDiv = document.getElementById('mol');
  const LOOP = '#9aa3ad', PALE = '#dde2e7', LOOP_ON = '#4a4f57', HET = '#e69f00';
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
    let ready = false, busy = false, pending = null;
    async function apply(params) {  // Mol* calls are asynchronous: keep only the latest request while one runs
      if (!ready || busy) { pending = params; return; }
      busy = true;
      try { await plugin.visual.select(params); } catch (err) { console.warn(err); }
      busy = false;
      if (pending) { const next = pending; pending = null; apply(next); }
    }
    const base = () => apply({ data: baseData, nonSelectedColor: LOOP });
    plugin.render(molDiv, {
      customData: { url, format: 'mmcif', binary: false }, bgColor: '#ffffff', hideControls: true,
      sequencePanel: false, leftPanel: false, pdbeLink: false, expanded: false, landscape: false,
      loadingOverlay: false, visualStyle: 'cartoon', hideStructure: ['water'],
      hideCanvasControls: ['expand', 'selection', 'animation', 'controlToggle', 'controlInfo'],
    });
    plugin.events.loadComplete.subscribe(ok => { if (ok) { ready = true; base(); } });
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
      on(residues, het) {
        if (!residues.length && !het.length) return base();
        const parts = runs(residues).map(r => Object.assign(query(r), { color: r.colour === LOOP ? LOOP_ON : r.colour,
                                                                         sideChain: true }));
        het.forEach(([chain, seq]) => parts.push({ auth_asym_id: chain, auth_seq_id: seq, color: HET,
                                                   representation: 'ball-and-stick', representationColor: HET }));
        apply({ data: parts, nonSelectedColor: PALE });
      },
      async fly(residues) {  // focus, then step back so the element sits in its surroundings
        if (!ready || !residues.length) return;
        await plugin.visual.focus(runs(residues).map(query));
        setTimeout(() => {
          const c = plugin.plugin && plugin.plugin.canvas3d;
          if (!c) return;
          const snap = c.camera.getSnapshot(), k = 2.4;
          const position = [0, 1, 2].map(i => snap.target[i] + (snap.position[i] - snap.target[i]) * k);
          c.camera.setState({ position, radius: snap.radius * k }, 300);
        }, 450);
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
    const sel = residues => ({ or: residues.map(k => ({ chain: data.residues[k].c, resi: data.residues[k].n })) });
    return {
      name: '3dmol',
      on(residues, het) {  // like the topology: the selection in its figure colours, everything else pale
        if (!residues.length && !het.length) { base(); viewer.render(); return; }
        viewer.setStyle({}, { cartoon: { color: PALE, opacity: 0.85 } });
        if (residues.length) {
          viewer.setStyle(sel(residues), { cartoon: { colorfunc: colourOf }, stick: { radius: 0.18, colorscheme: 'grayCarbon' } });
        }
        het.forEach(([chain, seq]) => viewer.setStyle({ chain, resi: seq }, { stick: { radius: 0.25, color: HET },
                                                                              sphere: { scale: 0.3, color: HET } }));
        viewer.render();
      },
      fly(residues) { if (residues.length) { viewer.zoomTo(sel(residues)); viewer.zoom(0.35); viewer.render(); } },
      reset() { viewer.zoomTo({}, 400); viewer.render(); },
      viewer,
    };
  }

  const mol = data.viewer === '3dmol' ? threeDmolViewer() : molstarViewer();
  window.foldmapViewer = mol;
  document.getElementById('molreset').onclick = () => mol.reset();
  const molOn = (residues, het) => mol.on(residues, het || []);
  const fly = residues => mol.fly(residues);

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
    if (item.type === 'disulfide' || item.type === 'glycan' || item.type === 'ligand') {
      const list = data.links[item.type === 'disulfide' ? 'disulfides' : item.type + 's'];
      const x = list.find(d => d.key === item.key);
      const keys = [item.type + ' ' + item.key].concat(x.residues.map(k => data.residues[k].e).filter(Boolean));
      return { keys: item.type === 'ligand' ? [keys[0]] : keys, residues: x.residues, het: x.het || [],
               label: item.type === 'ligand' ? x.text.split(',')[0] : item.type, text: x.text };
    }
    return { keys: [], residues: [item.k], label: name(item.k), text: name(item.k) + ' (loop)' };
  }
  let pinned = [], hovered = null, shiftDown = false;
  function render() {
    const items = hovered ? pinned.concat([hovered]) : pinned, states = items.map(stateOf);
    const keys = new Set(), res = new Set();
    let pair = null;
    const het = [];
    states.forEach(st => { st.keys.forEach(k => keys.add(k)); st.residues.forEach(k => res.add(k)); if (st.pair) pair = st.pair;
                           (st.het || []).forEach(h => het.push(h)); });
    svgOn(Array.from(keys));
    states.forEach(st => { if (st.chain) (chainParts[st.chain] || []).forEach(n => n.classList.add('on')); });
    mapOn(runsOf(Array.from(res)), pair);
    molOn(Array.from(res), het);
    let text = info.dataset.idle;
    if (states.length === 1) text = states[0].text;
    else if (states.length > 1) text = states.length + ' selected: ' + states.map(st => st.label).join(', ') + ' — ' + res.size + ' residues';
    if (pinned.length) text += '   ·   click again or Esc to clear, shift-click to add';
    info.textContent = text;
    info.classList.toggle('pinned', pinned.length > 0);
  }
  function hover(item) { hovered = item; render(); }
  function unhover() { if (hovered) { hovered = null; render(); } }
  function pick(item, additive) {
    pinned = lib.toggle(pinned, item, additive);
    hovered = null;
    render();
    if (pinned.length) {
      const ks = new Set();
      pinned.forEach(p => stateOf(p).residues.forEach(k => ks.add(k)));
      fly(Array.from(ks));
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
        const body = n.id.startsWith('loop:') ? n.id.slice(5) : n.id.slice(9).replace(/:\d+$/, '');
        const [a, b] = body.split('>'), s = lib.loopSpan(data, a, b);
        if (s && s[1] >= s[0]) return { type: 'loop', a, b };
      }
      const c = n.id.match(/^(legend-chain|legend-chain-label|terminus:[NC]|stub:[NC]):(.+)$/);
      if (c) return { type: 'chain', chain: c[2] };
      let l = n.id.match(/^(?:disulfide|ss-dot):(.+)$/);
      if (l) return { type: 'disulfide', key: l[1] };
      if ((l = n.id.match(/^glycan:(\d+):/))) return { type: 'glycan', key: l[1] };
      if ((l = n.id.match(/^ligand(?:-label|-tether)?:(.+?)(?::\d+)?$/))) return { type: 'ligand', key: l[1] };
    }
    return null;
  }
  svg.addEventListener('mouseover', ev => { const item = itemAt(ev.target); if (item) hover(item); else unhover(); });
  svg.addEventListener('mouseleave', unhover);
  svg.addEventListener('click', ev => {
    if (drag && drag.moved) return;
    const item = itemAt(ev.target);
    if (item) pick(item, ev.shiftKey);
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
  over.addEventListener('click', ev => pick(pairAt(ev), ev.shiftKey));

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
#stage svg.focus [id^="glycan:"]:not(.on), #stage svg.focus [id^="ligand"]:not(.on) { opacity: .2; transition: opacity .12s }
#stage svg [id^="disulfide:"], #stage svg [id^="ss-dot:"], #stage svg [id^="glycan:"], #stage svg [id^="ligand"] { cursor: pointer }
#molpanel { flex: 1.15 }
#mol { flex: 1; min-height: 200px; position: relative; border-radius: 6px; overflow: hidden; background: #fff }
#mappanel { flex: 1 }
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


def _model_cif(path, assembly: str, axes=None, centre=None) -> str:
    """The model as mmCIF; with axes (rows: page-right, page-up, toward viewer) turned into page coordinates, so
    the 3D viewer's default camera sees it the way the figure shows it."""
    model, _ = read_model(path, assembly)
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


def _links(layout, loops, sses, bb, look) -> dict:
    """Disulfides, glycans and ligands as drawn, each with its protein residues (Backbone indices) and, for
    glycans and ligands, their own residues (chain, number) for the 3D view."""
    from .render import ligand_marks

    links = layout.links
    out = {"disulfides": [], "glycans": [], "ligands": []}
    if links is None:
        return out
    shown = lambda r: not layout.partial or layout.res_chain[r] not in layout.partial  # noqa: E731
    label = lambda r: f"{bb.labels[r].chain} {bb.labels[r].name.title()}{bb.labels[r].seq}"  # noqa: E731
    if look.disulfides:
        for i, j in links.disulfides:
            if shown(i) and shown(j):
                out["disulfides"].append(
                    {"key": f"{i}-{j}", "residues": [i, j], "text": f"disulfide {label(i)} – {label(j)}"}
                )
    if look.glycans:
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
    for m in ligand_marks(layout, loops, sses, look):
        g = m["ligand"]
        key = f"{g.name}:{g.chain}{g.seq}"
        held = ", ".join(label(r) for r in g.contacts[:6]) + (" …" if len(g.contacts) > 6 else "")
        out["ligands"].append(
            {
                "key": key,
                "residues": list(g.contacts),
                "het": [list(h) for h in g.members or [(g.chain, g.seq)]],
                "text": f"{'metal' if g.metal else 'ligand'} {g.symbol}, held by {held}",
            }
        )
    return out


def _viewer_head(viewer: str) -> str:
    if viewer == "3dmol":
        return f'<script src="{THREEDMOL}"></script>'
    return f'<link rel="stylesheet" href="{MOLSTAR_CSS}">\n<script src="{MOLSTAR_JS}"></script>'


def build_page(
    path,
    mode: str = "projected",
    look: Style | None = None,
    title: str | None = None,
    viewer: str = "molstar",
    **layout_options,
) -> str:
    """The HTML page as a string. viewer: 'molstar' (Mol*, PDBe build) or '3dmol' for the 3D panel.
    layout_options go to cli.make_layout (symmetry, assembly, swap, ...)."""
    from .cli import make_layout

    if viewer not in VIEWERS:
        raise ValueError(f"unknown 3D viewer {viewer!r}; choose from {', '.join(VIEWERS)}")

    look = look or Style()
    layout, sses, bb = make_layout(path, mode, look=look, **layout_options)
    loops = route_loops(layout, sses, bb)
    svg = save_svg(draw(layout, loops, sses, title, look=look))
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
        "model": _model_cif(path, layout_options.get("assembly", "auto"), layout.page_axes, layout.page_centre),
        "viewer": viewer,
        "links": _links(layout, loops, sses, bb, look),
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
{_viewer_head(viewer)}
</head><body>
<header><h1>{name}</h1><p>{chains} chain(s) · {n_res} residues · {n_el} helices and strands · hover anything: the
same residues light up in all three views · click an element to fly the 3D view to it</p></header>
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
<footer id="credit">{CREDIT}</footer>
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
