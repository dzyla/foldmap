(function () {
  const data = JSON.parse(document.getElementById('site-data').textContent);
  const $ = id => document.getElementById(id);
  const el = (tag, props = {}, kids = []) => {
    const n = Object.assign(document.createElement(tag), props);
    kids.forEach(k => n.append(k));
    return n;
  };
  const SWATCH = {  // a colour per theme for the list
    publication: '#0072b2', journal: '#332288', print: '#8c8c8c', minimal: '#ffffff', presentation: '#e69f00',
    cartoon: '#ee6677', richardson: '#c8414b', rainbow: 'linear-gradient(90deg,#30123b,#28bceb,#a4fc3c,#fb8022,#7a0403)',
    shaded: 'linear-gradient(90deg,#9fc9ef,#0b4f86)', trace: 'linear-gradient(90deg,#cfe4f6,#0072b2)',
    flexibility: 'linear-gradient(90deg,#4575b4,#ffffbf,#d73027)', hydropathy: 'linear-gradient(90deg,#01665e,#f5f5f5,#8c510a)',
    goodsell: '#a6cee3', alphafold: 'linear-gradient(90deg,#ff7d45,#ffdb13,#65cbf3,#0053d6)',
    conservation: 'linear-gradient(90deg,#10c8d1,#ffffff,#a02560)', blueprint: '#0d2a4a',
  };

  // ---- light / dark
  const root = document.documentElement;
  try { const saved = localStorage.getItem('foldmap-theme'); if (saved) root.dataset.theme = saved; } catch (e) { /* private mode */ }
  $('mode-toggle').onclick = () => {
    const dark = root.dataset.theme ? root.dataset.theme === 'dark' : matchMedia('(prefers-color-scheme: dark)').matches;
    root.dataset.theme = dark ? 'light' : 'dark';
    try { localStorage.setItem('foldmap-theme', root.dataset.theme); } catch (e) { /* private mode */ }
  };

  // ---- copy buttons
  document.addEventListener('click', ev => {
    const b = ev.target.closest('.copy');
    if (!b) return;
    const text = (b.dataset.copy ? $(b.dataset.copy) : b.previousElementSibling).textContent;
    navigator.clipboard && navigator.clipboard.writeText(text).then(() => {
      b.textContent = 'copied'; b.classList.add('done');
      setTimeout(() => { b.textContent = 'copy'; b.classList.remove('done'); }, 1200);
    });
  });

  // ---- theme playground
  const state = { structure: data.structures[0] && data.structures[0].id, theme: 'publication' };
  if (!data.themes.some(t => t.name === state.theme) && data.themes[0]) state.theme = data.themes[0].name;
  const chips = $('structure-chips'), list = $('theme-list');
  data.structures.forEach(s => {
    const c = el('button', { className: 'chip', role: 'tab' }, [s.title, el('small', { textContent: s.id })]);
    c.title = s.note;
    c.onclick = () => { state.structure = s.id; showTheme(); };
    c.dataset.id = s.id;
    chips.append(c);
  });
  data.themes.forEach(t => {
    const sw = el('span', { className: 'sw' });
    sw.style.background = SWATCH[t.name] || '#999';
    if (t.name === 'minimal') sw.style.border = '1px solid #999';
    const b = el('button', { className: 'theme-item', role: 'option' }, [sw, t.name]);
    b.dataset.name = t.name;
    b.onclick = () => { state.theme = t.name; showTheme(); };
    list.append(b);
  });
  function showTheme() {
    const s = data.structures.find(x => x.id === state.structure);
    if (!s) return;
    let theme = state.theme;
    if (!s.figures[theme]) theme = 'publication';  // e.g. conservation is only drawn for ubiquitin
    const t = data.themes.find(x => x.name === theme);
    chips.querySelectorAll('.chip').forEach(c => c.setAttribute('aria-selected', c.dataset.id === s.id));
    list.querySelectorAll('.theme-item').forEach(b => {
      b.setAttribute('aria-selected', b.dataset.name === theme);
      b.disabled = !s.figures[b.dataset.name];
      b.style.opacity = b.disabled ? .4 : 1;
    });
    const panel = $('theme-panel'), img = $('theme-img');
    panel.classList.add('loading');
    panel.classList.toggle('dark', !!t.dark);
    img.onload = () => panel.classList.remove('loading');
    img.src = s.figures[theme];
    img.alt = `${s.title} (${s.id}) in the ${theme} theme`;
    $('theme-name').textContent = theme;
    $('theme-note').textContent = t.note;
    const extra = theme === 'conservation' ? ' --msa family.fasta' : '';
    $('theme-cmd').textContent = `foldmap plot ${s.id} --theme ${theme}${extra} -o ${s.id}.svg`;
  }
  showTheme();

  // ---- modes
  const tabs = $('mode-tabs'), body = $('mode-body');
  function showMode(key) {
    const m = data.modes.find(x => x.key === key);
    tabs.querySelectorAll('.tab').forEach(t => t.setAttribute('aria-selected', t.dataset.key === key));
    body.replaceChildren(
      el('div', { className: 'mode-text' }, [el('h3', { textContent: m.title }), el('p', { textContent: m.text })]),
      el('div', { className: 'mode-panels' + (m.panels.length > 1 ? ' two' : '') }, m.panels.map(p => {
        const fp = el('div', { className: 'figure-panel' + (p.dark ? ' dark' : '') },
                      [el('img', { src: p.src, alt: p.caption, loading: 'lazy' })]);
        const cmd = el('div', { className: 'cmd' }, [el('code', { textContent: p.command }),
                                                   el('button', { className: 'copy', textContent: 'copy' })]);
        return el('figure', {}, [fp, el('figcaption', { textContent: p.caption }), cmd]);
      })),
    );
  }
  data.modes.forEach(m => {
    const t = el('button', { className: 'tab', role: 'tab', textContent: m.title });
    t.dataset.key = m.key;
    t.onclick = () => showMode(m.key);
    tabs.append(t);
  });
  if (data.modes.length) showMode(data.modes[0].key);

  // ---- sequence views
  const grid = $('sequence-grid');
  data.sequences.forEach(s => grid.append(el('figure', {}, [
    el('div', { className: 'figure-panel' }, [el('img', { src: s.src, alt: s.caption, loading: 'lazy' })]),
    el('figcaption', { textContent: s.caption }),
    el('div', { className: 'cmd' }, [el('code', { textContent: s.command }), el('button', { className: 'copy', textContent: 'copy' })]),
  ])));
  if (!data.sequences.length) $('sequence').hidden = true;

  // ---- explorers: loaded when the section is first seen
  const frame = $('explorer'), echips = $('explorer-chips');
  function showExplorer(id) {
    const x = data.explorers.find(e => e.id === id);
    echips.querySelectorAll('.chip').forEach(c => c.setAttribute('aria-selected', c.dataset.id === id));
    frame.src = x.src;
    $('explorer-open').href = x.src;
  }
  data.explorers.forEach(x => {
    const c = el('button', { className: 'chip' }, [x.title, el('small', { textContent: x.id })]);
    c.dataset.id = x.id;
    c.onclick = () => showExplorer(x.id);
    echips.append(c);
  });
  if (data.explorers.length) {
    new IntersectionObserver((entries, obs) => {
      if (entries.some(e => e.isIntersecting)) { showExplorer(data.explorers[0].id); obs.disconnect(); }
    }, { rootMargin: '300px' }).observe($('explore'));
  }
})();
