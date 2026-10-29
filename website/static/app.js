/* BrainCellReg front end: Vue 3 (global build) + ECharts (SVG renderer). Every figure is drawn here from the API answers. */
const { createApp, ref, reactive, computed, watch, onMounted, onBeforeUnmount, nextTick, h } = Vue;

/* ---------- shared constants (order and colours as in the paper's figures) ---------- */
const REG = ['FCTX', 'TCTX', 'OCTX', 'HIPP', 'PUTM', 'THAL', 'SNIG', 'MEDU', 'WHMT', 'CRBL'];
const CELL = ['ExN', 'InN', 'Astro', 'Olig', 'OPC', 'Micro', 'Endo', 'Perc', 'NPC', 'Gran', 'Purk'];
// same as PAL_REG / PAL_CELL in server_scripts/fig_theme.R
const PAL_REG = { CRBL: '#4E79A7', FCTX: '#F28E2B', HIPP: '#59A14F', MEDU: '#B07AA1', OCTX: '#EDC948', PUTM: '#76B7B2', SNIG: '#E15759', TCTX: '#FF9DA7', THAL: '#9C755F', WHMT: '#8C8C8C' };
const PAL_CELL = { ExN: '#C0504D', InN: '#E8A33D', Astro: '#4F81BD', Olig: '#1B9E77', OPC: '#7FC8A9', Micro: '#8064A2', Endo: '#4BACC6', Perc: '#A6761D', NPC: '#B3B3B3', Gran: '#666666', Purk: '#E78AC3' };
const C = { ink: '#141C2E', mut: '#5A6478', line: '#D9DFE7', acc: '#5B3FA0', up: '#C0392B', dn: '#2E6DA4', na: '#EEF1F5' };
const SEQ = ['#F3F0FA', '#D6CCEE', '#A893D9', '#7A5EC0', '#5B3FA0', '#2E1D5C'];
const DIV = [C.dn, '#9DBCDB', '#F4F4F4', '#E5A39C', C.up];
const LDCOL = [['#357EBD', 0.2], ['#46B8DA', 0.4], ['#5CB85C', 0.6], ['#EEA236', 0.8], ['#D43F3A', 1.01]];
const FONT = "'IBM Plex Sans', 'Segoe UI', sans-serif";
let META = null;

const api = async (p) => {
  const r = await fetch('/api' + p);
  if (!r.ok) { let m = r.statusText; try { m = (await r.json()).detail || m; } catch (e) {} throw new Error(m); }
  return r.json();
};
const fp = (p) => p == null ? '—' : (p === 0 ? '<1e-300' : (p < 0.001 ? p.toExponential(1).replace('e', 'e') : p.toFixed(3)));
const f2 = (x, d = 2) => x == null ? '—' : Number(x).toFixed(d);
const fint = (x) => x == null ? '—' : Number(x).toLocaleString('en-US');
const fsci = (x) => x == null ? '—' : (x !== 0 && Math.abs(x) < 0.001 ? Number(x).toExponential(1) : Number(x).toFixed(3));  // cell fractions span 1e-10 to 0.5
const nlp = (p) => p == null ? null : (p <= 0 ? 300 : -Math.log10(p));
const go = (path) => { location.hash = '#' + path; };
const ctxLabel = (ctx) => ctx.replace('.', ' · ');

/* ---------- chart wrapper: SVG renderer, download as SVG or PNG ---------- */
const Chart = {
  props: { option: Object, height: { type: Number, default: 360 }, name: { type: String, default: 'BrainCellReg' } },
  emits: ['pick'],
  setup(props, { emit }) {
    const el = ref(null); let ch = null; let ro = null;
    const draw = () => { if (ch && props.option) ch.setOption(props.option, true); };
    onMounted(() => {
      ch = echarts.init(el.value, null, { renderer: 'svg' });
      ch.on('click', (e) => emit('pick', e));
      draw();
      ro = new ResizeObserver(() => ch && ch.resize()); ro.observe(el.value);
    });
    onBeforeUnmount(() => { ro && ro.disconnect(); ch && ch.dispose(); });
    watch(() => props.option, draw, { deep: false });
    const save = (type) => {
      const svg = ch.renderToSVGString();
      const fn = props.name.replace(/[^\w.-]+/g, '_');
      if (type === 'svg') {
        const a = document.createElement('a'); a.href = URL.createObjectURL(new Blob([svg], { type: 'image/svg+xml' })); a.download = fn + '.svg'; a.click(); return;
      }
      const img = new Image(); const s = 3;
      img.onload = () => {
        const cv = document.createElement('canvas'); cv.width = img.width * s; cv.height = img.height * s;
        const cx = cv.getContext('2d'); cx.fillStyle = '#fff'; cx.fillRect(0, 0, cv.width, cv.height); cx.scale(s, s); cx.drawImage(img, 0, 0);
        const a = document.createElement('a'); a.href = cv.toDataURL('image/png'); a.download = fn + '.png'; a.click();
      };
      img.src = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(svg);
    };
    return () => h('div', { class: 'fig' }, [
      h('div', { class: 'dl' }, [h('button', { class: 'btn', onClick: () => save('svg') }, 'SVG'), h('button', { class: 'btn', onClick: () => save('png') }, 'PNG')]),
      h('div', { ref: el, style: { height: props.height + 'px', width: '100%' } }),
    ]);
  },
};

const baseText = { fontFamily: FONT, color: C.ink, fontSize: 12 };

/* The 110-context grid: rows = regions (cortex at the top, cerebellum at the bottom), columns = cell types.
   values: {ctx: number|null}; text: {ctx: tooltip lines}; scale: 'seq' | 'div'; marks: Set of ctx outlined. */
function gridOption({ values, scale = 'seq', min, max, title = '', tip, marks = new Set(), selected = null, legend = '', compact = window.innerWidth < 700 }) {
  const data = [];
  REG.forEach((r, yi) => CELL.forEach((c, xi) => {
    const k = r + '.' + c; const v = values[k];
    data.push({ value: [xi, yi, v == null ? '-' : v], ctx: k,
      itemStyle: { borderColor: k === selected ? C.ink : (marks.has(k) ? C.ink : '#fff'), borderWidth: k === selected ? 2.5 : (marks.has(k) ? 1.5 : 1) } });
  }));
  const vals = data.map(d => d.value[2]).filter(v => v !== '-');
  const lo = min != null ? min : Math.min(...vals, 0), hi = max != null ? max : Math.max(...vals, 1);
  return {
    textStyle: baseText, animation: false,
    tooltip: { formatter: (p) => tip ? tip(p.data.ctx, p.data.value[2]) : `${ctxLabel(p.data.ctx)}<br>${p.data.value[2]}`, textStyle: { fontFamily: FONT } },
    grid: { left: compact ? 52 : 165, right: 16, top: 26, bottom: 54 },
    xAxis: { type: 'category', data: CELL, position: 'top', axisTick: { show: false }, axisLine: { show: false },
      axisLabel: { color: C.ink, fontSize: 12, fontWeight: 600, interval: 0, fontFamily: FONT } },
    yAxis: { type: 'category', data: REG, inverse: true, axisTick: { show: false }, axisLine: { show: false },
      axisLabel: { formatter: (r) => compact ? `{s${r}|■} ${r}` : `{s${r}|■} ${META ? META.regions.find(x => x.id === r).name : r}`,
        rich: Object.fromEntries(REG.map(r => ['s' + r, { color: PAL_REG[r], fontSize: 14 }])), color: C.ink, fontSize: 12 } },
    visualMap: { seriesIndex: 1, min: lo, max: hi, calculable: false, orient: 'horizontal', left: compact ? 52 : 165, bottom: 4, itemHeight: 160, itemWidth: 10,
      text: [Number.isInteger(hi) ? fint(hi) : f2(hi, 1), (Number.isInteger(lo) ? fint(lo) : f2(lo, 1)) + (legend ? '   ' + legend : '')], textStyle: { fontFamily: FONT, color: C.mut, fontSize: 11 },
      inRange: { color: scale === 'div' ? DIV : SEQ }, outOfRange: { color: C.na } },
    series: [
      // background: every context, so contexts without a value read as "nothing here" rather than as missing cells
      { type: 'heatmap', silent: true, data: data.map(x => ({ value: [x.value[0], x.value[1], 0], itemStyle: { color: C.na, borderColor: '#fff', borderWidth: 1 } })) },
      { type: 'heatmap', data: data.filter(x => x.value[2] !== '-'), label: { show: false }, emphasis: { itemStyle: { borderColor: C.ink, borderWidth: 2 } } },
      // outline of marked / selected contexts drawn on top
      { type: 'heatmap', silent: true, data: data.filter(x => x.itemStyle.borderWidth > 1).map(x => ({ value: [x.value[0], x.value[1], 0],
        itemStyle: { color: 'transparent', borderColor: C.ink, borderWidth: x.itemStyle.borderWidth } })) }],
  };
}

// note: a line written inside the figure (e.g. β and P of a genotype plot), so it is kept when the figure is saved
function boxOption(groups, { xname = '', yname = 'Expression', colorOf = () => C.acc, labelOf = (g) => g.group, note = '' } = {}) {
  return {
    textStyle: baseText, animation: false,
    title: note ? { text: note, left: 'center', top: 0, textStyle: { fontSize: 12, fontWeight: 600, color: C.ink, fontFamily: FONT } } : undefined,
    tooltip: { trigger: 'item', textStyle: { fontFamily: FONT },
      formatter: (p) => { const g = groups[p.dataIndex]; return `${labelOf(g)}<br>${g.n != null ? `n = ${g.n}<br>` : ''}median ${f2(g.median)}<br>IQR ${f2(g.q1)}–${f2(g.q3)}`; } },
    grid: { left: 52, right: 12, top: note ? 34 : 18, bottom: 46 },
    xAxis: { type: 'category', data: groups.map(g => labelOf(g)), name: xname, nameLocation: 'middle', nameGap: 30,
      axisLabel: { color: C.ink, interval: 0, fontSize: 12 }, axisLine: { lineStyle: { color: C.mut } }, axisTick: { show: false } },
    yAxis: { type: 'value', name: yname, nameLocation: 'middle', nameGap: 38, scale: true, splitLine: { lineStyle: { color: '#EEF1F5' } }, axisLabel: { color: C.mut } },
    series: [{ type: 'boxplot', boxWidth: [14, 46],
      data: groups.map((g, i) => ({ value: [g.lo ?? g.q1, g.q1, g.median, g.q3, g.hi ?? g.q3],
        itemStyle: { color: colorOf(g, i) + '33', borderColor: colorOf(g, i), borderWidth: 1.4 } })) },
      { type: 'scatter', symbolSize: 4, itemStyle: { color: C.mut, opacity: .6 }, data: groups.flatMap((g, i) => (g.out || []).map(v => [i, v])) }],
  };
}

/* ---------- search box ---------- */
const Search = {
  props: { big: Boolean, placeholder: { type: String, default: 'Gene, rs number or disease' } },
  setup(props) {
    const qtext = ref(''); const items = ref([]); const hi = ref(0); const open = ref(false); let t = null;
    watch(qtext, (v) => { clearTimeout(t); t = setTimeout(async () => { items.value = v.trim() ? await api('/search?q=' + encodeURIComponent(v)) : []; hi.value = 0; open.value = true; }, 120); });
    const pick = (it) => { open.value = false; qtext.value = ''; go(`/${it.type}/${it.type === 'gene' ? it.label : it.id}`); };
    const key = (e) => {
      if (e.key === 'ArrowDown') { hi.value = Math.min(hi.value + 1, items.value.length - 1); e.preventDefault(); }
      else if (e.key === 'ArrowUp') { hi.value = Math.max(hi.value - 1, 0); e.preventDefault(); }
      else if (e.key === 'Enter' && items.value[hi.value]) pick(items.value[hi.value]);
      else if (e.key === 'Escape') open.value = false;
    };
    const blur = () => setTimeout(() => { open.value = false; }, 150);
    return { qtext, items, hi, open, pick, key, blur };
  },
  template: `<div class="search">
    <input v-model="qtext" @keydown="key" @focus="open=true" @blur="blur" :placeholder="placeholder" aria-label="Search genes, variants and diseases">
    <ul class="sugg" v-if="open && items.length">
      <li v-for="(it,i) in items" :key="it.type+it.id" :class="{hi:i===hi}" @mousedown.prevent="pick(it)">
        <span class="t">{{it.type}}</span><span>{{it.label}}</span><span class="s" v-if="it.sub">{{it.sub}}</span></li>
    </ul></div>`,
};

/* ---------- searchable drop-down (PI 10-19: long lists such as 277 eSNPs need a search box) ----------
   options: [{ value, label, sub }]; v-model is the value. Typing filters by label / sub (case-insensitive); ↑ ↓ Enter Esc. */
const Combo = {
  props: { modelValue: { default: null }, options: { type: Array, default: () => [] }, label: { type: String, default: '' }, width: { type: String, default: '260px' } },
  emits: ['update:modelValue'],
  setup(props, { emit }) {
    const open = ref(false); const text = ref(''); const hi = ref(0); const box = ref(null); const inp = ref(null);
    const cur = computed(() => props.options.find(o => o.value === props.modelValue));
    const view = computed(() => {
      const s = text.value.trim().toLowerCase();
      const r = s ? props.options.filter(o => (o.label + ' ' + (o.sub || '')).toLowerCase().includes(s)) : props.options;
      return r.slice(0, 300);
    });
    const nMatch = computed(() => { const s = text.value.trim().toLowerCase(); return s ? props.options.filter(o => (o.label + ' ' + (o.sub || '')).toLowerCase().includes(s)).length : props.options.length; });
    const show = async () => { open.value = true; text.value = ''; hi.value = Math.max(0, view.value.findIndex(o => o.value === props.modelValue)); await nextTick(); inp.value && inp.value.focus(); };
    const pick = (o) => { emit('update:modelValue', o.value); open.value = false; };
    const key = (e) => {
      if (e.key === 'ArrowDown') { hi.value = Math.min(hi.value + 1, view.value.length - 1); e.preventDefault(); }
      else if (e.key === 'ArrowUp') { hi.value = Math.max(hi.value - 1, 0); e.preventDefault(); }
      else if (e.key === 'Enter' && view.value[hi.value]) pick(view.value[hi.value]);
      else if (e.key === 'Escape') open.value = false;
    };
    watch(text, () => { hi.value = 0; });
    const away = (e) => { if (open.value && box.value && !box.value.contains(e.target)) open.value = false; };
    onMounted(() => document.addEventListener('mousedown', away));
    onBeforeUnmount(() => document.removeEventListener('mousedown', away));
    return { open, text, hi, box, inp, cur, view, nMatch, show, pick, key };
  },
  template: `<div class="combo" ref="box" :style="{width}">
    <button type="button" class="combo-btn" @click="open ? (open=false) : show()" :aria-label="label">
      <span>{{cur ? cur.label : '—'}}</span><span class="s" v-if="cur && cur.sub">{{cur.sub}}</span><span class="caret">▾</span></button>
    <div class="combo-pop" v-if="open">
      <input ref="inp" v-model="text" @keydown="key" :placeholder="'Type to search ' + options.length.toLocaleString() + ' items'" :aria-label="'Search ' + label">
      <ul><li v-for="(o,i) in view" :key="String(o.value)" :class="{hi:i===hi, on:o.value===modelValue}" @mousedown.prevent="pick(o)">
        <span>{{o.label}}</span><span class="s" v-if="o.sub">{{o.sub}}</span></li>
        <li class="none" v-if="!view.length">No match</li></ul>
      <p class="combo-n" v-if="nMatch > view.length">Showing {{view.length}} of {{nMatch.toLocaleString()}}; type more to narrow.</p>
    </div></div>`,
};

/* ---------- sortable table ---------- */
const DataTable = {
  props: { rows: Array, cols: Array, filter: { type: Boolean, default: false }, pick: { type: Function, default: null }, picked: { default: null }, max: { type: Number, default: 500 } },
  setup(props) {
    const key = ref(null); const asc = ref(true); const f = ref('');
    const view = computed(() => {
      let r = props.rows || [];
      if (f.value) { const s = f.value.toLowerCase(); r = r.filter(x => props.cols.some(c => String(x[c.k] ?? '').toLowerCase().includes(s))); }
      if (key.value) { const k = key.value; r = [...r].sort((a, b) => { const x = a[k], y = b[k]; if (x == null) return 1; if (y == null) return -1; return (x > y ? 1 : x < y ? -1 : 0) * (asc.value ? 1 : -1); }); }
      return r;
    });
    const sort = (k) => { if (key.value === k) asc.value = !asc.value; else { key.value = k; asc.value = true; } };
    return { key, asc, f, view, sort };
  },
  template: `<div>
    <input v-if="filter" class="filter" v-model="f" placeholder="Filter rows" aria-label="Filter rows">
    <div class="tw"><table><thead><tr><th v-for="c in cols" :class="{n:c.n}" @click="sort(c.k)">{{c.t}}<span v-if="key===c.k">{{asc?' ↑':' ↓'}}</span></th></tr></thead>
    <tbody><tr v-for="(r,i) in view.slice(0,max)" :key="i" :class="{pick: picked!=null && picked(r)}" @click="pick && pick(r)" :style="pick?'cursor:pointer':''">
      <td v-for="c in cols" :class="{n:c.n}"><span v-if="c.html" v-html="c.f(r)"></span><span v-else>{{ c.f ? c.f(r) : r[c.k] }}</span></td></tr></tbody></table></div>
    <p class="cap" v-if="view.length>max">Showing {{max}} of {{view.length.toLocaleString()}} rows; filter or download for the rest.</p></div>`,
};

const ctxOptions = () => REG.flatMap(r => CELL.map(c => r + '.' + c));
const chip = (ctx) => `<span class="ctxchip"><span class="sw" style="background:${PAL_REG[ctx.split('.')[0]]}"></span>${ctxLabel(ctx)}</span>`;

/* ---------- pages ---------- */
const Home = {
  components: { Search, Chart },
  setup() {
    const m = META;
    const values = Object.fromEntries(m.contexts.map(c => [c.ctx, c.egenes]));
    const opt = gridOption({ values, max: Math.max(...m.contexts.map(c => c.egenes)), legend: 'eGenes',
      tip: (k, v) => { const c = m.contexts.find(x => x.ctx === k); return `<b>${ctxLabel(k)}</b><br>${v} eGenes · ${fint(c.independent_esnps)} independent eSNPs<br>${c.n_samples} samples, mean cell fraction ${f2(c.mean_frac, 3)}`; } });
    const pick = (e) => e.data && go('/context/' + e.data.ctx);
    // estimated cell composition of the 110 contexts; a cell opens that cell type in all regions (PI 10-20「细胞比例就放首页，也像这个图可以点进去」)
    const fracOpt = gridOption({ values: Object.fromEntries(m.contexts.map(c => [c.ctx, c.mean_frac])), min: 0, max: Math.max(...m.contexts.map(c => c.mean_frac)), legend: 'mean cell fraction',
      tip: (k) => { const c = m.contexts.find(x => x.ctx === k); return `<b>${ctxLabel(k)}</b><br>mean fraction ${fsci(c.mean_frac)} · median ${fsci(c.median_frac)}<br>${c.n_samples} samples${c.median_frac < 0.001 ? '<br>almost no cells of this type estimated here' : ''}`; } });
    const pickFrac = (e) => { if (e.data) { const [r, c] = e.data.ctx.split('.'); go(`/celltype/${c}:${r}`); } };
    return { m, opt, pick, fracOpt, pickFrac, fint };
  },
  template: `<div class="home">
    <section>
      <h1>Gene regulation across ten brain regions and eleven cell types</h1>
      <p class="lede">Expression quantitative trait loci from 1,231 post-mortem brain samples of {{m.totals.individuals}} donors, split by estimated cell composition into 110 region × cell-type contexts and tested against ten neurological and psychiatric diseases.</p>
      <Search big placeholder="Search a gene, an rs number or a disease"/>
      <p class="try">Try <a href="#/gene/GPNMB">GPNMB</a><a href="#/variant/rs156429">rs156429</a><a href="#/disease/PD">Parkinson's disease</a></p>
      <ul class="facts">
        <li><span>Genes tested</span><span>{{fint(m.totals.genes_tested)}}</span></li>
        <li><span>eGenes (q &lt; 0.05 in at least one context)</span><span>{{fint(m.totals.egenes)}}</span></li>
        <li><span>Independent eSNPs</span><span>{{fint(m.totals.independent_esnps)}}</span></li>
        <li><span>eSNP–gene pairs</span><span>{{fint(m.totals.esnp_pairs)}}</span></li>
        <li><span>Contexts with at least one eGene</span><span>{{m.totals.contexts_with_egenes}} of {{m.totals.contexts}}</span></li>
      </ul>
      <h2 style="margin-top:26px">Diseases</h2>
      <div class="dz-row"><a v-for="d in m.diseases" :href="'#/disease/'+d.id">{{d.name}}</a></div>
    </section>
    <div>
    <section class="card">
      <div class="cg-h"><h2>eGenes in each context</h2><span class="cg-note">Select a cell to open the context</span></div>
      <Chart :option="opt" :height="430" name="BrainCellReg_eGenes_per_context" @pick="pick"/>
    </section>
    <section class="card">
      <div class="cg-h"><h2>Cell fractions in each context</h2><span class="cg-note">Select a cell to compare that cell type across regions</span></div>
      <Chart :option="fracOpt" :height="430" name="BrainCellReg_cell_fractions" @pick="pickFrac"/>
    </section></div></div>`,
};

/* one cell type's estimated fraction in the 10 regions; id = "Purk" or "Purk:CRBL" (region to highlight) */
const CellPage = {
  components: { Chart },
  props: ['id'],
  setup(props) {
    const d = ref(null); const err = ref(''); const scale = ref('lin');
    const [cell, hl] = props.id.split(':');
    (async () => { try { d.value = await api('/celltype/' + cell); } catch (e) { err.value = e.message; } })();
    const tr = (v) => scale.value === 'log' ? Math.log10(Math.max(v, 1e-12)) : v;
    const opt = computed(() => {
      if (!d.value) return null;
      const g = d.value.regions.map(x => ({ ...x, group: x.region, q1: tr(x.q1), median: tr(x.median), q3: tr(x.q3), lo: tr(x.lo), hi: tr(x.hi) }));
      const o = boxOption(g, { yname: scale.value === 'log' ? 'log10 estimated fraction' : 'Estimated fraction', colorOf: (x) => PAL_REG[x.region] });
      o.series[0].data.forEach((p, i) => { p.ctx = g[i].ctx; if (g[i].region === hl) p.itemStyle.borderWidth = 3; });
      o.grid.left = 78;
      o.yAxis.nameGap = 60;
      o.yAxis.axisLabel.formatter = (v) => scale.value === 'log' ? String(v) : (v !== 0 && Math.abs(v) < 0.001 ? v.toExponential(0) : String(v));
      if (scale.value === 'log') o.yAxis.max = 0;
      o.tooltip.formatter = (p) => { const x = d.value.regions[p.dataIndex]; return x ? `<b>${x.region_name}</b><br>n = ${x.n}<br>mean ${fsci(x.mean)} · median ${fsci(x.median)}<br>IQR ${fsci(x.q1)}–${fsci(x.q3)}` : ''; };
      return o;
    });
    const top = computed(() => d.value ? d.value.regions.slice().sort((a, b) => b.median - a.median)[0] : null);
    const pick = (e) => e.data && e.data.ctx && go('/context/' + e.data.ctx);
    return { d, err, scale, opt, top, pick, hl, cell, fsci, fint, PAL_REG, PAL_CELL };
  },
  template: `<div>
    <div v-if="err" class="err">{{err}}</div>
    <div v-else-if="!d" class="load">Loading…</div>
    <template v-else>
      <div class="ph"><span class="kind">Cell type</span><h1><span class="sw" :style="{background:PAL_CELL[cell],width:'18px',height:'18px',marginRight:'10px'}"></span>{{d.cell_name}}</h1></div>
      <p class="sum">Estimated fraction of this cell type in the samples of each region. Highest: <b>{{top.region_name}}</b> (median {{fsci(top.median)}}).</p>
      <section class="card"><div class="cg-h"><h2>Fraction by region</h2>
        <div class="sel" style="margin:0"><button class="btn" :class="{on:scale==='lin'}" @click="scale='lin'">Linear</button><button class="btn" :class="{on:scale==='log'}" @click="scale='log'">Log scale</button></div></div>
        <p class="cap">Fractions estimated by deconvolution (BayesPrism / InstaPrism with the STAB reference). Select a box to open that context.</p>
        <Chart :option="opt" :height="360" :name="'BrainCellReg_'+cell+'_fraction_by_region'" @pick="pick"/></section>
      <section class="card"><h2>Contexts</h2>
        <div class="tw" style="max-height:none"><table><thead><tr><th>Context</th><th class="n">Samples</th><th class="n">Mean</th><th class="n">Median</th><th class="n">IQR</th><th class="n">eGenes</th><th></th></tr></thead>
        <tbody><tr v-for="x in d.regions" :class="{on:x.region===hl}"><td><a :href="'#/context/'+x.ctx"><span class="sw" :style="{background:PAL_REG[x.region]}"></span>{{x.region_name}}</a></td>
          <td class="n">{{x.n}}</td><td class="n">{{fsci(x.mean)}}</td><td class="n">{{fsci(x.median)}}</td><td class="n">{{fsci(x.q1)}}–{{fsci(x.q3)}}</td><td class="n">{{fint(x.egenes)}}</td>
          <td><span v-if="x.median<0.001" class="tag">near zero</span></td></tr></tbody></table></div>
        <p class="cap" v-if="d.regions.some(x=>x.median<0.001)">"Near zero": median fraction below 0.001, i.e. the deconvolution assigns almost no cells of this type to the region; eQTL and disease results of such a context should be read with caution.</p></section>
    </template></div>`,
};

const GenePage = {
  components: { Chart, DataTable, Combo },
  props: ['id'],
  setup(props) {
    const d = ref(null); const err = ref(''); const ctx = ref(null); const groups = ref(null); const geno = ref(null); const genoRs = ref(''); const esnps = ref([]);
    const loc = ref(null); const locDz = ref('PD'); const locCtx = ref(null); const locErr = ref('');
    const load = async () => {
      d.value = null; err.value = '';
      try {
        d.value = await api('/gene/' + props.id);
        const best = d.value.eqtl.filter(x => x.is_egene).sort((a, b) => a.qval - b.qval)[0] || d.value.eqtl.sort((a, b) => a.pval_nominal - b.pval_nominal)[0];
        ctx.value = best ? best.ctx : 'FCTX.ExN';
        locCtx.value = ctx.value;
        const ev = d.value.evidence.slice().sort((a, b) => b.n_ctx_smr_heidi_h4 - a.n_ctx_smr_heidi_h4)[0];
        if (ev) locDz.value = ev.disease;
      } catch (e) { err.value = e.message; }
    };
    const eq = computed(() => d.value ? Object.fromEntries(d.value.eqtl.map(x => [x.ctx, x])) : {});
    const lead = computed(() => eq.value[ctx.value]);
    const gridEq = computed(() => {
      if (!d.value) return null;
      const vals = {}; const marks = new Set();
      d.value.eqtl.forEach(x => { if (x.is_egene) { vals[x.ctx] = Math.sign(x.slope) * Math.min(nlp(x.pval_nominal), 30); marks.add(x.ctx); } });
      const m = Math.max(5, ...Object.values(vals).map(Math.abs));
      return gridOption({ values: vals, scale: 'div', min: -m, max: m, selected: ctx.value, legend: 'signed −log10 P, eGene contexts',
        tip: (k) => { const x = eq.value[k]; if (!x) return `${ctxLabel(k)}<br>not tested`; return `<b>${ctxLabel(k)}</b><br>${x.is_egene ? 'eGene' : 'not an eGene'} · q ${fp(x.qval)}<br>lead ${x.lead_rsid}${x.effect_allele ? ` (effect allele ${x.effect_allele})` : ''}, β ${f2(x.slope)}, P ${fp(x.pval_nominal)}`; } });
    });
    const gridEx = computed(() => {
      if (!d.value || !d.value.expr.length) return null;
      const vals = Object.fromEntries(d.value.expr.map(x => [x.ctx, x.median]));
      return gridOption({ values: vals, selected: ctx.value, legend: 'median expression',
        tip: (k) => { const x = d.value.expr.find(e => e.ctx === k); return x ? `<b>${ctxLabel(k)}</b><br>median ${f2(x.median)} (IQR ${f2(x.q1)}–${f2(x.q3)}), n ${x.n}` : ctxLabel(k); } });
    });
    const loadCtx = async () => {
      groups.value = null; geno.value = null;
      if (!ctx.value) return;
      const x = eq.value[ctx.value];
      // every eSNP of the gene in this context can be chosen (PI 10-19「这里怎么选都是只有一个snp」); the lead eSNP first
      esnps.value = [];
      // both requests at once rather than one after the other (#146 speed-up)
      const [gr, es] = await Promise.all([
        api(`/gene/${props.id}/groups?ctx=${ctx.value}`).catch(e => ({ error: e.message })),
        x && x.is_egene ? api(`/gene/${props.id}/esnps?ctx=${ctx.value}`).then(r => r.rows).catch(() => []) : Promise.resolve([])]);
      groups.value = gr; esnps.value = es;
      const old = genoRs.value;
      genoRs.value = x && x.is_egene ? x.lead_rsid : '';
      if (old === genoRs.value) loadGeno();   // same lead in the new context: the watcher does not fire
    };
    const loadGeno = async () => {
      geno.value = null;
      if (genoRs.value) { try { geno.value = await api(`/genotype?gene=${props.id}&rsid=${genoRs.value}`); } catch (e) { geno.value = { error: e.message }; } }
    };
    watch(genoRs, loadGeno);
    const genoStat = computed(() => esnps.value.find(r => r.rsid === genoRs.value) || null);
    const ctxOpts = ctxOptions().map(c => ({ value: c, label: ctxLabel(c) }));
    const esnpOpts = computed(() => esnps.value.map(r => ({ value: r.rsid, sub: 'P ' + fp(r.pval),
      label: r.rsid + (lead.value && r.rsid === lead.value.lead_rsid ? ' (lead)' : (r.independent ? ' (independent)' : '')) })));
    const loadLocus = async () => {
      loc.value = null; locErr.value = '';
      if (!d.value || !d.value.locus) return;
      try { loc.value = await api(`/locus/${props.id}?ctx=${locCtx.value}&disease=${locDz.value}`); } catch (e) { locErr.value = e.message; }
    };
    watch(() => props.id, load, { immediate: true });
    watch(ctx, loadCtx);
    watch([locCtx, locDz, () => d.value], loadLocus);
    const pickGrid = (e) => { if (e.data) { ctx.value = e.data.ctx; locCtx.value = e.data.ctx; } };
    const sexOpt = computed(() => groups.value && groups.value.sex ? boxOption(groups.value.sex, { colorOf: (g) => g.group === 'Female' ? '#B5577E' : '#3E6E9E', labelOf: (g) => `${g.group} (${g.n})` }) : null);
    const AGE_ORDER = ['children and adolescents', 'young adults', 'adults', 'middle-aged', 'elderly'];
    const ageOpt = computed(() => {
      if (!groups.value || !groups.value.age) return null;
      const a = groups.value.age.slice().sort((x, y) => (AGE_ORDER.indexOf(x.group) + 99) % 99 - (AGE_ORDER.indexOf(y.group) + 99) % 99);
      const o = boxOption(a, { colorOf: () => '#6C7A93', labelOf: (g) => `${g.group}\n(${g.n})` });
      o.xAxis.axisLabel = { ...o.xAxis.axisLabel, width: 78, overflow: 'break', lineHeight: 14 };
      o.grid.bottom = 58;
      return o;
    });
    const genoOpt = computed(() => {
      if (!geno.value || !geno.value.boxes) return null;
      const b = geno.value.boxes.filter(x => x.ctx === ctx.value).sort((a, c) => a.n_effect - c.n_effect);
      if (!b.length) return null;
      const ea = b[0].effect_allele, oa = b[0].other_allele;
      const lab = (g) => ea && oa ? [oa + oa, oa + ea, ea + ea][g.n_effect] + ` (${g.n})` : `${g.n_effect} copies (${g.n})`;
      const st = genoStat.value || (geno.value.stats || []).find(s => s.ctx === ctx.value);
      const note = st ? `${ctxLabel(ctx.value)}: β = ${f2(st.slope)}${ea ? ` per copy of ${ea}` : ''}, P = ${fp(st.pval)}` : '';
      return boxOption(b.map(x => ({ ...x, lo: x.whisker_lo, hi: x.whisker_hi })), { colorOf: (g) => ['#9DBCDB', '#7A5EC0', '#2E1D5C'][g.n_effect], labelOf: lab, xname: genoRs.value, note });
    });
    const locOpt = computed(() => loc.value && loc.value.eqtl.length ? locusOption(loc.value) : null);
    const evCols = [
      { t: 'Disease', k: 'disease_name', f: r => r.disease_name, html: false },
      { t: 'Contexts with SMR q<0.05', k: 'n_ctx_smr', n: true },
      { t: 'HEIDI P≥0.01', k: 'n_ctx_heidi_ge_0_01', n: true },
      { t: 'Coloc PP.H4≥0.8', k: 'n_ctx_h4_ge_0_8', n: true },
      { t: 'All three', k: 'n_ctx_smr_heidi_h4', n: true, html: true, f: r => r.n_ctx_smr_heidi_h4 > 0 ? `<span class="yes">${r.n_ctx_smr_heidi_h4}</span>` : '<span class="no">0</span>' },
      { t: 'Best q (SMR)', k: 'min_q_smr', n: true, f: r => fp(r.min_q_smr) },
      { t: 'Max PP.H4', k: 'max_h4', n: true, f: r => f2(r.max_h4) },
      { t: 'scGWAS module', k: 'scgwas_module', html: true, f: r => r.scgwas_module ? `<span class="yes">yes</span> ${r.scgwas_ctx || ''}` : '<span class="no">no</span>' },
    ];
    const eqCols = [
      { t: 'Context', k: 'ctx', html: true, f: r => chip(r.ctx) },
      { t: 'eGene', k: 'is_egene', html: true, f: r => r.is_egene ? '<span class="yes">yes</span>' : '<span class="no">no</span>' },
      { t: 'q', k: 'qval', n: true, f: r => fp(r.qval) },
      { t: 'Lead eSNP', k: 'lead_rsid', html: true, f: r => `<a href="#/variant/${r.lead_rsid}">${r.lead_rsid}</a>` },
      { t: 'Effect allele', k: 'effect_allele', f: r => r.effect_allele ? `${r.effect_allele}/${r.other_allele}` : 'uncertain' },
      { t: 'β', k: 'slope', n: true, f: r => f2(r.slope) },
      { t: 'P', k: 'pval_nominal', n: true, f: r => fp(r.pval_nominal) },
    ];
    const summary = computed(() => {
      if (!d.value) return '';
      const n = d.value.eqtl.filter(x => x.is_egene); const regs = new Set(n.map(x => x.region)); const cells = new Set(n.map(x => x.cell));
      return { n: n.length, regs: regs.size, cells: cells.size, three: d.value.evidence.filter(e => e.n_ctx_smr_heidi_h4 > 0).map(e => e.disease_name) };
    });
    return { d, err, ctx, lead, gridEq, gridEx, pickGrid, sexOpt, ageOpt, genoOpt, genoRs, esnps, genoStat, groups, geno, locOpt, loc, locErr, locDz, locCtx,
      evCols, eqCols, summary, fint, fp, f2, chip, ctxOptions: ctxOptions(), ctxOpts, esnpOpts, diseases: META.diseases, ctxLabel, locusHeight };
  },
  template: `<div>
    <div v-if="err" class="err">{{err}}</div>
    <div v-else-if="!d" class="load">Loading {{id}}…</div>
    <template v-else>
      <div class="ph"><span class="kind">Gene</span><h1>{{d.info.symbol}}</h1>
        <span class="meta">{{d.info.gene}} · {{d.info.biotype.replace(/_/g,' ')}} · chr{{d.info.chr}}:{{fint(d.info.start)}}–{{fint(d.info.end)}} ({{d.info.strand}}) · GRCh38</span></div>
      <p class="sum" v-if="summary.n">eGene in <b>{{summary.n}}</b> of 110 contexts ({{summary.regs}} regions, {{summary.cells}} cell types), {{fint(d.n_esnp_pairs)}} eSNP–gene pairs.
        <template v-if="summary.three.length"> SMR, HEIDI and colocalisation all support a link with <b>{{summary.three.join(', ')}}</b>.</template></p>
      <p class="sum" v-else>Tested in {{d.eqtl.length}} contexts; not an eGene in any (q ≥ 0.05).</p>

      <div class="cols" style="margin-bottom:18px">
        <section class="card"><h2>eQTL across the 110 contexts</h2><p class="cap">Lead eSNP of each context; colour is −log10 P signed by the direction of effect, shown where the gene is an eGene. Select a context to update the panels below.</p>
          <Chart :option="gridEq" :height="400" :name="d.info.symbol+'_eQTL_contexts'" @pick="pickGrid"/></section>
        <section class="card" v-if="gridEx"><h2>Expression across the 110 contexts</h2><p class="cap">Median deconvolved expression in each context.</p>
          <Chart :option="gridEx" :height="400" :name="d.info.symbol+'_expression_contexts'" @pick="pickGrid"/></section>
      </div>

      <section class="card">
        <div class="sel"><span>Context</span><Combo v-model="ctx" :options="ctxOpts" label="Context" width="230px"/>
          <span v-if="lead">{{lead.is_egene ? 'eGene' : 'not an eGene'}}, q {{fp(lead.qval)}}</span></div>
        <div class="cols3">
          <div><h2>By genotype<span v-if="genoRs"> of {{genoRs}}</span></h2>
            <div class="sel" v-if="esnps.length > 1"><span>eSNP</span><Combo v-model="genoRs" :options="esnpOpts" label="eSNP" width="300px"/>
              <span class="cap">{{esnps.length.toLocaleString()}} eSNPs, strongest first</span></div>
            <p class="cap" v-if="!genoRs">Not an eGene in this context, so it has no eSNPs.</p>
            <p class="cap" v-else-if="geno && geno.error">{{geno.error}}</p>
            <template v-else-if="genoOpt"><p class="cap">Genotype groups of {{genoRs}}<template v-if="genoStat">; β {{f2(genoStat.slope)}} per copy of {{genoStat.effect_allele || 'the effect allele'}}, P {{fp(genoStat.pval)}}</template>.</p>
              <Chart :option="genoOpt" :height="300" :name="d.info.symbol+'_'+genoRs+'_'+ctx"/></template>
            <p class="cap" v-else>Loading…</p></div>
          <div><h2>By sex</h2><p class="cap" v-if="groups && groups.error">{{groups.error}}</p><Chart v-else-if="sexOpt" :option="sexOpt" :height="300" :name="d.info.symbol+'_sex_'+ctx"/></div>
          <div><h2>By age group</h2><Chart v-if="ageOpt" :option="ageOpt" :height="300" :name="d.info.symbol+'_age_'+ctx"/></div>
        </div>
      </section>

      <section class="card" v-if="d.evidence.length"><h2>Disease evidence</h2><p class="cap">Number of contexts in which each test supports the gene for each disease (SMR with the corrected alleles; colocalisation by coloc).</p>
        <DataTable :rows="d.evidence" :cols="evCols"/></section>

      <section class="card" v-if="d.locus"><h2>Locus</h2>
        <p class="cap" v-if="!loc || loc.have_gwas">Every variant tested for {{d.info.symbol}} in the selected context (±1 Mb of the TSS) above, the disease GWAS below; colour is LD r² (1000 Genomes European) with the context's lead variant (purple diamond).</p>
        <p class="cap" v-else>Every variant tested for {{d.info.symbol}} in the selected context (±1 Mb of the TSS, tensorQTL nominal P). GWAS and LD are not in this build.</p>
        <div class="sel"><span>Context</span><Combo v-model="locCtx" :options="ctxOpts" label="Locus context" width="230px"/>
          <span>GWAS</span><select v-model="locDz" aria-label="Locus GWAS"><option v-for="z in diseases" :value="z.id">{{z.name}}</option></select></div>
        <p class="cap" v-if="loc && !loc.eqtl.length">{{d.info.symbol}} was not tested in {{ctxLabel(locCtx)}} (not expressed in this context).</p>
        <div v-if="locErr" class="err">{{locErr}}</div>
        <Chart v-else-if="locOpt" :option="locOpt" :height="locusHeight(loc)" :name="d.info.symbol+'_locus_'+locCtx+(loc.have_gwas ? '_'+locDz : '')"/></section>

      <section class="card"><h2>All contexts</h2><DataTable :rows="d.eqtl" :cols="eqCols" filter/></section>
    </template></div>`,
};

// rows of the gene track: label width estimated from the plot width (11 px font ≈ 6.2 px per character) so that names do not
// overlap; returns the genes with .row and the number of rows used
const MAX_GENE_ROWS = 6;
const LOCUS_PLOT_PX = () => Math.max(400, Math.min(window.innerWidth, 1240) - 120);
function geneRows(L, x0, x1) {
  const genes = L.genes.filter(g => g.end / 1e6 > x0 && g.start / 1e6 < x1).sort((a, b) => a.start - b.start);
  const bpPerPx = (x1 - x0) * 1e6 / LOCUS_PLOT_PX();
  const ends = [];
  genes.forEach(g => {
    const half = ((String(g.symbol).length + 2) * 6.2 / 2 + 4) * bpPerPx;
    const c = (g.start + g.end) / 2;
    const lo = Math.min(g.start, c - half), hi = Math.max(g.end, c + half);
    let k = 0; while (ends[k] != null && ends[k] > lo) k++;
    if (k < MAX_GENE_ROWS) { ends[k] = hi; g.row = k; g.hide = false; }
    else {   // no free row: bar in the row whose last gene ends first, name only in the tooltip
      let best = 0; for (let j = 1; j < MAX_GENE_ROWS; j++) if (ends[j] < ends[best]) best = j;
      g.row = best; g.hide = true; ends[best] = Math.max(ends[best], g.end);
    }
  });
  return { genes, nrows: Math.max(1, Math.min(ends.length, MAX_GENE_ROWS)) };
}
function locusHeight(L) {
  if (!L || !L.eqtl.length) return 400;
  const xs = L.eqtl.map(x => x.pos / 1e6);
  const { nrows } = geneRows(L, Math.min(...xs), Math.max(...xs));
  return (L.have_gwas === false ? 312 : 492) + nrows * 17 + 12 + 40;
}

function locusOption(L) {
  const col = (r2) => r2 == null ? '#C9CED8' : LDCOL.find(([, t]) => r2 < t)[0];
  const pt = (x) => ({ value: [x.pos / 1e6, nlp(x.p)], rsid: x.rsid, r2: x.r2, p: x.p,
    itemStyle: { color: x.rsid === L.index ? '#7A2C9E' : col(x.r2), borderColor: x.rsid === L.index ? C.ink : 'transparent' },
    symbol: x.rsid === L.index ? 'diamond' : 'circle', symbolSize: x.rsid === L.index ? 13 : 6, z: x.rsid === L.index ? 10 : 1 });
  const ord = (a) => a.slice().sort((x, y) => (x.r2 ?? -1) - (y.r2 ?? -1));
  const xs = L.eqtl.map(x => x.pos / 1e6); const x0 = xs.reduce((a, b) => Math.min(a, b), Infinity), x1 = xs.reduce((a, b) => Math.max(a, b), -Infinity);
  const gw = L.have_gwas !== false;  // false: eQTL track and genes only (full results, no GWAS / LD for this gene)
  // gene rows: a gene and its label must not overlap any other gene or label in the same row (PI 10-19「这些名字有时相互有遮挡」)
  const { genes, nrows } = geneRows(L, x0, x1);
  const TRACK = nrows * 17 + 12;
  const tip = (p) => Array.isArray(p.data) ? `<b>${p.data[3]}</b> (${p.data[4]} strand)` :
    (p.data.rsid ? `<b>${p.data.rsid}</b><br>P ${fp(p.data.p)}<br>r² ${p.data.r2 == null ? 'n/a' : f2(p.data.r2)}` : '');
  const geneItem = (params, api) => {
    const s = api.coord([api.value(0), api.value(2)]), e = api.coord([api.value(1), api.value(2)]);
    const me = api.value(3) === L.gene;
    const kids = [{ type: 'rect', shape: { x: s[0], y: s[1] - 3, width: Math.max(e[0] - s[0], 2), height: 6 }, style: { fill: me ? C.acc : '#7D879A' } }];
    if (me || !api.value(5)) kids.push({ type: 'text', style: { x: (s[0] + e[0]) / 2, y: s[1] + 6, text: api.value(3) + (api.value(4) === '+' ? ' →' : ' ←'),
      fill: me ? C.acc : C.mut, font: `${me ? 600 : 400} 11px ${FONT}`, align: 'center', verticalAlign: 'top' } });
    return { type: 'group', children: kids };
  };
  const ax = (i, name) => ({ gridIndex: i, type: 'value', name, nameLocation: 'middle', nameGap: 34, min: 0, splitLine: { lineStyle: { color: '#EEF1F5' } }, axisLabel: { color: C.mut } });
  const xa = (i, show) => ({ gridIndex: i, type: 'value', min: x0, max: x1, axisLabel: { show, color: C.mut, formatter: (v) => v.toFixed(2) },
    name: show ? `chr${L.chr} position (Mb, GRCh38)` : '', nameLocation: 'middle', nameGap: 28, splitLine: { show: false } });
  if (!gw) {
    const sig = (x) => ({ value: [x.pos / 1e6, nlp(x.p)], rsid: x.rsid, p: x.p, r2: null,
      itemStyle: { color: x.rsid === L.index ? '#7A2C9E' : (x.slope >= 0 ? C.up : C.dn), opacity: x.rsid === L.index ? 1 : 0.55, borderColor: x.rsid === L.index ? C.ink : 'transparent' },
      symbol: x.rsid === L.index ? 'diamond' : 'circle', symbolSize: x.rsid === L.index ? 13 : 5, z: x.rsid === L.index ? 10 : 1 });
    const tip2 = (p) => p.data.rsid ? `<b>${p.data.rsid}</b><br>P ${fp(p.data.p)}` : '';
    return {
      textStyle: baseText, animation: false,
      tooltip: { trigger: 'item', formatter: tip2, textStyle: { fontFamily: FONT } },
      title: [{ text: `eQTL · ${ctxLabel(L.ctx)} · ${fint(L.eqtl.length)} variants`, left: 64, top: 0, textStyle: { fontSize: 13, fontWeight: 600, fontFamily: FONT } }],
      legend: { right: 10, top: 0, data: [{ name: 'β > 0', itemStyle: { color: C.up } }, { name: 'β < 0', itemStyle: { color: C.dn } }], textStyle: { fontSize: 11, color: C.mut }, itemWidth: 10, itemHeight: 10 },
      grid: [{ left: 64, right: 16, top: 26, height: 250 }, { left: 64, right: 16, top: 312, height: TRACK }],
      xAxis: [xa(0, false), xa(1, true)],
      yAxis: [ax(0, '−log10 P'), { gridIndex: 1, type: 'value', min: -0.5, max: nrows, show: false, inverse: true }],
      series: [
        { type: 'scatter', xAxisIndex: 0, yAxisIndex: 0, large: L.eqtl.length > 5000, data: L.eqtl.map(sig) },
        { type: 'custom', xAxisIndex: 1, yAxisIndex: 1, data: genes.map(g => [g.start / 1e6, g.end / 1e6, g.row, g.symbol, g.strand, g.hide ? 1 : 0]), renderItem: geneItem },
        { type: 'scatter', name: 'β > 0', data: [], itemStyle: { color: C.up } }, { type: 'scatter', name: 'β < 0', data: [], itemStyle: { color: C.dn } },
      ],
    };
  }
  return {
    textStyle: baseText, animation: false,
    tooltip: { trigger: 'item', formatter: tip, textStyle: { fontFamily: FONT } },
    title: [{ text: `eQTL · ${ctxLabel(L.ctx)}`, left: 64, top: 0, textStyle: { fontSize: 13, fontWeight: 600, fontFamily: FONT } },
            { text: `GWAS · ${META.diseases.find(z => z.id === L.disease).name}`, left: 64, top: 238, textStyle: { fontSize: 13, fontWeight: 600, fontFamily: FONT } }],
    legend: { right: 10, top: 0, data: LDCOL.map(([c], i) => ({ name: ['r² < 0.2', '0.2–0.4', '0.4–0.6', '0.6–0.8', '≥ 0.8'][i], itemStyle: { color: c } })), textStyle: { fontSize: 11, color: C.mut }, itemWidth: 10, itemHeight: 10 },
    grid: [{ left: 64, right: 16, top: 26, height: 190 }, { left: 64, right: 16, top: 262, height: 190 }, { left: 64, right: 16, top: 492, height: TRACK }],
    xAxis: [xa(0, false), xa(1, false), xa(2, true)],
    yAxis: [ax(0, '−log10 P'), ax(1, '−log10 P'), { gridIndex: 2, type: 'value', min: -0.5, max: nrows, show: false, inverse: true }],
    series: [
      { type: 'scatter', xAxisIndex: 0, yAxisIndex: 0, data: ord(L.eqtl).map(pt) },
      { type: 'scatter', xAxisIndex: 1, yAxisIndex: 1, data: ord(L.gwas).map(pt),
        markLine: { silent: true, symbol: 'none', lineStyle: { color: C.up, type: 'dashed', width: 1 }, label: { show: false }, data: [{ yAxis: -Math.log10(5e-8) }] } },
      { type: 'custom', xAxisIndex: 2, yAxisIndex: 2, data: genes.map(g => [g.start / 1e6, g.end / 1e6, g.row, g.symbol, g.strand, g.hide ? 1 : 0]), renderItem: geneItem },
      ...LDCOL.map(([c], i) => ({ type: 'scatter', name: ['r² < 0.2', '0.2–0.4', '0.4–0.6', '0.6–0.8', '≥ 0.8'][i], data: [], itemStyle: { color: c } })),
    ],
  };
}

const VariantPage = {
  components: { Chart, DataTable, Combo },
  props: ['id'],
  setup(props) {
    const d = ref(null); const err = ref(''); const gene = ref(null); const ctx = ref(null); const geno = ref(null);
    const tested = ref(null); const testedErr = ref('');
    const load = async () => {
      d.value = null; err.value = ''; tested.value = null; testedErr.value = '';
      try {
        // the full-results request starts together with the variant request instead of after it (#146 speed-up)
        const tp = META.full ? api(`/variant/${props.id}/tested`).then(r => r.rows).catch(e => { testedErr.value = e.message; return null; }) : null;
        d.value = await api('/variant/' + props.id);
        if (tp) tested.value = await tp;
        const first = d.value.pairs[0] || (tested.value || [])[0];
        if (first) { gene.value = first.symbol || first.gene; ctx.value = first.ctx; }
      } catch (e) { err.value = e.message; }
    };
    watch(() => props.id, load, { immediate: true });
    // significant pairs (esnp_pairs) and, with the full results, every tested pair
    const egenes = computed(() => d.value ? [...new Set(d.value.pairs.map(p => p.symbol))] : []);
    const genes = computed(() => {
      if (!d.value) return [];
      const t = tested.value ? [...new Set(tested.value.map(p => p.symbol || p.gene))] : [];
      return [...new Set([...egenes.value, ...t])];
    });
    const rows = computed(() => d.value ? d.value.pairs.filter(p => p.symbol === gene.value) : []);
    const trows = computed(() => tested.value ? tested.value.filter(p => (p.symbol || p.gene) === gene.value) : []);
    const ea = computed(() => (d.value && d.value.info && d.value.info.effect_allele) || (rows.value[0] && rows.value[0].effect_allele) || 'the effect allele');
    const grid = computed(() => {
      if (!d.value) return null;
      const full = tested.value != null;
      const src = full ? trows.value : rows.value;
      const vals = Object.fromEntries(src.map(p => [p.ctx, p.slope]));
      const m = Math.max(0.5, ...Object.values(vals).filter(v => v != null).map(Math.abs));
      const lead = new Set(d.value.lead_of.map(x => x.ctx));
      const marks = full ? new Set(trows.value.filter(p => p.esnp).map(p => p.ctx)) : lead;
      return gridOption({ values: vals, scale: 'div', min: -m, max: m, marks, selected: ctx.value, legend: `β per copy of ${ea.value}`,
        tip: (k) => {
          const p = src.find(x => x.ctx === k);
          if (!p) return full ? `${ctxLabel(k)}<br>${gene.value} not tested here (not expressed in this context)` : `${ctxLabel(k)}<br>not significant for ${gene.value} here (this version holds significant pairs only)`;
          const s = full ? (p.esnp ? 'eSNP (significant)' : 'not significant') : (p.independent ? 'independent eSNP' : '');
          return `<b>${ctxLabel(k)}</b><br>β ${f2(p.slope)}, P ${fp(p.p ?? p.pval)}${s ? '<br>' + s : ''}${lead.has(k) ? '<br>lead eSNP' : ''}${full && p.independent ? '<br>independent' : ''}`;
        } });
    });
    const loadGeno = async () => {
      geno.value = null;
      if (!d.value || !d.value.geno_genes.length) return;
      const pr = d.value.pairs.find(p => p.symbol === gene.value);
      if (!pr) return;
      const gid = pr.gene;
      if (!d.value.geno_genes.includes(gid)) return;
      try { geno.value = await api(`/genotype?gene=${gid}&rsid=${props.id}`); } catch (e) { geno.value = { error: e.message }; return; }
      // open on the context where the pair is strongest among those with genotype groups
      const have = new Set(geno.value.boxes.map(b => b.ctx));
      if (!have.has(ctx.value)) { const best = rows.value.filter(r => have.has(r.ctx)).sort((a, b) => a.pval - b.pval)[0]; if (best) ctx.value = best.ctx; }
    };
    watch([gene, () => d.value], loadGeno);
    watch(gene, () => { const r = rows.value[0] || trows.value.slice().sort((a, b) => a.p - b.p)[0]; if (r) ctx.value = r.ctx; });
    const genoOpt = computed(() => {
      if (!geno.value || !geno.value.boxes) return null;
      const b = geno.value.boxes.filter(x => x.ctx === ctx.value).sort((a, c) => a.n_effect - c.n_effect);
      if (!b.length) return null;
      const ea = b[0].effect_allele, oa = b[0].other_allele;
      const st = (geno.value.stats || []).find(s => s.ctx === ctx.value);
      const note = st ? `${ctxLabel(ctx.value)}: β = ${f2(st.slope)}${ea ? ` per copy of ${ea}` : ''}, P = ${fp(st.pval)}` : '';
      return boxOption(b.map(x => ({ ...x, lo: x.whisker_lo, hi: x.whisker_hi })), { colorOf: (g) => ['#9DBCDB', '#7A5EC0', '#2E1D5C'][g.n_effect],
        labelOf: (g) => ea && oa ? [oa + oa, oa + ea, ea + ea][g.n_effect] + ` (${g.n})` : `${g.n_effect} (${g.n})`, xname: props.id, note });
    });
    const genoCtx = computed(() => geno.value && geno.value.boxes ? [...new Set(geno.value.boxes.map(b => b.ctx))] : []);
    const pick = (e) => { if (e.data) ctx.value = e.data.ctx; };
    const cols = [
      { t: 'Context', k: 'ctx', html: true, f: r => chip(r.ctx) },
      { t: 'Gene', k: 'symbol', html: true, f: r => `<a href="#/gene/${r.symbol}">${r.symbol}</a>` },
      { t: 'Distance to TSS', k: 'dist', n: true, f: r => fint(r.dist) },
      { t: 'β', k: 'slope', n: true, f: r => f2(r.slope) },
      { t: 'P', k: 'pval', n: true, f: r => fp(r.pval) },
      { t: 'Independent', k: 'independent', html: true, f: r => r.independent ? '<span class="yes">yes</span>' : (notIn1000G.value ? '<span class="no">not clumped</span>' : '<span class="no">no</span>') },
    ];
    // LD clumping used 1000G EUR: a variant absent from 1000G could not be clumped and is never marked independent (PI 10-19 #133)
    const notIn1000G = computed(() => !!(d.value && d.value.info && (d.value.info.in_1000g === false || d.value.info.allele_status === 'not_in_1000G')));
    const tcols = [
      { t: 'Context', k: 'ctx', html: true, f: r => chip(r.ctx) },
      { t: 'Gene', k: 'symbol', html: true, f: r => `<a href="#/gene/${r.symbol || r.gene}">${r.symbol || r.gene}</a>` },
      { t: 'β', k: 'slope', n: true, f: r => f2(r.slope) },
      { t: 'SE', k: 'se', n: true, f: r => f2(r.se, 3) },
      { t: 'P', k: 'p', n: true, f: r => fp(r.p) },
      { t: 'eSNP', k: 'esnp', html: true, f: r => r.esnp ? '<span class="yes">yes</span>' : '<span class="no">no</span>' },
    ];
    // header line: chr / position / alleles from the variant index when present, otherwise from the first significant pair
    const p0 = computed(() => {
      if (!d.value) return null;
      const i = d.value.info, p = d.value.pairs[0] || {};
      const maf = p.MAF != null ? p.MAF : (tested.value && tested.value[0] && tested.value[0].af != null ? Math.min(tested.value[0].af, 1 - tested.value[0].af) : null);
      return { chr: i ? String(i.chr).replace('chr', '') : p.chr, pos: i ? i.pos : p.pos, effect_allele: i ? i.effect_allele : p.effect_allele, other_allele: i ? i.other_allele : p.other_allele, MAF: maf };
    });
    const geneOpts = computed(() => genes.value.map(g => ({ value: g, label: g, sub: egenes.value.includes(g) ? 'eSNP' : 'not significant' })));
    const nTestedGenes = computed(() => tested.value ? new Set(tested.value.map(r => r.gene)).size : 0);
    const nTestedCtx = computed(() => tested.value ? new Set(tested.value.map(r => r.ctx)).size : 0);
    return { d, err, gene, genes, egenes, rows, trows, tested, testedErr, grid, pick, ctx, genoOpt, genoCtx, geno, cols, tcols, p0, nTestedGenes, nTestedCtx,
      fint, fp, f2, ctxLabel, notIn1000G, geneOpts, hasGeno: META.geno_box };
  },
  template: `<div>
    <div v-if="err" class="err">{{err}}</div>
    <div v-else-if="!d" class="load">Loading {{id}}…</div>
    <template v-else>
      <div class="ph"><span class="kind">Variant</span><h1>{{d.rsid}}</h1>
        <span class="meta">chr{{p0.chr}}:{{fint(p0.pos)}} · {{p0.effect_allele ? p0.effect_allele+' (effect) / '+p0.other_allele : 'alleles uncertain'}} · MAF {{f2(p0.MAF,3)}}</span></div>
      <p class="sum"><template v-if="d.pairs.length">eSNP of <b>{{egenes.length}}</b> gene{{egenes.length>1?'s':''}} ({{egenes.join(', ')}}) in {{new Set(d.pairs.map(p=>p.ctx)).size}} contexts; lead eSNP in {{d.lead_of.length}}.</template>
        <template v-else>Not an eSNP: no gene passes its significance threshold for this variant in any context.</template>
        <template v-if="tested"> Tested with {{nTestedGenes}} gene{{nTestedGenes>1?'s':''}} in {{nTestedCtx}} contexts ({{fint(tested.length)}} tests).</template>
        <template v-if="d.smr.length"> Top SNP of {{d.smr.length}} significant SMR tests.</template></p>
      <div class="cols" style="margin-bottom:18px">
        <section class="card"><div class="sel"><span>Gene</span><Combo v-model="gene" :options="geneOpts" label="Gene" width="260px"/></div>
          <h2>Effect on {{gene}} across contexts</h2><p class="cap" v-if="tested">Every context where the pair was tested; outlined: significant (eSNP). Grey: gene not tested in that context.</p>
          <p class="cap" v-else>Outlined: contexts where {{d.rsid}} is the lead eSNP.</p>
          <Chart :option="grid" :height="400" :name="d.rsid+'_'+gene+'_contexts'" @pick="pick"/></section>
        <section class="card"><h2>{{gene}} expression by genotype</h2>
          <div class="sel" v-if="genoCtx.length"><span>Context</span><Combo v-model="ctx" :options="genoCtx.map(c => ({ value: c, label: ctxLabel(c) }))" label="Context" width="230px"/></div>
          <p class="cap" v-if="!hasGeno">Genotype groups are not in this build.</p>
          <p class="cap" v-else-if="!genoCtx.length">No genotype groups for {{d.rsid}} and {{gene}} in this build.</p>
          <p class="cap" v-else-if="!genoOpt">{{ctxLabel(ctx)}} has no genotype groups for this pair; choose a context from the list.</p>
          <Chart v-else :option="genoOpt" :height="340" :name="d.rsid+'_'+gene+'_'+ctx"/></section>
      </div>
      <section class="card" v-if="d.pairs.length"><h2>All eSNP–gene pairs</h2>
        <p class="cap">Significant pairs only<template v-if="tested">; every tested pair is in the next table</template>. Independent: left after LD clumping, done in each context on its own (PLINK, 250 kb, r² 0.8, 1000 Genomes European).<template v-if="notIn1000G"> <b>{{d.rsid}} is not in 1000 Genomes, so it could not be clumped ("not clumped").</b></template> eSNPs in LD form one group and the one with the smallest P in that context is kept, so the same variant can be independent in one context and merged into a stronger variant in another.</p>
        <DataTable :rows="d.pairs" :cols="cols" filter/></section>
      <section class="card" v-if="tested"><h2>All tested variant–gene pairs</h2>
        <p class="cap">Every gene whose ±1 Mb cis window contains {{d.rsid}}, in every context where the gene was tested (tensorQTL nominal results). eSNP: P below the gene's permutation-derived threshold in that context.</p>
        <DataTable :rows="tested" :cols="tcols" filter/></section>
      <div v-if="testedErr" class="err">{{testedErr}}</div>
    </template></div>`,
};

const DiseasePage = {
  components: { Chart, DataTable },
  props: ['id'],
  setup(props) {
    const d = ref(null); const err = ref(''); const gsel = ref(null);
    const load = async () => { d.value = null; err.value = ''; try { d.value = await api('/disease/' + props.id); const g = d.value.genes[0]; gsel.value = g ? g.gene : null; } catch (e) { err.value = e.message; } };
    watch(() => props.id, load, { immediate: true });
    const fgrid = computed(() => {
      if (!d.value) return null;
      const vals = {}; const marks = new Set();
      d.value.fisher.forEach(x => { if (x.OR > 0 && x.E_N > 0) vals[x.ctx] = Math.log2(x.OR); if (x.significant) marks.add(x.ctx); });
      const m = Math.max(1, ...Object.values(vals).map(Math.abs));
      return gridOption({ values: vals, scale: 'div', min: -m, max: m, marks, legend: 'log2 odds ratio',
        tip: (k) => { const x = d.value.fisher.find(f => f.ctx === k); return x ? `<b>${ctxLabel(k)}</b><br>OR ${f2(x.OR)} (${f2(x.OR_lower95)}–${f2(x.OR_upper95)})<br>q ${fp(x.q_value)}${x.significant ? ' · significant' : ''}<br>eSNPs at GWAS loci ${x.E_N}` : ctxLabel(k); } });
    });
    const sgrid = computed(() => {
      if (!d.value || !gsel.value) return null;
      const rows = d.value.smr.filter(x => x.gene === gsel.value);
      const vals = Object.fromEntries(rows.map(x => [x.ctx, nlp(x.p_SMR)]));
      const marks = new Set(rows.filter(x => x.smr_heidi_coloc).map(x => x.ctx));
      return gridOption({ values: vals, min: 0, marks, legend: '−log10 P (SMR), q<0.05 shown',
        tip: (k) => { const x = rows.find(r => r.ctx === k); return x ? `<b>${ctxLabel(k)}</b><br>P ${fp(x.p_SMR)}, q ${fp(x.q_SMR)}<br>β ${f2(x.b_SMR)} · HEIDI P ${fp(x.p_HEIDI)} · PP.H4 ${f2(x.PP_H4_abf)}` : `${ctxLabel(k)}<br>SMR q ≥ 0.05`; } });
    });
    const gsym = computed(() => d.value && gsel.value ? d.value.genes.find(g => g.gene === gsel.value)?.symbol : '');
    const gcols = [
      { t: 'Gene', k: 'symbol', html: true, f: r => `<a href="#/gene/${r.symbol}">${r.symbol}</a>${r.MHC ? '<span class="tag">MHC</span>' : ''}` },
      { t: 'SMR q<0.05', k: 'n_ctx_smr', n: true },
      { t: 'HEIDI P≥0.01', k: 'n_ctx_heidi_ge_0_01', n: true },
      { t: 'PP.H4≥0.8', k: 'n_ctx_h4_ge_0_8', n: true },
      { t: 'All three', k: 'n_ctx_smr_heidi_h4', n: true, html: true, f: r => r.n_ctx_smr_heidi_h4 > 0 ? `<span class="yes">${r.n_ctx_smr_heidi_h4}</span>` : '<span class="no">0</span>' },
      { t: 'Regions', k: 'n_regions', n: true },
      { t: 'Cell types', k: 'n_cells', n: true },
      { t: 'Best q', k: 'min_q_smr', n: true, f: r => fp(r.min_q_smr) },
      { t: 'scGWAS', k: 'scgwas_module', html: true, f: r => r.scgwas_module ? '<span class="yes">yes</span>' : '<span class="no">no</span>' },
    ];
    const mcols = [
      { t: 'Context', k: 'ctx', html: true, f: r => chip(r.ctx) },
      { t: 'Seed gene', k: 'seed_symbol', html: true, f: r => `<a href="#/gene/${r.seed_symbol}">${r.seed_symbol}</a>` },
      { t: 'Module genes', k: 'module_symbols', f: r => r.module_symbols.join(', ') },
      { t: 'Score', k: 'module_score', n: true, f: r => f2(r.module_score) },
      { t: 'P (permutation)', k: 'p', n: true, f: r => r.p === 0 ? '0 *' : fp(r.p) },
    ];
    const nsig = computed(() => d.value ? d.value.fisher.filter(x => x.significant).length : 0);
    const n3 = computed(() => d.value ? d.value.genes.filter(x => x.n_ctx_smr_heidi_h4 > 0) : []);
    const pickGene = (r) => { gsel.value = r.gene; };
    return { d, err, fgrid, sgrid, gsel, gsym, gcols, mcols, nsig, n3, pickGene, fint };
  },
  template: `<div>
    <div v-if="err" class="err">{{err}}</div>
    <div v-else-if="!d" class="load">Loading…</div>
    <template v-else>
      <div class="ph"><span class="kind">Disease</span><h1>{{d.name}}</h1><span class="meta">GWAS {{d.gwas}}</span></div>
      <p class="sum">eSNPs enriched at GWAS loci in <b>{{nsig}}</b> of 110 contexts. <b>{{d.genes.length}}</b> genes pass SMR (q &lt; 0.05);
        <b>{{n3.length}}</b> also pass HEIDI and colocalisation<template v-if="n3.length"> ({{n3.map(g=>g.symbol).join(', ')}})</template>.
        {{d.scgwas.length}} significant scGWAS modules.</p>
      <div class="cols" style="margin-bottom:18px">
        <section class="card"><h2>Enrichment at GWAS loci</h2><p class="cap">Fisher's exact test, eSNPs vs non-eSNPs at loci of "{{d.trait}}" in the GWAS Catalog; outlined cells are significant.</p>
          <Chart :option="fgrid" :height="400" :name="d.id+'_Fisher_contexts'"/></section>
        <section class="card"><h2>SMR for {{gsym || '—'}}</h2><p class="cap">Select a gene in the table below. Outlined: SMR, HEIDI and colocalisation all pass.</p>
          <Chart v-if="sgrid" :option="sgrid" :height="400" :name="d.id+'_'+gsym+'_SMR_contexts'"/><p v-else class="empty">No gene passes SMR for this disease.</p></section>
      </div>
      <section class="card"><h2>Genes</h2><p class="cap">Counts are numbers of contexts. Select a row to draw its SMR results above.</p>
        <DataTable :rows="d.genes" :cols="gcols" filter :pick="pickGene" :picked="(r)=>r.gene===gsel"/></section>
      <section class="card" v-if="d.scgwas.length"><h2>scGWAS modules</h2><p class="cap">* No permuted module scored as high as the observed one.</p><DataTable :rows="d.scgwas" :cols="mcols" filter/></section>
    </template></div>`,
};

const ContextPage = {
  components: { Chart, DataTable },
  props: ['id'],
  setup(props) {
    const d = ref(null); const err = ref('');
    const load = async () => { d.value = null; err.value = ''; try { d.value = await api('/context/' + props.id); } catch (e) { err.value = e.message; } };
    watch(() => props.id, load, { immediate: true });
    const frOpt = computed(() => d.value ? boxOption(d.value.fractions.map(f => ({ ...f, group: f.cell })), { yname: 'Estimated fraction', colorOf: (g) => PAL_CELL[g.cell] }) : null);
    const shOpt = computed(() => {
      if (!d.value) return null;
      const s = d.value.sharing.filter(x => x.jaccard != null).slice(0, 15).reverse();
      return { textStyle: baseText, animation: false, tooltip: { textStyle: { fontFamily: FONT } }, grid: { left: 110, right: 30, top: 8, bottom: 36 },
        xAxis: { type: 'value', name: 'Jaccard index of eGenes', nameLocation: 'middle', nameGap: 24, splitLine: { lineStyle: { color: '#EEF1F5' } }, axisLabel: { color: C.mut } },
        yAxis: { type: 'category', data: s.map(x => ctxLabel(x.ctx)), axisTick: { show: false }, axisLabel: { color: C.ink } },
        series: [{ type: 'bar', barWidth: '62%', data: s.map(x => ({ value: x.jaccard, ctx: x.ctx, itemStyle: { color: PAL_REG[x.ctx.split('.')[0]] } })),
          label: { show: true, position: 'right', formatter: (p) => f2(p.value), fontSize: 11, color: C.mut } }] };
    });
    const cols = [
      { t: 'Gene', k: 'symbol', html: true, f: r => `<a href="#/gene/${r.symbol}">${r.symbol}</a>${r.biotype !== 'protein_coding' ? '<span class="tag v">' + r.biotype.replace(/_/g, ' ') + '</span>' : ''}` },
      { t: 'Lead eSNP', k: 'lead_rsid', html: true, f: r => `<a href="#/variant/${r.lead_rsid}">${r.lead_rsid}</a>` },
      { t: 'β', k: 'slope', n: true, f: r => f2(r.slope) },
      { t: 'P', k: 'pval_nominal', n: true, f: r => fp(r.pval_nominal) },
      { t: 'q', k: 'qval', n: true, f: r => fp(r.qval) },
    ];
    const pickBar = (e) => e.data && e.data.ctx && go('/context/' + e.data.ctx);
    return { d, err, frOpt, shOpt, cols, pickBar, fint, f2, PAL_REG };
  },
  template: `<div>
    <div v-if="err" class="err">{{err}}</div>
    <div v-else-if="!d" class="load">Loading…</div>
    <template v-else>
      <div class="ph"><span class="kind">Context</span><h1><span class="sw" :style="{background:PAL_REG[d.context.region],width:'18px',height:'18px',marginRight:'10px'}"></span>{{d.region_name}} · {{d.cell_name}}</h1></div>
      <p class="sum"><b>{{fint(d.context.egenes)}}</b> eGenes ({{fint(d.context.egenes_lncrna)}} non-coding) and <b>{{fint(d.context.independent_esnps)}}</b> independent eSNPs among {{fint(d.context.genes_tested)}} genes tested in {{d.context.n_samples}} samples. Mean estimated {{d.cell_name.toLowerCase()}} fraction {{f2(d.context.mean_frac,3)}}.</p>
      <div class="cols" style="margin-bottom:18px">
        <section class="card"><h2>Cell composition of {{d.region_name.toLowerCase()}} samples</h2><p class="cap">Fractions estimated by deconvolution (BayesPrism / InstaPrism with the STAB reference).</p>
          <Chart :option="frOpt" :height="320" :name="d.context.ctx+'_fractions'"/></section>
        <section class="card"><h2>Contexts sharing the most eGenes</h2><p class="cap">Select a bar to open that context.</p>
          <Chart :option="shOpt" :height="320" :name="d.context.ctx+'_sharing'" @pick="pickBar"/></section>
      </div>
      <section class="card"><h2>eGenes</h2><DataTable :rows="d.egenes" :cols="cols" filter/></section>
    </template></div>`,
};

const Browse = {
  components: { Chart },
  setup() {
    const metric = ref('egenes');
    const opt = computed(() => {
      const k = metric.value; const lab = { egenes: 'eGenes', independent_esnps: 'independent eSNPs', mean_frac: 'mean cell fraction', n_samples: 'samples' }[k];
      const values = Object.fromEntries(META.contexts.map(c => [c.ctx, c[k]]));
      return gridOption({ values, legend: lab, tip: (c, v) => `<b>${ctxLabel(c)}</b><br>${k === 'mean_frac' ? f2(v, 3) : fint(v)} ${lab}` });
    });
    return { metric, opt, pick: (e) => e.data && go('/context/' + e.data.ctx) };
  },
  template: `<section class="card"><div class="cg-h"><h2>The 110 contexts</h2>
    <div class="sel" style="margin:0"><button v-for="(l,k) in {egenes:'eGenes',independent_esnps:'Independent eSNPs',mean_frac:'Cell fraction',n_samples:'Samples'}" class="btn" :class="{on:metric===k}" @click="metric=k">{{l}}</button></div></div>
    <Chart :option="opt" :height="460" name="BrainCellReg_contexts" @pick="pick"/></section>`,
};

const Diseases = {
  setup() { return { ds: META.diseases }; },
  template: `<section class="card prose"><h2>The ten diseases</h2><p class="cap">GWAS used for SMR, colocalisation and scGWAS.</p>
    <div class="tw" style="max-height:none"><table><tbody><tr v-for="d in ds"><td><a :href="'#/disease/'+d.id">{{d.name}}</a></td><td class="n">{{d.gwas}}</td></tr></tbody></table></div></section>`,
};

const Download = {
  setup() {
    const rows = [
      ['contexts', 'The 110 contexts: samples, cell fraction, genes tested, eGenes, independent eSNPs'],
      ['egene_stats', 'Permutation result of every gene in every context (lead variant, β, q)'],
      ['esnp_pairs', 'All eSNP–gene pairs of eGenes'],
      ['genes', 'Genes tested, with GRCh38 coordinates'],
      ['samples', 'Samples: region, sex, age group, estimated cell fractions'],
      ['fisher', 'Fisher enrichment at GWAS Catalog loci'],
      ['smr', 'SMR, HEIDI and colocalisation for the ten diseases'],
      ['gene_evidence', 'Disease × gene summary of the integration'],
      ['scgwas_modules', 'Significant scGWAS modules'],
    ];
    return { rows };
  },
  template: `<section class="card prose" style="max-width:none"><h2>Download</h2><p class="cap">Parquet files; open with DuckDB, pandas (pyarrow) or R (arrow).</p>
    <div class="tw" style="max-height:none"><table><tbody><tr v-for="r in rows"><td><a :href="'/api/download/'+r[0]">{{r[0]}}.parquet</a></td><td style="white-space:normal">{{r[1]}}</td></tr></tbody></table></div></section>`,
};

const About = {
  template: `<section class="card prose"><h2>About</h2>
    <p><b>BrainCellReg</b> (Brain Cell-type-specific Gene Regulation Resource) holds the results of a study of gene regulation in the human brain. Expression from 1,231 post-mortem samples of 134 donors (UK Brain Expression Consortium, ten regions) was split into eleven cell types by deconvolution with a single-cell reference, giving 110 region × cell-type contexts. Each context was tested for cis-eQTL with tensorQTL, and the results were integrated with genome-wide association studies of ten diseases by Fisher enrichment, SMR with HEIDI, colocalisation and scGWAS.</p>
    <p>Cell-type contexts are defined by estimated cell composition, not by sorted cells or single-nucleus sequencing.</p>
    <p>Citation and contact will be added on publication.</p></section>`,
};

/* ---------- Help (PI 10-16 #93 「好像还差一个使用说明的页面」) ---------- */
const Help = {
  setup() { return { full: META.full }; },
  template: `<section class="card prose help"><h2>How to use BrainCellReg</h2>
    <h3>Search</h3>
    <p>Type in the search box on the home page or in the top bar. It accepts a <b>gene</b> (symbol such as GPNMB, or Ensembl id), an <b>rs number</b> (such as rs156429) or a <b>disease</b> (such as Parkinson). Suggestions appear as you type; press Enter or pick one.</p>
    <p v-if="full">Any of the 6.4 million genotyped variants can be searched by rs number, eSNP or not. The site holds the full results of every tested variant–gene pair (about 7.8 billion tests: every gene × every variant within 1 Mb of its TSS × every context where the gene was tested), so non-significant pairs can be looked up too.</p>
    <p v-else>An rs number is found when the variant is an eSNP of at least one eGene (see Terms). Variants that were tested but are not eSNPs are not in this build; the full results of every tested variant–gene pair (about 7.8 billion tests, 176 GB) are on the online version.</p>
    <h3>Pages</h3>
    <ul>
      <li><b>Home</b>: eGenes in each of the 110 contexts. Click a cell to open that context.</li>
      <li><b>Gene</b>: the gene's lead eSNP in every context (colour: signed −log10 P, only where the gene is an eGene) and its expression in every context. Click a cell to choose a context; the panels below then show expression by genotype of the lead eSNP, by sex and by age group, the disease evidence (SMR, HEIDI, colocalisation, scGWAS) and, where available, the locus plot of the eQTL above the disease GWAS. The table at the end lists all 110 contexts.</li>
      <li v-if="full"><b>Variant</b>: the effect of the variant on each gene it was tested with, in each context (β per copy of the effect allele; outlined: significant, i.e. an eSNP), expression by genotype for significant pairs, its eSNP–gene pairs and every tested pair. Genes it is not significant for are marked "(not significant)" in the gene list.</li>
      <li v-else><b>Variant</b>: the effect of the variant on each gene in each context (β per copy of the effect allele; outlined: contexts where it is the lead eSNP), expression by genotype, and all of its eSNP–gene pairs.</li>
      <li><b>Disease</b>: Fisher enrichment of eSNPs at the disease's GWAS loci in each context (outlined: significant), SMR results of a selected gene (outlined: SMR, HEIDI and colocalisation all pass), the evidence table of every gene, and the scGWAS modules.</li>
      <li><b>Context</b>: estimated cell composition of the region's samples, the contexts sharing the most eGenes (Jaccard index), and the eGenes of the context.</li>
      <li><b>Download</b>: every table as a Parquet file (DuckDB, pandas with pyarrow, or R with arrow).</li>
    </ul>
    <h3>Figures and tables</h3>
    <p>Every figure is drawn in the browser from the data. The <b>SVG</b> and <b>PNG</b> buttons above a figure save it. Tables can be filtered with the box above them; long tables show the first 500 rows, use Download for the rest.</p>
    <p>Colours: red = higher expression with the effect allele (or enrichment), blue = lower (or depletion). Regions are listed cortex first and cerebellum last, in the same order and colours as in the paper.</p>
    <h3>Terms</h3>
    <ul>
      <li><b>Context</b>: one brain region × one cell type (10 × 11 = 110). Cell-type expression was estimated from bulk tissue by deconvolution (BayesPrism / InstaPrism with a single-cell reference), not from sorted cells.</li>
      <li><b>eGene</b>: a gene with a significant cis-eQTL in a context (tensorQTL permutation, q &lt; 0.05).</li>
      <li><b>eSNP</b>: a variant whose nominal P for an eGene is below that gene's permutation-derived threshold in the context. <b>Lead eSNP</b>: the most significant one. <b>Independent eSNPs</b>: eSNPs left after LD clumping in each context (PLINK, 250 kb, r² 0.8, 1000 Genomes European); in each group of eSNPs in LD the one with the smallest P in that context is kept, so whether a variant is independent can differ between contexts. Every significant pair is in the database, independent or not. Clumping needs the variant in the 1000 Genomes reference; eSNPs absent from it could not be clumped and are labelled "not clumped".</li>
      <li><b>β</b>: effect per copy of the effect allele on deconvolved expression. The effect allele is the allele counted in the genotype dosages; its letter was matched to European reference frequencies (1000 Genomes, then dbSNP ALFA; either strand), or to LD with neighbouring variants. "uncertain" marks variants whose letter could not be matched: absent from both references, letters that differ from the reference alleles (including A/T and C/G variants, whose strand cannot be told), or allele frequencies too close to 0.5 to decide.</li>
      <li><b>SMR / HEIDI</b>: summary-data-based Mendelian randomisation between expression and disease (q &lt; 0.05), and the HEIDI test of a shared causal variant (P ≥ 0.01 passes). <b>Colocalisation</b>: coloc posterior PP.H4 ≥ 0.8.</li>
      <li><b>Genotype boxes</b>: summary statistics of each genotype group; no individual-level genotype is stored on this site.</li>
    </ul>
    <h3>Examples</h3>
    <p><a href="#/gene/GPNMB">GPNMB</a> (Parkinson's disease) · <a href="#/variant/rs156429">rs156429</a> · <a href="#/disease/PD">Parkinson's disease</a> · <a href="#/context/FCTX.ExN">Frontal cortex · Excitatory neurons</a></p>
    </section>`,
};

/* ---------- app shell with hash routing ---------- */
const App = {
  components: { Search, Home, GenePage, VariantPage, DiseasePage, ContextPage, CellPage, Browse, Diseases, Download, Help, About },
  setup() {
    const route = ref(location.hash.slice(1) || '/'); const ready = ref(false); const err = ref('');
    window.addEventListener('hashchange', () => { route.value = location.hash.slice(1) || '/'; window.scrollTo(0, 0); });
    onMounted(async () => { try { META = await api('/meta'); ready.value = true; } catch (e) { err.value = 'The database is not reachable: ' + e.message; } });
    const view = computed(() => {
      const [, kind, id] = route.value.split('/');
      const t = { gene: 'GenePage', variant: 'VariantPage', disease: 'DiseasePage', context: 'ContextPage', celltype: 'CellPage', browse: 'Browse', diseases: 'Diseases', download: 'Download', help: 'Help', about: 'About' }[kind] || 'Home';
      return { t, id: id ? decodeURIComponent(id) : null, kind: kind || '' };
    });
    watch(view, (v) => { document.title = v.id ? `${v.id} · BrainCellReg` : 'BrainCellReg'; });
    return { view, ready, err };
  },
  template: `<header class="bar"><div class="bar-in">
      <a class="mark" href="#/">BrainCell<b>Reg</b></a>
      <nav class="nav"><a href="#/browse" :class="{on:view.kind==='browse'}">Contexts</a><a href="#/diseases" :class="{on:view.kind==='disease'||view.kind==='diseases'}">Diseases</a>
        <a href="#/download" :class="{on:view.kind==='download'}">Download</a><a href="#/help" :class="{on:view.kind==='help'}">Help</a><a href="#/about" :class="{on:view.kind==='about'}">About</a></nav>
      <Search v-if="ready && view.t!=='Home'"/></div></header>
    <main><div v-if="err" class="err">{{err}}</div><div v-else-if="!ready" class="load">Loading…</div>
      <component v-else :is="view.t" :id="view.id" :key="view.t+(view.id||'')"/></main>
    <footer class="foot">BrainCellReg · Brain Cell-type-specific Gene Regulation Resource · results of the cell-type-resolved brain eQTL study · draft build for review</footer>`,
};

createApp(App).mount('#app');
