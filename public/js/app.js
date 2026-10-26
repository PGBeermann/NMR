/* Analizador RMN — interfaz web */
'use strict';

const $ = (sel) => document.querySelector(sel);
const esc = (s) =>
  String(s ?? '').replace(/[&<>"']/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' })[c]);
const fmt = (v, d = 2) => (v === null || v === undefined || v === '' ? '—' : Number(v).toFixed(d));
const NUC = { '1H': '¹H', '13C': '¹³C' };

let jsmeApplet = null;
let pendingStructure = null;

// JSME llama a esta función cuando termina de cargarse
window.jsmeOnLoad = function () {
  const el = $('#jsme_container');
  el.innerHTML = '';
  const w = Math.max(320, el.parentElement.clientWidth - 2);
  jsmeApplet = new JSApplet.JSME('jsme_container', `${w}px`, '400px', { options: 'oldlook,star,newlook' });
  if (pendingStructure) loadIntoEditor(pendingStructure);
};

function loadIntoEditor({ smiles, molfile }) {
  if (!jsmeApplet) {
    pendingStructure = { smiles, molfile };
    return;
  }
  try {
    if (molfile) jsmeApplet.readMolFile(molfile);
    else if (smiles) jsmeApplet.readGenericMolecularInput(smiles);
  } catch (e) {
    console.warn('JSME', e);
  }
}

async function api(url, opts = {}) {
  const res = await fetch(url, {
    ...opts,
    headers: { 'Content-Type': 'application/json', ...(opts.headers || {}) },
  });
  const data = await res.json().catch(() => ({}));
  if (!res.ok) throw new Error(data.error || `Error ${res.status}`);
  return data;
}

function showLoading(text) {
  $('#loadingText').textContent = text || 'Analizando…';
  $('#loading').classList.remove('hidden');
  $('#error').classList.add('hidden');
}
function hideLoading() {
  $('#loading').classList.add('hidden');
}
function showError(msg) {
  hideLoading();
  const e = $('#error');
  e.textContent = msg;
  e.classList.remove('hidden');
}

const options = () => ({ frequency: Number($('#frequency').value), solvent: $('#solvent').value });

// ---------- Pestañas ----------
document.querySelectorAll('.tab').forEach((b) =>
  b.addEventListener('click', () => {
    document.querySelectorAll('.tab').forEach((t) => t.classList.toggle('active', t === b));
    document.querySelectorAll('.tabpanel').forEach((p) => p.classList.add('hidden'));
    $(`#tab-${b.dataset.tab}`).classList.remove('hidden');
    if (b.dataset.tab === 'saved') {
      loadSaved();
      loadUploads();
    }
  }),
);

// ---------- Estado ----------
async function refreshStatus() {
  try {
    const s = await api('/api/status');
    const n = s.nmrshiftdb;
    let txt = `NMRShiftDB2: ${n.state}`;
    if (n.state === 'lista') txt = `NMRShiftDB2: ${n.records.toLocaleString('es')} compuestos`;
    if (n.state === 'construyendo' && n.progress) {
      const p = n.progress;
      txt =
        p.phase === 'descargando'
          ? `NMRShiftDB2: descargando ${(p.bytes / 1e6).toFixed(0)}${p.total ? ` / ${(p.total / 1e6).toFixed(0)}` : ''} MB`
          : `NMRShiftDB2: indexando ${p.records?.toLocaleString('es') || ''}`;
    }
    const pred = n.predictor === 'listo' ? '' : n.predictor === 'construyendo' ? ' · predictor local: construyendo…' : '';
    $('#status').textContent = `${txt}${pred} · ${s.store.compounds} compuestos guardados`;
    const sel = $('#solvent');
    if (!sel.options.length) {
      for (const [k, v] of Object.entries(s.solvents)) sel.add(new Option(v.name, k, k === 'CDCl3', k === 'CDCl3'));
    }
    if (n.state === 'construyendo' || n.predictor === 'construyendo') setTimeout(refreshStatus, 5000);
  } catch {
    $('#status').textContent = 'Servidor no disponible';
  }
}

// ---------- Modo 1: dibujo ----------
$('#examples').addEventListener('change', (e) => {
  if (!e.target.value) return;
  $('#smilesInput').value = e.target.value;
  loadIntoEditor({ smiles: e.target.value });
});
$('#smilesInput').addEventListener('change', (e) => e.target.value && loadIntoEditor({ smiles: e.target.value }));
$('#btnClear').addEventListener('click', () => {
  jsmeApplet?.reset();
  $('#smilesInput').value = '';
  $('#examples').value = '';
});
$('#btnAnalyze').addEventListener('click', () => {
  let smiles = '';
  let molfile = '';
  if (jsmeApplet) {
    smiles = jsmeApplet.smiles();
    molfile = jsmeApplet.molFile();
  }
  if (!smiles) smiles = $('#smilesInput').value.trim();
  if (!smiles && !molfile) return showError('Dibuje una estructura o escriba un SMILES.');
  runAnalysis({ smiles, molfile: smiles ? molfile : '', mode: 'estructura', query: smiles });
});

/**
 * Análisis completo de una estructura. Por defecto se muestra en #results; con `target`
 * se muestra en otra zona (p. ej. debajo del espectro experimental cargado, sin borrarlo).
 */
async function runAnalysis(payload, { target = '#results', header = '' } = {}) {
  const box = $(target);
  const inline = target !== '#results';
  if (inline) {
    $('#error').classList.add('hidden');
    box.innerHTML = '<div class="loading"><div class="spinner"></div><span>Prediciendo los espectros del candidato y buscando en las bases de datos…</span></div>';
    box.scrollIntoView({ behavior: 'smooth', block: 'start' });
  } else {
    showLoading('Prediciendo espectros y buscando en las bases de datos…');
    box.innerHTML = '';
  }
  try {
    const r = await api('/api/analyze', { method: 'POST', body: JSON.stringify({ ...payload, ...options() }) });
    if (!inline) hideLoading();
    renderResults(r, target, header);
    if (payload.mode !== 'estructura') loadIntoEditor({ molfile: r.compound.molfile });
    refreshStatus();
    box.scrollIntoView({ behavior: 'smooth', block: 'start' });
  } catch (e) {
    if (inline) box.innerHTML = `<div class="error">${esc(e.message)}</div>`;
    else showError(e.message);
  }
}

// ---------- Modo 2: fórmula ----------
let formulaData = null;
let formulaPage = 0;
const PAGE = 25;

$('#formFormula').addEventListener('submit', async (e) => {
  e.preventDefault();
  showLoading('Buscando compuestos con esa fórmula en PubChem y NMRShiftDB2…');
  $('#formulaResults').innerHTML = '';
  try {
    formulaData = await api(`/api/search/formula?formula=${encodeURIComponent($('#formulaInput').value)}`);
    formulaPage = 0;
    hideLoading();
    renderFormula();
  } catch (err) {
    showError(err.message);
  }
});

function renderFormula(filter = $('#formulaFilter')?.value || '') {
  const d = formulaData;
  const f = filter.toLowerCase();
  const items = d.items.filter((i) => !f || `${i.name} ${i.iupac} ${i.smiles}`.toLowerCase().includes(f));
  const pages = Math.max(1, Math.ceil(items.length / PAGE));
  formulaPage = Math.min(formulaPage, pages - 1);
  const slice = items.slice(formulaPage * PAGE, formulaPage * PAGE + PAGE);
  $('#formulaResults').innerHTML = `
    <p><b>${d.total}</b> compuestos con fórmula <b>${esc(d.formula)}</b>
      (${d.withSpectra} con espectros asignados en NMRShiftDB2)${d.truncated ? ' · se muestran los primeros 500 de PubChem' : ''}.</p>
    <div class="searchbar"><input id="formulaFilter" placeholder="Filtrar resultados por nombre" value="${esc(filter)}" /></div>
    <div class="table-wrap"><table>
      <thead><tr><th></th><th>Nombre</th><th>Nombre IUPAC</th><th>PM</th><th>Espectros</th><th></th></tr></thead>
      <tbody>${slice
        .map(
          (i, k) => `<tr>
          <td>${i.cid ? `<img class="pc-img" loading="lazy" alt="" src="https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${i.cid}/PNG?image_size=100x100">` : ''}</td>
          <td><b>${esc(i.name)}</b><br><span class="muted">${i.cid ? `CID ${i.cid}` : ''}</span></td>
          <td>${esc(i.iupac || '—')}<br><code class="muted">${esc(i.smiles || '')}</code></td>
          <td class="num">${fmt(i.mw, 2)}</td>
          <td>${i.nmrshiftdb.length ? [...new Set(i.nmrshiftdb.flatMap((n) => n.nuclei))].map((n) => `<span class="badge exp">${NUC[n] || n} exp.</span>`).join('') : '<span class="muted">solo predicción</span>'}</td>
          <td><button class="small primary" data-k="${formulaPage * PAGE + k}">Analizar</button></td>
        </tr>`,
        )
        .join('')}</tbody>
    </table></div>
    <div class="pager">
      <button class="small" id="pgPrev" ${formulaPage === 0 ? 'disabled' : ''}>‹ Anterior</button>
      <span>Página ${formulaPage + 1} de ${pages}</span>
      <button class="small" id="pgNext" ${formulaPage >= pages - 1 ? 'disabled' : ''}>Siguiente ›</button>
    </div>`;
  $('#formulaFilter').addEventListener('input', (e) => {
    formulaPage = 0;
    renderFormula(e.target.value);
    const inp = $('#formulaFilter');
    inp.focus();
    inp.setSelectionRange(inp.value.length, inp.value.length);
  });
  $('#pgPrev').onclick = () => {
    formulaPage--;
    renderFormula(filter);
  };
  $('#pgNext').onclick = () => {
    formulaPage++;
    renderFormula(filter);
  };
  $('#formulaResults')
    .querySelectorAll('button[data-k]')
    .forEach((b) =>
      b.addEventListener('click', () => {
        const i = items[Number(b.dataset.k)];
        runAnalysis({ smiles: i.smiles, cid: i.cid, name: i.name, mode: 'formula', query: d.formula });
      }),
    );
}

// ---------- Modo 3: nombre ----------
$('#formName').addEventListener('submit', async (e) => {
  e.preventDefault();
  const q = $('#nameInput').value.trim();
  showLoading('Resolviendo el nombre…');
  $('#nameResults').innerHTML = '';
  try {
    const d = await api(`/api/search/name?q=${encodeURIComponent(q)}`);
    hideLoading();
    if (d.items.length === 1) {
      const i = d.items[0];
      $('#nameResults').innerHTML = `<p class="muted">Encontrado en ${esc(d.source)}: <b>${esc(i.name)}</b> ${i.cid ? `(CID ${i.cid})` : ''}</p>`;
      return runAnalysis({ smiles: i.smiles, cid: i.cid, name: i.cid ? i.name : null, mode: 'nombre', query: q });
    }
    $('#nameResults').innerHTML = `<p>${d.items.length} coincidencias en ${esc(d.source)}:</p>
      <div class="table-wrap"><table><thead><tr><th>Nombre</th><th>IUPAC</th><th>Fórmula</th><th></th></tr></thead><tbody>
      ${d.items.map((i, k) => `<tr><td>${esc(i.name)}<br><span class="muted">${i.cid ? `CID ${i.cid}` : ''}</span></td><td>${esc(i.iupac || '')}</td><td>${esc(i.formula || '')}</td><td><button class="small primary" data-k="${k}">Analizar</button></td></tr>`).join('')}
      </tbody></table></div>`;
    $('#nameResults')
      .querySelectorAll('button[data-k]')
      .forEach((b) =>
        b.addEventListener('click', () => {
          const i = d.items[Number(b.dataset.k)];
          runAnalysis({ smiles: i.smiles, cid: i.cid, name: i.name, mode: 'nombre', query: q });
        }),
      );
  } catch (err) {
    showError(err.message);
  }
});

// ---------- Guardados ----------
async function loadSaved() {
  const q = $('#savedFilter').value.trim();
  try {
    const d = await api(`/api/compounds?q=${encodeURIComponent(q)}&limit=200`);
    if (!d.items.length) {
      $('#savedList').innerHTML = '<p class="muted">Aún no hay espectros guardados. Cada análisis se guarda automáticamente.</p>';
      return;
    }
    $('#savedList').innerHTML = `<p class="muted">${d.total} compuestos guardados.</p>
      <div class="table-wrap"><table><thead><tr><th>Compuesto</th><th>Fórmula</th><th>Espectros</th><th>Actualizado</th><th></th></tr></thead><tbody>
      ${d.items
        .map(
          (c) => `<tr><td><b>${esc(c.name || c.iupac || c.smiles)}</b><br><span class="muted">${esc(c.iupac || '')}</span></td>
        <td>${esc(c.formula)}</td>
        <td><span class="badge">${c.n_spectra} total</span>${c.n_exp ? `<span class="badge exp">${c.n_exp} exp.</span>` : ''}</td>
        <td class="num">${new Date(c.updated_at).toLocaleString('es')}</td>
        <td class="num"><button class="small primary" data-view="${c.id}">Ver</button>
          <a class="small" href="/api/compounds/${c.id}/export">JSON</a>
          <button class="small danger" data-del="${c.id}">Eliminar</button></td></tr>`,
        )
        .join('')}
      </tbody></table></div>`;
    $('#savedList').querySelectorAll('[data-view]').forEach((b) =>
      b.addEventListener('click', async () => {
        showLoading('Cargando…');
        try {
          const r = await api(`/api/compounds/${b.dataset.view}`);
          hideLoading();
          renderResults({ ...r, summary: null, warnings: [] });
          $('#results').scrollIntoView({ behavior: 'smooth' });
        } catch (e) {
          showError(e.message);
        }
      }),
    );
    $('#savedList').querySelectorAll('[data-del]').forEach((b) =>
      b.addEventListener('click', async () => {
        if (!confirm('¿Eliminar este compuesto y todos sus espectros guardados?')) return;
        await api(`/api/compounds/${b.dataset.del}`, { method: 'DELETE' });
        loadSaved();
        refreshStatus();
      }),
    );
  } catch (e) {
    showError(e.message);
  }
}
$('#btnSavedRefresh').addEventListener('click', loadSaved);
$('#savedFilter').addEventListener('input', () => {
  clearTimeout(loadSaved.t);
  loadSaved.t = setTimeout(loadSaved, 300);
});

// ---------- Resultados ----------
function renderResults(r, target = '#results', header = '') {
  const c = r.compound;
  const info = c.info || {};
  const s = r.summary;
  const parts = [];
  parts.push(`<section class="panel compound">
    <div class="mol">${c.svg || ''}</div>
    <div>
      <h2>${esc(c.name || c.iupac || c.formula)}</h2>
      ${c.iupac && c.iupac !== c.name ? `<p class="muted" style="margin:0 0 8px">IUPAC: ${esc(c.iupac)}</p>` : ''}
      <dl class="kv">
        <dt>Fórmula</dt><dd>${esc(c.formula)}</dd>
        <dt>Masa molar</dt><dd>${fmt(c.mw, 3)} g/mol (exacta ${fmt(c.exactMass, 4)})</dd>
        <dt>Insaturaciones (IDH)</dt><dd>${info.unsaturation ?? '—'}</dd>
        <dt>SMILES</dt><dd><code>${esc(c.smiles)}</code></dd>
        <dt>InChIKey</dt><dd>${esc(c.inchikey || '—')}</dd>
        <dt>PubChem</dt><dd>${c.cid ? `<a href="${c.pubchemUrl}" target="_blank" rel="noopener">CID ${c.cid}</a>` : 'no encontrado'}</dd>
        ${info.nmrshiftdbIds?.length ? `<dt>NMRShiftDB2</dt><dd>${info.nmrshiftdbIds.map((id) => `<a target="_blank" rel="noopener" href="https://nmrshiftdb.nmr.uni-koeln.de/portal/js_pane/P-Results/nmrshiftdbaction/showDetailsFromHome/molNumber/${id}">#${id}</a>`).join(', ')}</dd>` : ''}
      </dl>
      ${info.synonyms?.length ? `<details><summary>Sinónimos</summary><p class="explain">${info.synonyms.map(esc).join(' · ')}</p></details>` : ''}
      ${
        s
          ? `<div class="chips">${s.functionalGroups.map((g) => `<span class="badge ok">${esc(g)}</span>`).join('')}</div>
        <ul class="summary-list">
          ${s.notes.map((n) => `<li>${esc(n)}</li>`).join('')}
          ${s.protonEnvironments !== null ? `<li>${s.protonEnvironments} tipos de H químicamente distintos (${s.nH} H en total) y ${s.carbonEnvironments} señales de ¹³C (${s.nC} C).</li>` : ''}
          ${s.symmetry ? `<li>${esc(s.symmetry)}</li>` : ''}
          <li>${s.experimentalCount} espectro(s) encontrados en bases de datos.</li>
        </ul>`
          : ''
      }
      <p class="muted" style="font-size:.82rem">Los átomos se numeran como en la estructura. Todo el análisis se guardó automáticamente (pestaña «Espectros guardados»).
      ${c.id ? ` <a href="/api/compounds/${c.id}/export">Descargar análisis completo (JSON)</a>` : ''}</p>
    </div>
  </section>`);
  for (const w of r.warnings || []) parts.push(`<div class="warn">${esc(w)}</div>`);

  const predicted = r.spectra.filter((x) => x.kind === 'predicho');
  const others = r.spectra.filter((x) => x.kind !== 'predicho');
  if (predicted.length) parts.push('<h2 class="section-title">Espectros predichos y análisis</h2>');
  for (const sp of predicted) parts.push(renderSpectrum(sp));
  // Espectros experimentales de las bases de datos: no se muestran hasta que el usuario los pide,
  // y se ofrecen por separado para ¹H y ¹³C
  const uid = `g${++expGroupSeq}`;
  const groups = { '1H': others.filter((x) => x.nucleus === '1H'), '13C': others.filter((x) => x.nucleus === '13C') };
  expGroups.set(uid, groups);
  parts.push(`<h2 class="section-title">Espectros experimentales en bases de datos (${others.length})</h2>`);
  if (!others.length) parts.push('<p class="muted">No se encontraron espectros experimentales para este compuesto en PubChem ni en NMRShiftDB2.</p>');
  else {
    const btn = (nuc) =>
      groups[nuc].length
        ? `<button class="exp-toggle" data-uid="${uid}" data-nuc="${nuc}" aria-expanded="false">Mostrar espectros experimentales de ${NUC[nuc]} RMN (${groups[nuc].length})</button>`
        : `<button disabled title="No hay espectros de ${NUC[nuc]} en las bases de datos">${NUC[nuc]} RMN: sin espectros</button>`;
    parts.push(`<div class="btns exp-buttons">${btn('1H')}${btn('13C')}</div>
      <p class="hint">Los espectros experimentales (PubChem/HMDB, SpectraBase, NMRShiftDB2) se muestran sólo al pulsar el botón del núcleo correspondiente.</p>
      <div class="exp-group hidden" id="exp-${uid}-1H"></div>
      <div class="exp-group hidden" id="exp-${uid}-13C"></div>`);
  }
  if (info.pubchemLinks?.length)
    parts.push(
      `<p class="muted">Otras referencias: ${info.pubchemLinks
        .map((l) => (l.url ? `<a href="${esc(l.url)}" target="_blank" rel="noopener">${esc(l.text)}</a>` : esc(l.text)))
        .join(' · ')}</p>`,
    );
  $(target).innerHTML = header + parts.join('\n');
  bindSpectrumControls();
  bindExperimentalToggles($(target));
}

const expGroups = new Map();
let expGroupSeq = 0;

function bindExperimentalToggles(root) {
  root.querySelectorAll('.exp-toggle').forEach((b) =>
    b.addEventListener('click', () => {
      const { uid, nuc } = b.dataset;
      const box = document.getElementById(`exp-${uid}-${nuc}`);
      const open = box.classList.contains('hidden');
      if (open && !box.dataset.rendered) {
        // Se dibujan (y se piden las imágenes) sólo la primera vez que se solicitan
        const list = expGroups.get(uid)?.[nuc] || [];
        box.innerHTML = `<h3 class="section-title" style="font-size:1rem">Espectros experimentales de ${NUC[nuc]} RMN (${list.length})</h3>${list.map((sp) => renderSpectrum(sp)).join('\n')}`;
        box.dataset.rendered = '1';
        bindSpectrumControls();
      }
      box.classList.toggle('hidden', !open);
      b.setAttribute('aria-expanded', String(open));
      b.classList.toggle('primary', open);
      const n = expGroups.get(uid)?.[nuc]?.length ?? 0;
      b.textContent = `${open ? 'Ocultar' : 'Mostrar'} espectros experimentales de ${NUC[nuc]} RMN (${n})`;
      if (open) box.scrollIntoView({ behavior: 'smooth', block: 'start' });
    }),
  );
}

/** Letra de señal resaltada. */
function letterBadge(l) {
  return l ? `<span class="sig-letter">${esc(l)}</span>` : '';
}

/**
 * Estructura con las letras de las señales de ¹H junto a los átomos que portan esos H,
 * acompañada de la leyenda letra → desplazamiento.
 */
function assignmentBlock(svg, signals, title = 'Asignación de los hidrógenos') {
  if (!svg || !signals?.some((g) => g.label)) return '';
  const legend = signals
    .filter((g) => g.label)
    .map(
      (g) =>
        `<li>${letterBadge(g.label)} <b>δ ${fmt(g.delta, 2)}</b> · ${esc(g.multiplicity || '')}${g.nH ? ` · ${g.nH}H` : ''}${
          g.environment ? ` <span class="muted">— ${esc(g.environment)}</span>` : g.note ? ` <span class="muted">— ${esc(g.note)}</span>` : ''
        }</li>`,
    )
    .join('');
  return `<div class="assign">
    <div class="assign-mol">${svg}</div>
    <div><h4 style="margin:0 0 6px">${esc(title)}</h4>
      <p class="hint" style="margin:0 0 6px">Cada letra señala el átomo que porta esos hidrógenos y el pico correspondiente del espectro.</p>
      <ul class="assign-legend">${legend}</ul></div>
  </div>`;
}

function downloads(sp) {
  const b = `/api/spectra/${sp.id}`;
  return `<div class="downloads">
    <a href="${b}/image.png?download=1">Imagen PNG</a>
    <a href="${b}/image.svg?download=1">Imagen SVG</a>
    <a href="${b}/export?format=jdx">JCAMP-DX</a>
    <a href="${b}/export?format=jdx-peaks">JCAMP picos</a>
    <a href="${b}/export?format=csv">CSV (señales + intensidades)</a>
    <a href="${b}/export?format=xy">CSV x/y</a>
    <a href="${b}/export?format=json">JSON</a>
  </div>`;
}

function zoomControls(sp) {
  return `<div class="zoom" data-id="${sp.id}">
    Zoom (ppm): <input type="number" step="0.1" class="z-from" placeholder="desde"> –
    <input type="number" step="0.1" class="z-to" placeholder="hasta">
    <button class="small z-apply">Aplicar</button><button class="small z-reset">Espectro completo</button>
    <span class="muted">También puede pulsar una fila de la tabla para ampliar esa señal.</span>
  </div>`;
}

function comparisonTable(cmp) {
  if (!cmp) return '';
  return `<details><summary>Comparación con la predicción: error medio ${fmt(cmp.mae, 3)} ppm (máx. ${fmt(cmp.max, 3)}), ${cmp.n} señales, ${esc(cmp.method || '')}</summary>
    <div class="table-wrap"><table><thead><tr><th>Átomo</th><th>δ experimental</th><th>δ predicho</th><th>Diferencia</th></tr></thead><tbody>
    ${cmp.rows.map((x) => `<tr><td>${esc(x.atom)}</td><td class="num">${fmt(x.experimental, 2)}</td><td class="num">${fmt(x.predicted, 2)}</td><td class="num">${x.diff > 0 ? '+' : ''}${fmt(x.diff, 2)}</td></tr>`).join('')}
    </tbody></table></div></details>`;
}

function renderSpectrum(sp) {
  const kindBadge = sp.kind === 'predicho' ? '<span class="badge">predicho</span>' : sp.kind === 'calculado' ? '<span class="badge calc">calculado</span>' : '<span class="badge exp">experimental</span>';
  const meta = [sp.frequency ? `${sp.frequency} MHz` : null, sp.solvent, sp.meta?.temperature && sp.meta.temperature !== 'Unreported' ? `${sp.meta.temperature} K` : null]
    .filter(Boolean)
    .join(' · ');
  const hasImage = sp.kind === 'predicho' || (sp.peaks && sp.peaks.length);
  let body = '';
  if (hasImage) {
    body += `<img class="spec-img" id="img-${sp.id}" src="${sp.image}" alt="Espectro ${esc(sp.nucleus)}">`;
    body += zoomControls(sp);
  } else if (sp.meta?.thumbnail) {
    body += `<p class="muted">Esta fuente sólo publica la imagen del espectro (sin lista de picos).</p>
      <img class="thumb" loading="lazy" src="${esc(sp.meta.thumbnail)}" alt="Espectro de referencia">`;
  }
  if (sp.meta?.url) body += `<p class="muted" style="font-size:.85rem">Fuente: <a href="${esc(sp.meta.url)}" target="_blank" rel="noopener">${esc(sp.meta.url)}</a></p>`;
  if (sp.meta?.note) body += `<p class="muted" style="font-size:.85rem">${esc(sp.meta.note)}</p>`;
  if (sp.meta?.method) body += `<p class="muted" style="font-size:.85rem">${esc(sp.meta.method)}</p>`;
  if (sp.meta?.fields) {
    const f = Object.entries(sp.meta.fields).filter(([k]) => !/Shifts/i.test(k));
    if (f.length) body += `<dl class="kv">${f.map(([k, v]) => `<dt>${esc(k)}</dt><dd>${esc(v)}</dd>`).join('')}</dl>`;
  }
  if (sp.analysis?.literature) body += `<p><b>Formato de publicación:</b></p><div class="literature">${esc(sp.analysis.literature)}</div>`;

  if (sp.kind === 'predicho' && sp.nucleus === '1H') {
    body += assignmentBlock(sp.analysis.assignmentSvg, sp.signals);
    body += `<p><b>Relación de integrales:</b> ${esc(sp.analysis.integrationRatio)} (total ${sp.analysis.totalH} H, ${sp.analysis.environments} tipos de H)</p>
    <div class="table-wrap"><table><thead><tr><th>Señal</th><th>δ (ppm)</th><th>Multiplicidad</th><th>J (Hz)</th><th>Integración</th><th>Átomo nº</th><th>Entorno / interpretación</th></tr></thead><tbody>
    ${sp.signals
      .map(
        (g) => `<tr class="clickable" data-zoom="${sp.id}" data-ppm="${g.delta}">
      <td>${letterBadge(g.label)}</td>
      <td class="num"><b>${fmt(g.delta, 2)}</b></td>
      <td>${esc(g.multiplicity)}<br><span class="muted">${esc(g.multiplicityName)}</span></td>
      <td class="num">${g.js.length ? g.js.map((j) => `${fmt(j.J, 1)}${j.bonds ? ` <span class="muted">(${j.bonds}J)</span>` : ''}`).join('<br>') : '—'}</td>
      <td class="num">${g.nH}H</td>
      <td>${g.atoms.join(', ')}</td>
      <td class="explain">${g.exchangeable ? '<span class="badge calc">intercambiable · estimado</span> ' : ''}${esc(g.explanation)}<br><span class="muted">${esc(g.region)}</span></td></tr>`,
      )
      .join('')}
    </tbody></table></div>
    <details><summary>Lista de líneas del espectro simulado (${sp.peaks.length}) con intensidades</summary>
      <div class="table-wrap"><table><thead><tr><th>δ (ppm)</th><th>ν (Hz)</th><th>Intensidad relativa (%)</th><th>Señal</th></tr></thead><tbody>
      ${sp.peaks.map((p) => `<tr><td class="num">${fmt(p.ppm, 4)}</td><td class="num">${fmt(p.hz, 2)}</td><td class="num">${fmt(p.intensity, 1)}</td><td class="num">${p.label ? letterBadge(p.label) + ' ' : ''}${fmt(p.signal, 2)}</td></tr>`).join('')}
      </tbody></table></div></details>`;
  } else if (sp.kind === 'predicho') {
    const d = sp.analysis.dept;
    body += `<p><b>DEPT:</b> ${d.CH3} CH₃ · ${d.CH2} CH₂ · ${d.CH} CH · ${d.C} C cuaternarios (total ${sp.analysis.totalC} C en ${sp.analysis.environments} señales)</p>
    <div class="table-wrap"><table><thead><tr><th>δ (ppm)</th><th>Tipo</th><th>DEPT-135</th><th>DEPT-90</th><th>Intensidad</th><th>Átomos</th><th>Entorno / interpretación</th></tr></thead><tbody>
    ${sp.signals
      .map(
        (g) => `<tr class="clickable" data-zoom="${sp.id}" data-ppm="${g.delta}" data-width="6">
      <td class="num"><b>${fmt(g.delta, 1)}</b></td><td>${esc(g.type)}</td><td>${esc(g.dept135)}</td><td>${esc(g.dept90)}</td>
      <td class="num">${fmt(g.intensity, 0)}%</td><td>${g.atoms.map((a) => `C${a}`).join(', ')}</td>
      <td class="explain">${esc(g.explanation)}<br><span class="muted">${esc(g.region)}</span></td></tr>`,
      )
      .join('')}
    </tbody></table></div>`;
  } else if (sp.signals?.length) {
    body += `<details ${sp.signals.length <= 30 ? 'open' : ''}><summary>Señales / picos (${sp.signals.length})</summary>
      <div class="table-wrap"><table><thead><tr><th>δ (ppm)</th><th>Tipo / mult.</th><th>Intensidad</th><th>Átomos</th></tr></thead><tbody>
      ${sp.signals
        .map(
          (g) => `<tr class="clickable" data-zoom="${sp.id}" data-ppm="${g.delta}" data-width="${sp.nucleus === '1H' ? 0.4 : 6}"><td class="num">${fmt(g.delta, 2)}</td><td>${esc(g.type || g.multiplicity || '—')}</td>
          <td class="num">${g.intensity !== undefined ? fmt(g.intensity, 1) : g.count ? `${g.count} át.` : '—'}</td><td>${(g.atoms || []).map((a) => `${a}`).join(', ') || '—'}</td></tr>`,
        )
        .join('')}
      </tbody></table></div></details>`;
    body += comparisonTable(sp.analysis?.comparison);
  }
  body += downloads(sp);
  return `<section class="panel spectrum ${sp.kind}">
    <h3>${NUC[sp.nucleus] || sp.nucleus} RMN ${kindBadge}</h3>
    <div class="meta">${esc(sp.source)}${meta ? ` · ${esc(meta)}` : ''}</div>
    ${body}
  </section>`;
}

function setZoom(id, from, to) {
  const img = $(`#img-${id}`);
  if (!img) return;
  const box = document.querySelector(`.zoom[data-id="${id}"]`);
  if (from === undefined) {
    img.src = `/api/spectra/${id}/image.svg`;
    box.querySelector('.z-from').value = '';
    box.querySelector('.z-to').value = '';
    return;
  }
  img.src = `/api/spectra/${id}/image.svg?from=${from}&to=${to}`;
  box.querySelector('.z-from').value = from;
  box.querySelector('.z-to').value = to;
}

function bindSpectrumControls() {
  document.querySelectorAll('.zoom[data-id]').forEach((z) => {
    const id = z.dataset.id;
    z.querySelector('.z-apply').onclick = () => {
      const a = parseFloat(z.querySelector('.z-from').value);
      const b = parseFloat(z.querySelector('.z-to').value);
      if (Number.isFinite(a) && Number.isFinite(b) && a !== b) setZoom(id, Math.min(a, b), Math.max(a, b));
    };
    z.querySelector('.z-reset').onclick = () => setZoom(id);
  });
  // Cada fila se enlaza una sola vez (la función se vuelve a llamar al desplegar espectros)
  document.querySelectorAll('tr[data-zoom]:not([data-bound])').forEach((tr) => {
    tr.dataset.bound = '1';
    tr.addEventListener('click', () => {
      const ppm = Number(tr.dataset.ppm);
      const w = Number(tr.dataset.width) || 0.25;
      setZoom(tr.dataset.zoom, +(ppm - w).toFixed(3), +(ppm + w).toFixed(3));
      $(`#img-${tr.dataset.zoom}`)?.scrollIntoView({ behavior: 'smooth', block: 'center' });
    });
  });
}

refreshStatus();
