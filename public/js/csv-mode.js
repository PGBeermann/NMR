/* Modo 4: espectro experimental en CSV → estructuras candidatas con % de concordancia */
'use strict';

let currentUpload = null; // { id, csv }

$('#csvFile').addEventListener('change', async (e) => {
  const f = e.target.files[0];
  if (!f) return;
  if (f.size > 5e6) return showError('El archivo es demasiado grande (máx. 5 MB).');
  $('#csvText').value = await f.text();
  if (!$('#csvName').value) $('#csvName').value = f.name.replace(/\.[^.]+$/, '');
});

$('#btnCsvExample').addEventListener('click', async () => {
  $('#csvText').value = await (await fetch('examples/fenol-cdcl3-400MHz.tsv')).text();
  $('#csvName').value = 'Ejemplo (aromático + OH, CDCl₃, 400 MHz)';
  $('#csvFormula').value = 'C6H6O';
});

$('#btnCsv').addEventListener('click', () => analyzeCsv());

// Al elegir ¹³C se propone una DE máxima acorde (1.5 ppm)
$('#csvNucleus').addEventListener('change', (e) => {
  if (e.target.value === '13C' && Number($('#csvSdMax').value) < 0.5) $('#csvSdMax').value = 1.5;
  if (e.target.value === '1H' && Number($('#csvSdMax').value) > 0.5) $('#csvSdMax').value = 0.1;
});

async function analyzeCsv({ integrals = null, uploadId = null } = {}) {
  const csv = currentUpload && uploadId ? currentUpload.csv : $('#csvText').value;
  if (!csv.trim()) return showError('Cargue un archivo o pegue la tabla de picos.');
  showLoading('Analizando el espectro y buscando estructuras compatibles (si indicó la fórmula puede tardar hasta un minuto)…');
  $('#results').innerHTML = '';
  $('#uploadView').innerHTML = '';
  try {
    const r = await api('/api/spectrum/analyze', {
      method: 'POST',
      body: JSON.stringify({
        csv,
        name: $('#csvName').value.trim(),
        nucleus: $('#csvNucleus').value,
        formula: $('#csvFormula').value.trim(),
        excludeSolvent: $('#csvExcludeSolvent').checked,
        rangeOnly: $('#csvRangeOnly').checked,
        sdMax: Number($('#csvSdMax').value) || null,
        minScore: $('#csvMinScore').value === '' ? null : Number($('#csvMinScore').value),
        integrals,
        uploadId,
        ...options(),
      }),
    });
    hideLoading();
    currentUpload = { id: r.id, csv };
    renderUploadResult(r);
    refreshStatus();
    $('#uploadView').scrollIntoView({ behavior: 'smooth', block: 'start' });
  } catch (e) {
    showError(e.message);
  }
}

function scoreColor(p) {
  return p >= 85 ? '#2e7d4f' : p >= 70 ? '#5f8f1f' : p >= 55 ? '#b07d00' : '#b02a37';
}

function signalRow(u, s, i, isH) {
  if (!isH) {
    return `<tr><td class="num"><b>${fmt(s.delta, 2)}</b></td><td class="num">${fmt(s.intensity, 1)}</td>
      <td class="explain">${esc(s.region)}${s.solvent ? ' <span class="badge calc">disolvente</span>' : ''}</td></tr>`;
  }
  const J = s.J.length ? s.J.join(', ') : s.apparentJ?.length ? `<span class="muted">aparentes: ${s.apparentJ.join(', ')}</span>` : '—';
  return `<tr class="clickable" data-uzoom="${u.id}" data-from="${s.from}" data-to="${s.to}">
    <td>${letterBadge(s.label)}</td>
    <td class="num"><b>${fmt(s.delta, 3)}</b></td>
    <td class="num">${fmt(s.to, 3)}–${fmt(s.from, 3)}</td>
    <td>${esc(s.multiplicity)}<br><span class="muted">${esc(s.multiplicityName)}</span></td>
    <td class="num">${J}</td>
    <td class="num">${s.nLines}</td>
    <td class="num">${fmt((s.fraction ?? 0) * 100, 1)}</td>
    <td><input class="nh" type="number" min="0.5" step="0.5" value="${s.nH}" data-i="${i}"></td>
    <td class="explain">${esc(s.region)}${s.solvent ? ' <span class="badge calc">disolvente</span>' : ''}</td></tr>`;
}

function candidateCard(x, k, isH) {
  const badges = x.sources
    .map((s) => `<span class="badge ${/predic/i.test(s) ? '' : /calcul/i.test(s) ? 'calc' : 'exp'}">${esc(s)}</span>`)
    .join(' ');
  const rows = x.matches
    .map((m) => {
      const cand = m.candidate === null ? '—' : `${fmt(m.candidate, 2)}${m.n ? ` <span class="muted">(${m.n}${isH ? 'H' : 'C'})</span>` : ''}`;
      const diff =
        m.diff === null
          ? `<span class="muted">${esc(m.note || '')}</span>`
          : `${m.diff > 0 ? '+' : ''}${fmt(m.diff, 3)}${m.note ? `<br><span class="muted">${esc(m.note)}</span>` : ''}`;
      return `<tr><td>${letterBadge(m.label)}</td><td class="num">${fmt(m.experimental, 2)}</td><td class="num">${cand}</td><td class="num">${diff}</td></tr>`;
    })
    .join('');
  return `<div class="cand">
    <div class="rank">#${k + 1} · ${badges}</div>
    <div class="mol">${x.svg || ''}</div>
    ${x.labeled ? '<div class="hint" style="text-align:center;margin-top:-4px">Las letras indican el multiplete experimental asignado a los H de cada átomo.</div>' : ''}
    <h4>${esc(x.name || x.iupac || x.smiles)}</h4>
    ${x.iupac && x.iupac !== x.name ? `<div class="muted" style="font-size:.8rem">${esc(x.iupac)}</div>` : ''}
    <div class="muted" style="font-size:.8rem">${esc(x.formula || '')}${x.solvent ? ` · espectro en ${esc(x.solvent)}` : ''}</div>
    <div><span class="pct" style="color:${scoreColor(x.score)}">${fmt(x.score, 1)} %</span> de concordancia</div>
    <div class="meter"><span style="width:${Math.max(2, x.score)}%;background:${scoreColor(x.score)}"></span></div>
    ${
      x.errorSd !== null && x.errorSd !== undefined
        ? `<div style="font-size:.85rem"><b>DE del error de δ: ${fmt(x.errorSd, 3)} ppm</b> ${x.sdOk ? '<span class="badge ok">≤ límite</span>' : '<span class="badge exp">> límite</span>'}
           <span class="muted">· error medio ${x.errorMean > 0 ? '+' : ''}${fmt(x.errorMean, 3)} ppm · ${x.nErrors} señales</span></div>`
        : ''
    }
    <div class="breakdown">Desplazamientos ${fmt(x.shift, 0)} % · Cobertura ${fmt(x.coverage, 0)} %${
      x.integration !== null ? ` · Integración ${fmt(x.integration, 0)} %` : ''
    } · Nº de señales ${fmt(x.signalCount, 0)} %${x.multiplicity !== null && x.multiplicity !== undefined ? ` · Multiplicidad ${fmt(x.multiplicity, 0)} %` : ''}${
      x.exchangeable ? ' · incluye OH/NH intercambiable' : ''
    }</div>
    <details><summary>Señales comparadas</summary>
      <table><thead><tr><th>Señal</th><th>δ exp.</th><th>δ candidato</th><th>Error (ppm)</th></tr></thead><tbody>${rows}</tbody></table>
      ${x.unexplained?.length ? `<p class="muted" style="font-size:.8rem">Señales del candidato sin correspondencia: ${x.unexplained.join(', ')} ppm</p>` : ''}
    </details>
    <div class="btns" style="margin-top:auto">
      <button class="small primary" data-cand="${k}">Análisis completo</button>
      ${x.url ? `<a class="small" target="_blank" rel="noopener" href="${esc(x.url)}">NMRShiftDB2</a>` : ''}
      ${x.cid ? `<a class="small" target="_blank" rel="noopener" href="https://pubchem.ncbi.nlm.nih.gov/compound/${x.cid}">PubChem</a>` : ''}
    </div>
  </div>`;
}

function renderUploadResult(u) {
  const a = u.analysis;
  const isH = u.nucleus === '1H';
  const c = u.candidates || { candidates: [] };
  const head = isH
    ? '<th>Señal</th><th>δ (ppm)</th><th>Intervalo</th><th>Multiplicidad</th><th>J (Hz)</th><th>Líneas</th><th>Área (%)</th><th>nH</th><th>Región</th>'
    : '<th>δ (ppm)</th><th>Intensidad por área (%)</th><th>Región</th>';
  const parts = [];
  parts.push(`<section class="panel spectrum experimental">
    <h3>${NUC[u.nucleus]} RMN — ${esc(u.name || 'Espectro cargado')} <span class="badge exp">experimental (CSV)</span>
      <button class="small" id="btnCloseUpload" style="margin-left:auto" title="Quitar este espectro de la pantalla (sigue guardado)">Cerrar espectro ✕</button></h3>
    <div class="meta">${u.frequency} MHz (${esc(a.frequencySource || '')}) · ${esc(u.solvent)} · ${u.peaks.length} picos${
      u.formula ? ` · ${esc(u.formula)} (IDH ${a.unsaturation})` : ''
    }</div>
    <img class="spec-img" id="img-u${u.id}" src="${u.image}" alt="Espectro cargado">
    <div class="zoom" id="zoom-u${u.id}">
      Zoom (ppm): <input type="number" step="0.1" class="z-from" placeholder="desde"> –
      <input type="number" step="0.1" class="z-to" placeholder="hasta">
      <button class="small z-apply">Aplicar</button><button class="small z-reset">Espectro completo</button>
      <span class="muted">Pulse una fila de la tabla para ampliar ese multiplete.</span>
    </div>
    ${a.notes.map((n) => `<div class="warn">${esc(n)}</div>`).join('')}
    <p><b>Formato de publicación:</b></p><div class="literature">${esc(a.literature)}</div>
    <div class="table-wrap"><table><thead><tr>${head}</tr></thead><tbody>
      ${a.signals.map((s, i) => signalRow(u, s, i, isH)).join('')}
    </tbody></table></div>
    ${isH ? '<div class="btns"><button id="btnRecalc">Recalcular candidatos con estos nH</button><span class="hint">Corrija el nº de H de cada multiplete si conoce las integrales reales.</span></div>' : ''}
    ${
      a.interpretation?.length
        ? `<h4 style="margin:14px 0 6px">Interpretación estructural (fragmentos probables)</h4>
      <ul class="summary-list">${a.interpretation.map((x) => `<li><b>${esc(x.fragment)}</b> — <span class="explain">${esc(x.evidence)}</span></li>`).join('')}</ul>`
        : ''
    }
    <div class="downloads">
      <a href="/api/uploads/${u.id}/image.png?download=1">Imagen PNG</a>
      <a href="/api/uploads/${u.id}/image.svg?download=1">Imagen SVG</a>
      <a href="/api/uploads/${u.id}/export?format=csv">CSV (multipletes + candidatos)</a>
      <a href="/api/uploads/${u.id}/export">JSON</a>
    </div>
  </section>`);

  const minScore = c.params?.minScore ?? 90;
  const sdTxt = c.params?.sdMax ?? (isH ? 0.1 : 1.5);
  parts.push(`<h2 class="section-title">Estructuras candidatas: DE del error ≤ ${sdTxt} ppm y concordancia ≥ ${minScore} % (${c.candidates.length})</h2>
    <p class="muted" style="font-size:.85rem">${esc(c.method || '')} Se compararon ${c.searched?.nmrshiftdb ?? 0} compuestos de NMRShiftDB2${
      c.searched?.predicted ? ` y ${c.searched.predicted} espectros predichos de compuestos de PubChem` : ''
    }. ${u.formula ? '' : '<b>Sin fórmula molecular la búsqueda es muy amplia: indíquela para afinar los resultados.</b>'}</p>`);
  for (const w of c.warnings || []) parts.push(`<div class="warn">${esc(w)}</div>`);
  if (!c.candidates.length) {
    const b = c.bestBelow;
    parts.push(
      `<div class="warn">Ninguna estructura cumple a la vez DE del error ≤ ${sdTxt} ppm y concordancia ≥ ${minScore} %.${
        b ? ` La mejor encontrada fue <b>${esc(b.name || b.smiles)}</b> con ${fmt(b.score, 1)} % (DE ${fmt(b.errorSd, 3)} ppm).` : ''
      } Puede indicar la fórmula molecular, corregir los nH, ampliar la DE máxima o bajar la concordancia mínima.</div>`,
    );
  }
  parts.push(`<div class="cands">${c.candidates.map((x, k) => candidateCard(x, k, isH)).join('')}</div>`);
  parts.push('<div id="candAnalysis"></div>');
  $('#uploadView').innerHTML = parts.join('\n');

  // Zoom
  const img = $(`#img-u${u.id}`);
  const box = $(`#zoom-u${u.id}`);
  const setZ = (f, t) => {
    img.src = f === undefined ? u.image : `${u.image}?from=${f}&to=${t}`;
    box.querySelector('.z-from').value = f ?? '';
    box.querySelector('.z-to').value = t ?? '';
  };
  box.querySelector('.z-apply').onclick = () => {
    const f = parseFloat(box.querySelector('.z-from').value);
    const t = parseFloat(box.querySelector('.z-to').value);
    if (Number.isFinite(f) && Number.isFinite(t) && f !== t) setZ(Math.min(f, t), Math.max(f, t));
  };
  box.querySelector('.z-reset').onclick = () => setZ();
  document.querySelectorAll('tr[data-uzoom]').forEach((tr) =>
    tr.addEventListener('click', (ev) => {
      if (ev.target.tagName === 'INPUT') return;
      const f = Number(tr.dataset.from);
      const t = Number(tr.dataset.to);
      const pad = Math.max(0.03, (t - f) * 0.4);
      setZ(+(f - pad).toFixed(3), +(t + pad).toFixed(3));
      img.scrollIntoView({ behavior: 'smooth', block: 'center' });
    }),
  );
  $('#btnRecalc')?.addEventListener('click', () => {
    const ints = [...document.querySelectorAll('input.nh')].map((i) => Number(i.value));
    if (ints.some((v) => !(v > 0))) return showError('Todos los nH deben ser mayores que cero.');
    analyzeCsv({ integrals: ints, uploadId: u.id });
  });
  $('#btnCloseUpload').addEventListener('click', () => {
    $('#uploadView').innerHTML = '';
    currentUpload = null;
  });
  // Al elegir un candidato, su análisis completo aparece debajo: el espectro experimental sigue arriba
  document.querySelectorAll('[data-cand]').forEach((b) =>
    b.addEventListener('click', () => {
      const x = c.candidates[Number(b.dataset.cand)];
      document.querySelectorAll('.cand').forEach((el) => el.classList.remove('selected'));
      b.closest('.cand').classList.add('selected');
      const header = `<div class="panel selected-cand">
        <div><b>Estructura candidata seleccionada:</b> ${esc(x.name || x.iupac || x.smiles)} —
        <span style="color:${scoreColor(x.score)};font-weight:700">${fmt(x.score, 1)} % de concordancia</span>
        con el espectro experimental «${esc(u.name || 'cargado')}» (arriba).</div>
        <div class="btns"><button class="small" id="btnBackSpectrum">↑ Ver espectro experimental</button>
        <button class="small" id="btnCloseCand">Cerrar análisis del candidato</button></div>
      </div>`;
      runAnalysis(
        { smiles: x.smiles, cid: x.cid, name: x.name, mode: 'espectro', query: u.name || 'csv' },
        { target: '#candAnalysis', header },
      ).then(() => {
        $('#btnBackSpectrum')?.addEventListener('click', () => $('#uploadView').scrollIntoView({ behavior: 'smooth', block: 'start' }));
        $('#btnCloseCand')?.addEventListener('click', () => {
          $('#candAnalysis').innerHTML = '';
          document.querySelectorAll('.cand').forEach((el) => el.classList.remove('selected'));
        });
      });
    }),
  );
}

async function loadUploads() {
  try {
    const d = await api('/api/uploads');
    if (!d.items.length) {
      $('#uploadList').innerHTML = '<p class="muted">Aún no hay espectros cargados.</p>';
      return;
    }
    $('#uploadList').innerHTML = `<div class="table-wrap"><table><thead><tr><th>Muestra</th><th>Núcleo</th><th>Fórmula</th><th>Actualizado</th><th></th></tr></thead><tbody>
      ${d.items
        .map(
          (x) => `<tr><td><b>${esc(x.name || `Espectro #${x.id}`)}</b></td><td>${NUC[x.nucleus] || x.nucleus} · ${x.frequency} MHz</td><td>${esc(x.formula || '—')}</td>
        <td class="num">${new Date(x.updated_at).toLocaleString('es')}</td>
        <td class="num"><button class="small primary" data-uview="${x.id}">Ver</button>
        <a class="small" href="/api/uploads/${x.id}/export?format=csv">CSV</a>
        <button class="small danger" data-udel="${x.id}">Eliminar</button></td></tr>`,
        )
        .join('')}</tbody></table></div>`;
    $('#uploadList').querySelectorAll('[data-uview]').forEach((b) =>
      b.addEventListener('click', async () => {
        const u = await api(`/api/uploads/${b.dataset.uview}`);
        currentUpload = { id: u.id, csv: u.raw };
        $('#csvText').value = u.raw || '';
        $('#csvName').value = u.name || '';
        $('#csvFormula').value = u.formula || '';
        $('#results').innerHTML = '';
        renderUploadResult(u);
        $('#uploadView').scrollIntoView({ behavior: 'smooth' });
      }),
    );
    $('#uploadList').querySelectorAll('[data-udel]').forEach((b) =>
      b.addEventListener('click', async () => {
        if (!confirm('¿Eliminar este espectro cargado y su análisis?')) return;
        await api(`/api/uploads/${b.dataset.udel}`, { method: 'DELETE' });
        loadUploads();
      }),
    );
  } catch (e) {
    showError(e.message);
  }
}
