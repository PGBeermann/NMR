import { Router } from 'express';
import * as store from '../db.js';
import * as pubchem from '../services/pubchem.js';
import * as opsin from '../services/opsin.js';
import * as nmrshiftdb from '../services/nmrshiftdb.js';
import { analyze, publicSpectrum, savedAnalysis } from '../services/analysis.js';
import { elementCounts, normalizeFormula, unsaturation } from '../services/molecule.js';
import { SOLVENTS } from '../services/chemistry.js';
import { renderSpectrum, renderUpload } from '../services/render.js';
import { analyzePeaks, detectNucleus, interpret, literature, parsePeakCsv } from '../services/peaklist.js';
import { MIN_SCORE, SD_MAX, findCandidates } from '../services/candidates.js';
import { carbonFrequency } from '../services/spectrum.js';
import { svgToPng } from '../services/plot.js';
import { toCsv, toJcamp, toJcampPeaks, toXyCsv } from '../services/exporters.js';

export const api = Router();

const wrap = (fn) => (req, res, next) => Promise.resolve(fn(req, res, next)).catch(next);
const safeName = (s) => String(s || 'espectro').replace(/[^\w.-]+/g, '_').slice(0, 80);

api.get('/status', (req, res) => {
  res.json({ ok: true, nmrshiftdb: nmrshiftdb.status, store: store.stats(), solvents: SOLVENTS });
});

api.post(
  '/nmrshiftdb/reindex',
  wrap(async (req, res) => {
    res.json(nmrshiftdb.startIndexing({ force: true }));
  }),
);

// Modo 1: estructura dibujada en JSME (o cualquier SMILES/molfile/CID)
api.post(
  '/analyze',
  wrap(async (req, res) => {
    const { smiles, molfile, cid, name, frequency, solvent, mode = 'estructura', query } = req.body || {};
    const result = await analyze({ smiles, molfile, cid, name, frequency, solvent });
    store.logSearch(mode, query || smiles || cid || 'molfile', result.spectra.length);
    res.json(result);
  }),
);

// Modo 2: búsqueda por fórmula molecular (todos los compuestos con sus nombres)
api.get(
  '/search/formula',
  wrap(async (req, res) => {
    const formula = normalizeFormula(req.query.formula);
    const [cids, local] = await Promise.all([
      pubchem.cidsByFormula(formula).catch(() => []),
      Promise.resolve(nmrshiftdb.findByFormula(formula)),
    ]);
    const props = cids.length ? await pubchem.properties(cids) : [];
    const localByKey = new Map();
    for (const e of local) {
      if (!e.inchikey) continue;
      if (!localByKey.has(e.inchikey)) localByKey.set(e.inchikey, []);
      localByKey.get(e.inchikey).push(e);
    }
    const seen = new Set();
    const items = props
      .filter((p) => p.formula === formula)
      .map((p) => {
        const hits = localByKey.get(p.inchikey) || [];
        seen.add(p.inchikey);
        return {
          cid: p.cid,
          name: p.title,
          iupac: p.iupac,
          formula: p.formula,
          mw: p.mw,
          smiles: p.smiles,
          inchikey: p.inchikey,
          nmrshiftdb: hits.map((h) => ({ id: h.id, nuclei: h.nuclei })),
        };
      });
    // Compuestos con espectros en NMRShiftDB2 que no aparecieron en la lista de PubChem
    for (const e of local) {
      if (e.inchikey && seen.has(e.inchikey)) continue;
      if (e.inchikey) seen.add(e.inchikey);
      items.push({
        cid: null,
        name: e.name || `NMRShiftDB2 #${e.id}`,
        iupac: null,
        formula: e.formula,
        smiles: e.smiles,
        inchikey: e.inchikey,
        nmrshiftdb: [{ id: e.id, nuclei: e.nuclei }],
      });
    }
    // Primero los que tienen espectros experimentales indexados
    items.sort((a, b) => (b.nmrshiftdb.length > 0) - (a.nmrshiftdb.length > 0));
    store.logSearch('formula', formula, items.length);
    res.json({
      formula,
      total: items.length,
      pubchemTotal: cids.length,
      truncated: cids.length >= 500,
      withSpectra: items.filter((i) => i.nmrshiftdb.length).length,
      items,
    });
  }),
);

// Modo 3: búsqueda por nombre IUPAC (o común)
api.get(
  '/search/name',
  wrap(async (req, res) => {
    const q = String(req.query.q || '').trim();
    if (!q) return res.status(400).json({ error: 'Escriba un nombre' });
    let items = [];
    let source = 'PubChem';
    const cids = await pubchem.cidsByName(q).catch(() => []);
    if (cids.length) {
      const props = await pubchem.properties(cids.slice(0, 20));
      items = props.map((p) => ({ cid: p.cid, name: p.title, iupac: p.iupac, formula: p.formula, mw: p.mw, smiles: p.smiles, inchikey: p.inchikey }));
    } else {
      const s = await opsin.nameToStructure(q).catch(() => null);
      if (s) {
        source = 'OPSIN (nombre → estructura)';
        items = [{ cid: null, name: q, iupac: q, smiles: s.smiles, inchikey: s.inchikey }];
      }
    }
    store.logSearch('nombre', q, items.length);
    if (!items.length) return res.status(404).json({ error: `No se encontró ningún compuesto llamado «${q}»` });
    res.json({ query: q, source, items });
  }),
);

// Modo 4: espectro experimental en CSV → análisis + estructuras candidatas con % de concordancia
api.post(
  '/spectrum/analyze',
  wrap(async (req, res) => {
    const b = req.body || {};
    const peaks = parsePeakCsv(b.csv);
    if (peaks.length > 20000) throw Object.assign(new Error('La tabla tiene demasiados picos (máx. 20 000).'), { status: 400 });
    const nucleus = b.nucleus === '1H' || b.nucleus === '13C' ? b.nucleus : detectNucleus(peaks);
    const frequencyFromFile = peaks.frequency || null;
    const userFreq = Number(b.frequency) || 400;
    const frequency = nucleus === '1H' ? frequencyFromFile || userFreq : frequencyFromFile || carbonFrequency(userFreq);
    const solvent = SOLVENTS[b.solvent] ? b.solvent : 'CDCl3';
    const formula = b.formula && String(b.formula).trim() ? normalizeFormula(b.formula) : null;
    const totalH = formula ? elementCounts(formula).H || 0 : null;
    const analysis = analyzePeaks(peaks, {
      nucleus,
      frequency,
      solvent,
      excludeSolvent: Boolean(b.excludeSolvent),
      totalH,
      integrals: Array.isArray(b.integrals) ? b.integrals : null,
      gapHz: Number(b.gapHz) || 15,
    });
    analysis.literature = literature(analysis, solvent);
    analysis.interpretation = interpret(analysis);
    analysis.frequencySource = frequencyFromFile ? 'calculada de la columna Hz del archivo' : 'indicada por el usuario';
    if (formula) analysis.unsaturation = unsaturation(formula);
    // DE máxima de los errores de desplazamiento (ppm) y concordancia mínima (%) para mostrar candidatos
    const sdCap = nucleus === '1H' ? 1 : 20;
    const sdMax = Number(b.sdMax) > 0 ? Math.min(sdCap, Number(b.sdMax)) : SD_MAX[nucleus];
    const minScore = Number.isFinite(Number(b.minScore)) && b.minScore !== '' && b.minScore !== null ? Math.max(0, Math.min(100, Number(b.minScore))) : MIN_SCORE;
    const cands = await findCandidates(analysis, { formula, rangeOnly: Boolean(b.rangeOnly), sdMax, minScore });
    const record = {
      name: String(b.name || '').slice(0, 200) || null,
      nucleus,
      frequency,
      solvent,
      formula,
      raw: String(b.csv).slice(0, 2_000_000),
      peaks: peaks.map((p) => ({ ppm: p.ppm, intensity: p.intensity })),
      analysis,
      candidates: cands,
    };
    let id = Number(b.uploadId) || null;
    if (id && store.getUpload(id)) store.updateUpload(id, record);
    else id = store.saveUpload(record);
    store.logSearch('espectro', record.name || `${nucleus} ${peaks.length} picos`, cands.candidates.length);
    const saved = store.getUpload(id);
    res.json({ ...saved, raw: undefined, image: `/api/uploads/${id}/image.svg` });
  }),
);

api.get('/uploads', (req, res) => res.json({ items: store.listUploads() }));

api.get('/uploads/:id', (req, res) => {
  const u = store.getUpload(Number(req.params.id));
  if (!u) return res.status(404).json({ error: 'No encontrado' });
  res.json({ ...u, image: `/api/uploads/${u.id}/image.svg` });
});

api.delete('/uploads/:id', (req, res) => res.json({ deleted: store.deleteUpload(Number(req.params.id)) }));

api.get('/uploads/:id/image.:ext', (req, res) => {
  const u = store.getUpload(Number(req.params.id));
  if (!u) return res.status(404).json({ error: 'No encontrado' });
  const svg = renderUpload(u, {
    from: num(req.query.from),
    to: num(req.query.to),
    width: Math.min(3000, num(req.query.width) || 1100),
    height: Math.min(2000, num(req.query.height) || 500),
  });
  const base = `${safeName(u.name || 'espectro')}_${u.nucleus}`;
  if (req.params.ext === 'png') {
    res.type('png');
    if (req.query.download) res.attachment(`${base}.png`);
    return res.send(svgToPng(svg));
  }
  res.type('image/svg+xml');
  if (req.query.download) res.attachment(`${base}.svg`);
  res.send(svg);
});

api.get('/uploads/:id/export', (req, res) => {
  const u = store.getUpload(Number(req.params.id));
  if (!u) return res.status(404).json({ error: 'No encontrado' });
  const base = `${safeName(u.name || 'espectro')}_${u.nucleus}`;
  if (req.query.format === 'csv') {
    const rows = [
      ['# Espectro', u.name || ''],
      ['# Nucleo', u.nucleus],
      ['# Frecuencia (MHz)', u.frequency],
      ['# Disolvente', u.solvent],
      ['# Formula', u.formula || ''],
      [],
      ['MULTIPLETES'],
      ['delta_ppm', 'desde', 'hasta', 'multiplicidad', 'J_Hz', 'nH', 'fraccion_integral', 'n_lineas', 'region'],
      ...(u.analysis?.signals || []).map((s) => [s.delta, s.from ?? '', s.to ?? '', s.multiplicity ?? '', (s.J || []).join(' / '), s.nH ?? '', s.fraction ?? '', s.nLines ?? '', s.region]),
      [],
      ['CANDIDATOS'],
      ['concordancia_%', 'DE_error_ppm', 'error_medio_ppm', 'nombre', 'formula', 'SMILES', 'fuente', 'desplazamientos_%', 'cobertura_%', 'integracion_%'],
      ...(u.candidates?.candidates || []).map((c) => [
        c.score,
        c.errorSd ?? '',
        c.errorMean ?? '',
        c.name || '',
        c.formula || '',
        c.smiles,
        c.source,
        c.shift,
        c.coverage,
        c.integration ?? '',
      ]),
    ];
    const csv = rows.map((r) => r.map((v) => (/[",\n]/.test(String(v)) ? `"${String(v).replace(/"/g, '""')}"` : v)).join(',')).join('\n');
    return res.attachment(`${base}_analisis.csv`).type('text/csv').send('﻿' + csv + '\n');
  }
  res.attachment(`${base}_analisis.json`).json(u);
});

// Espectros guardados
api.get('/compounds', (req, res) => {
  res.json(store.listCompounds({ q: String(req.query.q || ''), limit: Math.min(500, Number(req.query.limit) || 100), offset: Number(req.query.offset) || 0 }));
});

api.get('/compounds/:id', (req, res) => {
  const c = store.getCompound(Number(req.params.id));
  if (!c) return res.status(404).json({ error: 'No encontrado' });
  res.json(savedAnalysis(c));
});

api.get('/compounds/:id/export', (req, res) => {
  const c = store.getCompound(Number(req.params.id));
  if (!c) return res.status(404).json({ error: 'No encontrado' });
  const full = { compound: c, spectra: store.getCompoundSpectra(c.id) };
  res.attachment(`${safeName(c.name || c.formula)}_RMN.json`).json(full);
});

api.delete('/compounds/:id', (req, res) => {
  res.json({ deleted: store.deleteCompound(Number(req.params.id)) });
});

api.get('/spectra/:id', (req, res) => {
  const s = store.getSpectrum(Number(req.params.id));
  if (!s) return res.status(404).json({ error: 'No encontrado' });
  res.json(publicSpectrum(s));
});

api.delete('/spectra/:id', (req, res) => {
  res.json({ deleted: store.deleteSpectrum(Number(req.params.id)) });
});

function loadSpectrum(req, res) {
  const s = store.getSpectrum(Number(req.params.id));
  if (!s) {
    res.status(404).json({ error: 'Espectro no encontrado' });
    return null;
  }
  return { s, c: store.getCompound(s.compoundId) };
}

const num = (v) => (v === undefined || v === '' ? undefined : Number(v));

api.get('/spectra/:id/image.:ext', (req, res) => {
  const r = loadSpectrum(req, res);
  if (!r) return;
  const svg = renderSpectrum(r.s, r.c, {
    from: num(req.query.from),
    to: num(req.query.to),
    width: Math.min(3000, num(req.query.width) || 1100),
    height: Math.min(2000, num(req.query.height) || 500),
  });
  const base = `${safeName(r.c?.name || r.c?.formula)}_${r.s.nucleus}_${r.s.kind}`;
  if (req.params.ext === 'png') {
    res.type('png');
    if (req.query.download) res.attachment(`${base}.png`);
    return res.send(svgToPng(svg));
  }
  res.type('image/svg+xml');
  if (req.query.download) res.attachment(`${base}.svg`);
  res.send(svg);
});

api.get('/spectra/:id/export', (req, res) => {
  const r = loadSpectrum(req, res);
  if (!r) return;
  const fmt = String(req.query.format || 'json');
  const base = `${safeName(r.c?.name || r.c?.formula)}_${r.s.nucleus}_${r.s.kind}`;
  if (fmt === 'jdx') return res.attachment(`${base}.jdx`).type('text/plain').send(toJcamp(r.s, r.c));
  if (fmt === 'jdx-peaks') return res.attachment(`${base}_picos.jdx`).type('text/plain').send(toJcampPeaks(r.s, r.c));
  if (fmt === 'csv') return res.attachment(`${base}.csv`).type('text/csv').send('﻿' + toCsv(r.s, r.c));
  if (fmt === 'xy') return res.attachment(`${base}_xy.csv`).type('text/csv').send(toXyCsv(r.s));
  res.attachment(`${base}.json`).json({ compound: r.c, spectrum: r.s });
});
