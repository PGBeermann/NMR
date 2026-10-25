import { Router } from 'express';
import * as store from '../db.js';
import * as pubchem from '../services/pubchem.js';
import * as opsin from '../services/opsin.js';
import * as nmrshiftdb from '../services/nmrshiftdb.js';
import { analyze, publicSpectrum, savedAnalysis } from '../services/analysis.js';
import { normalizeFormula } from '../services/molecule.js';
import { SOLVENTS } from '../services/chemistry.js';
import { renderSpectrum } from '../services/render.js';
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
