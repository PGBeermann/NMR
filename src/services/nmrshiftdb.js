// Consultas a la copia local indexada de NMRShiftDB2.
import fs from 'node:fs';
import path from 'node:path';
import { Worker } from 'node:worker_threads';
import { DatabaseSync } from 'node:sqlite';
import { config } from '../config.js';

const dbFile = path.join(config.dataDir, 'nmrshiftdb.db');
let db = null;

export const status = { state: 'no disponible', detail: null, records: 0, builtAt: null, progress: null };

function open() {
  if (db) return db;
  if (!fs.existsSync(dbFile)) return null;
  db = new DatabaseSync(dbFile, { readOnly: true });
  const meta = Object.fromEntries(db.prepare('SELECT key, value FROM meta').all().map((r) => [r.key, r.value]));
  status.state = 'lista';
  status.records = Number(meta.records) || 0;
  status.builtAt = meta.built_at || null;
  status.progress = null;
  return db;
}

export function isReady() {
  return Boolean(open());
}

export function startIndexing({ force = false } = {}) {
  if (status.state === 'construyendo') return status;
  if (!force && open()) return status;
  status.state = 'construyendo';
  status.detail = 'Descargando e indexando NMRShiftDB2 (puede tardar varios minutos)…';
  const worker = new Worker(new URL('./nmrshiftdb-worker.js', import.meta.url), {
    workerData: { url: config.nmrshiftdbUrl, dataDir: config.dataDir },
  });
  worker.on('message', (m) => {
    if (m.type === 'progress') status.progress = m;
    if (m.type === 'done') {
      if (db) db.close();
      db = null;
      open();
      status.detail = `Índice construido: ${m.records} compuestos`;
    }
    if (m.type === 'error') {
      status.state = 'error';
      status.detail = m.message;
    }
  });
  worker.on('error', (e) => {
    status.state = 'error';
    status.detail = e.message;
  });
  return status;
}

function rowToEntry(r) {
  return {
    id: r.id,
    inchikey: r.inchikey,
    formula: r.formula,
    smiles: r.smiles,
    name: r.name,
    molfile: r.molfile,
    spectra: JSON.parse(r.spectra || '[]'),
    url: `https://nmrshiftdb.nmr.uni-koeln.de/portal/js_pane/P-Results/nmrshiftdbaction/showDetailsFromHome/molNumber/${r.id}`,
  };
}

/** Busca por estructura: idcode exacto, luego sin estereoquímica, luego InChIKey. */
export function findByStructure({ idcode, idcodeNoStereo, inchikey }) {
  const d = open();
  if (!d) return [];
  let rows = idcode ? d.prepare('SELECT * FROM entries WHERE idcode = ?').all(idcode) : [];
  if (!rows.length && inchikey) rows = d.prepare('SELECT * FROM entries WHERE inchikey = ?').all(inchikey);
  if (!rows.length && idcodeNoStereo) rows = d.prepare('SELECT * FROM entries WHERE idcode_ns = ?').all(idcodeNoStereo);
  if (!rows.length && inchikey) rows = d.prepare('SELECT * FROM entries WHERE ikey14 = ?').all(inchikey.slice(0, 14));
  return rows.slice(0, 10).map(rowToEntry);
}

export function findByFormula(formula) {
  const d = open();
  if (!d) return [];
  return d
    .prepare('SELECT id, inchikey, formula, smiles, name, spectra FROM entries WHERE formula = ? LIMIT 1000')
    .all(formula)
    .map((r) => {
      const spectra = JSON.parse(r.spectra || '[]');
      return {
        id: r.id,
        inchikey: r.inchikey,
        formula: r.formula,
        smiles: r.smiles,
        name: r.name,
        nuclei: [...new Set(spectra.map((s) => s.nucleus))],
      };
    });
}

// Índice en memoria de todas las señales por núcleo (para buscar por espectro).
const spectraCache = {};

function groupShifts(signals, nucleus) {
  // 1H: cada entrada es un H (con el átomo pesado que lo porta) → nº de H por señal
  const tol = nucleus === '1H' ? 0.015 : 0.05;
  const sorted = [...signals].sort((a, b) => a.delta - b.delta);
  const out = [];
  for (const s of sorted) {
    const g = out[out.length - 1];
    if (g && s.delta - g.last <= tol) {
      g.sum += s.delta;
      g.n++;
      g.last = s.delta;
      if (s.atom !== null && s.atom !== undefined) g.atoms.add(s.atom);
    } else out.push({ sum: s.delta, n: 1, last: s.delta, atoms: new Set(s.atom !== null && s.atom !== undefined ? [s.atom] : []) });
  }
  return out.map((g) => ({ delta: g.sum / g.n, n: g.n, atoms: [...g.atoms] }));
}

export function getAllSpectra(nucleus) {
  if (spectraCache[nucleus]) return spectraCache[nucleus];
  const d = open();
  if (!d) return null;
  const like = nucleus === '1H' ? '%"nucleus":"1H"%' : '%"nucleus":"13C"%';
  const list = [];
  for (const r of d.prepare('SELECT id, inchikey, formula, smiles, spectra FROM entries WHERE spectra LIKE ?').iterate(like)) {
    const specs = JSON.parse(r.spectra)
      .filter((s) => s.nucleus === nucleus && s.signals.length)
      .map((s) => ({ calculated: s.calculated, solvent: s.solvent, field: s.field, signals: groupShifts(s.signals, nucleus) }));
    if (specs.length) list.push({ id: r.id, inchikey: r.inchikey, formula: r.formula, smiles: r.smiles, spectra: specs });
  }
  spectraCache[nucleus] = list;
  return list;
}

export function getMolfile(id) {
  const d = open();
  return d ? d.prepare('SELECT molfile FROM entries WHERE id = ?').get(id)?.molfile || null : null;
}

export function getEntry(id) {
  const d = open();
  if (!d) return null;
  const r = d.prepare('SELECT * FROM entries WHERE id = ?').get(id);
  return r ? rowToEntry(r) : null;
}
