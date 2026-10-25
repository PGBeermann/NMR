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

export function getEntry(id) {
  const d = open();
  if (!d) return null;
  const r = d.prepare('SELECT * FROM entries WHERE id = ?').get(id);
  return r ? rowToEntry(r) : null;
}
