// Base de datos local (SQLite integrado en Node) donde se guardan todos los
// compuestos analizados y sus espectros (señales, intensidades, constantes J…).
import fs from 'node:fs';
import path from 'node:path';
import { DatabaseSync } from 'node:sqlite';
import { config } from './config.js';

fs.mkdirSync(config.dataDir, { recursive: true });

export const db = new DatabaseSync(path.join(config.dataDir, 'nmr.db'));
db.exec('PRAGMA journal_mode = WAL; PRAGMA foreign_keys = ON;');

db.exec(`
CREATE TABLE IF NOT EXISTS compounds (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  idcode TEXT UNIQUE NOT NULL,
  name TEXT,
  iupac TEXT,
  formula TEXT,
  mw REAL,
  exact_mass REAL,
  smiles TEXT,
  inchikey TEXT,
  cid INTEGER,
  molfile TEXT,
  info TEXT,
  created_at TEXT NOT NULL,
  updated_at TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS spectra (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  compound_id INTEGER NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
  nucleus TEXT NOT NULL,
  kind TEXT NOT NULL,
  source TEXT NOT NULL,
  source_id TEXT NOT NULL DEFAULT '',
  frequency REAL NOT NULL DEFAULT 0,
  solvent TEXT,
  signals TEXT,
  peaks TEXT,
  analysis TEXT,
  meta TEXT,
  created_at TEXT NOT NULL,
  updated_at TEXT NOT NULL,
  UNIQUE (compound_id, nucleus, kind, source, source_id, frequency)
);
CREATE INDEX IF NOT EXISTS spectra_compound ON spectra(compound_id);
CREATE TABLE IF NOT EXISTS prediction_cache (
  idcode TEXT PRIMARY KEY,
  proton TEXT,
  carbon TEXT,
  created_at TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS uploads (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  name TEXT,
  nucleus TEXT NOT NULL,
  frequency REAL,
  solvent TEXT,
  formula TEXT,
  raw TEXT,
  peaks TEXT,
  analysis TEXT,
  candidates TEXT,
  created_at TEXT NOT NULL,
  updated_at TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS searches (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  mode TEXT NOT NULL,
  query TEXT NOT NULL,
  results INTEGER,
  created_at TEXT NOT NULL
);
`);

// Migración: curva del espectro cargado (JCAMP-DX), reducida, para dibujarlo
if (!db.prepare('PRAGMA table_info(uploads)').all().some((c) => c.name === 'spectrum')) db.exec('ALTER TABLE uploads ADD COLUMN spectrum TEXT');

const now = () => new Date().toISOString();
const j = (v) => (v === undefined ? null : JSON.stringify(v));
const p = (v) => (v ? JSON.parse(v) : null);

export function logSearch(mode, query, results) {
  db.prepare('INSERT INTO searches (mode, query, results, created_at) VALUES (?, ?, ?, ?)').run(
    mode,
    String(query).slice(0, 2000),
    results ?? null,
    now(),
  );
}

export function upsertCompound(c) {
  const t = now();
  db.prepare(
    `INSERT INTO compounds (idcode, name, iupac, formula, mw, exact_mass, smiles, inchikey, cid, molfile, info, created_at, updated_at)
     VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
     ON CONFLICT(idcode) DO UPDATE SET
       name = COALESCE(excluded.name, compounds.name),
       iupac = COALESCE(excluded.iupac, compounds.iupac),
       formula = excluded.formula, mw = excluded.mw, exact_mass = excluded.exact_mass,
       smiles = excluded.smiles,
       inchikey = COALESCE(excluded.inchikey, compounds.inchikey),
       cid = COALESCE(excluded.cid, compounds.cid),
       molfile = excluded.molfile, info = excluded.info, updated_at = excluded.updated_at`,
  ).run(
    c.idcode,
    c.name ?? null,
    c.iupac ?? null,
    c.formula ?? null,
    c.mw ?? null,
    c.exactMass ?? null,
    c.smiles ?? null,
    c.inchikey ?? null,
    c.cid ?? null,
    c.molfile ?? null,
    j(c.info),
    t,
    t,
  );
  return db.prepare('SELECT id FROM compounds WHERE idcode = ?').get(c.idcode).id;
}

export function upsertSpectrum(compoundId, s) {
  const t = now();
  db.prepare(
    `INSERT INTO spectra (compound_id, nucleus, kind, source, source_id, frequency, solvent, signals, peaks, analysis, meta, created_at, updated_at)
     VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
     ON CONFLICT(compound_id, nucleus, kind, source, source_id, frequency) DO UPDATE SET
       solvent = excluded.solvent, signals = excluded.signals, peaks = excluded.peaks,
       analysis = excluded.analysis, meta = excluded.meta, updated_at = excluded.updated_at`,
  ).run(
    compoundId,
    s.nucleus,
    s.kind,
    s.source,
    String(s.sourceId ?? ''),
    Number(s.frequency) || 0,
    s.solvent ?? null,
    j(s.signals),
    j(s.peaks),
    j(s.analysis),
    j(s.meta),
    t,
    t,
  );
  return db
    .prepare(
      'SELECT id FROM spectra WHERE compound_id = ? AND nucleus = ? AND kind = ? AND source = ? AND source_id = ? AND frequency = ?',
    )
    .get(compoundId, s.nucleus, s.kind, s.source, String(s.sourceId ?? ''), Number(s.frequency) || 0).id;
}

function rowToSpectrum(r) {
  if (!r) return null;
  return {
    id: r.id,
    compoundId: r.compound_id,
    nucleus: r.nucleus,
    kind: r.kind,
    source: r.source,
    sourceId: r.source_id,
    frequency: r.frequency,
    solvent: r.solvent,
    signals: p(r.signals),
    peaks: p(r.peaks),
    analysis: p(r.analysis),
    meta: p(r.meta),
    createdAt: r.created_at,
    updatedAt: r.updated_at,
  };
}

function rowToCompound(r) {
  if (!r) return null;
  return {
    id: r.id,
    idcode: r.idcode,
    name: r.name,
    iupac: r.iupac,
    formula: r.formula,
    mw: r.mw,
    exactMass: r.exact_mass,
    smiles: r.smiles,
    inchikey: r.inchikey,
    cid: r.cid,
    molfile: r.molfile,
    info: p(r.info),
    createdAt: r.created_at,
    updatedAt: r.updated_at,
  };
}

export function getSpectrum(id) {
  return rowToSpectrum(db.prepare('SELECT * FROM spectra WHERE id = ?').get(id));
}

export function getCompound(id) {
  return rowToCompound(db.prepare('SELECT * FROM compounds WHERE id = ?').get(id));
}

export function getCompoundSpectra(compoundId) {
  return db
    .prepare('SELECT * FROM spectra WHERE compound_id = ? ORDER BY kind DESC, nucleus, id')
    .all(compoundId)
    .map(rowToSpectrum);
}

export function listCompounds({ q = '', limit = 100, offset = 0 } = {}) {
  const like = `%${q}%`;
  const rows = db
    .prepare(
      `SELECT c.id, c.name, c.iupac, c.formula, c.smiles, c.cid, c.inchikey, c.updated_at,
              (SELECT COUNT(*) FROM spectra s WHERE s.compound_id = c.id) AS n_spectra,
              (SELECT COUNT(*) FROM spectra s WHERE s.compound_id = c.id AND s.kind = 'experimental') AS n_exp
       FROM compounds c
       WHERE (? = '' OR c.name LIKE ? OR c.iupac LIKE ? OR c.formula LIKE ? OR c.smiles LIKE ?)
       ORDER BY c.updated_at DESC LIMIT ? OFFSET ?`,
    )
    .all(q, like, like, like, like, limit, offset);
  const total = db
    .prepare(
      `SELECT COUNT(*) AS n FROM compounds c
       WHERE (? = '' OR c.name LIKE ? OR c.iupac LIKE ? OR c.formula LIKE ? OR c.smiles LIKE ?)`,
    )
    .get(q, like, like, like, like).n;
  return { total, items: rows };
}

export function deleteCompound(id) {
  return db.prepare('DELETE FROM compounds WHERE id = ?').run(id).changes;
}

export function deleteSpectrum(id) {
  return db.prepare('DELETE FROM spectra WHERE id = ?').run(id).changes;
}

// Caché persistente de predicciones simplificadas (señales agrupadas) para candidatos
export function getCachedPrediction(idcode) {
  const r = db.prepare('SELECT proton, carbon FROM prediction_cache WHERE idcode = ?').get(idcode);
  return r ? { proton: p(r.proton), carbon: p(r.carbon) } : null;
}
export function setCachedPrediction(idcode, proton, carbon) {
  db.prepare('INSERT OR REPLACE INTO prediction_cache (idcode, proton, carbon, created_at) VALUES (?, ?, ?, ?)').run(idcode, j(proton), j(carbon), now());
}

// Espectros cargados por el usuario (CSV)
export function saveUpload(u) {
  const t = now();
  const info = db
    .prepare(
      `INSERT INTO uploads (name, nucleus, frequency, solvent, formula, raw, peaks, analysis, candidates, spectrum, created_at, updated_at)
       VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)`,
    )
    .run(u.name ?? null, u.nucleus, u.frequency ?? null, u.solvent ?? null, u.formula ?? null, u.raw ?? null, j(u.peaks), j(u.analysis), j(u.candidates), j(u.spectrum ?? null), t, t);
  return Number(info.lastInsertRowid);
}
export function updateUpload(id, fields) {
  const t = now();
  db.prepare('UPDATE uploads SET analysis = ?, candidates = ?, formula = ?, updated_at = ? WHERE id = ?').run(
    j(fields.analysis),
    j(fields.candidates),
    fields.formula ?? null,
    t,
    id,
  );
}
export function getUpload(id) {
  const r = db.prepare('SELECT * FROM uploads WHERE id = ?').get(id);
  if (!r) return null;
  return {
    id: r.id,
    name: r.name,
    nucleus: r.nucleus,
    frequency: r.frequency,
    solvent: r.solvent,
    formula: r.formula,
    raw: r.raw,
    peaks: p(r.peaks),
    analysis: p(r.analysis),
    candidates: p(r.candidates),
    spectrum: p(r.spectrum),
    createdAt: r.created_at,
    updatedAt: r.updated_at,
  };
}
export function listUploads() {
  return db
    .prepare('SELECT id, name, nucleus, frequency, solvent, formula, created_at, updated_at FROM uploads ORDER BY updated_at DESC LIMIT 200')
    .all();
}
export function deleteUpload(id) {
  return db.prepare('DELETE FROM uploads WHERE id = ?').run(id).changes;
}

export function stats() {
  return {
    compounds: db.prepare('SELECT COUNT(*) AS n FROM compounds').get().n,
    spectra: db.prepare('SELECT COUNT(*) AS n FROM spectra').get().n,
    searches: db.prepare('SELECT COUNT(*) AS n FROM searches').get().n,
  };
}
