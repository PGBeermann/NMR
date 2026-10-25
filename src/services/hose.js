// Códigos de entorno tipo HOSE (esferas 0–4) rápidos, calculados con OpenChemLib,
// y predictor local de desplazamientos basado en NMRShiftDB2 (sin servicios externos).
import fs from 'node:fs';
import path from 'node:path';
import { DatabaseSync } from 'node:sqlite';
import OCL from 'openchemlib';
import { config } from '../config.js';

export const MAX_SPHERE = 4;

function atomCode(mol, a) {
  const z = mol.getAtomicNo(a);
  const ar = mol.isAromaticAtom(a) ? 'a' : '';
  const ch = mol.getAtomCharge(a);
  return `${z}${ar}${mol.getAllHydrogens(a)}${ch ? `c${ch}` : ''}`;
}

/** FNV-1a de 32 bits → entero (claves compactas para la tabla). */
function hash(str) {
  let h = 0x811c9dc5;
  for (let i = 0; i < str.length; i++) {
    h ^= str.charCodeAt(i);
    h = Math.imul(h, 0x01000193);
  }
  return h >>> 0;
}

/**
 * Códigos (hash) por esfera para el átomo `a`. Esfera k = átomos a k enlaces,
 * con el orden de enlace hacia la capa anterior ("a" = aromático), ordenados.
 */
export function hoseKeys(mol, a, maxSphere = MAX_SPHERE) {
  const keys = [];
  let code = atomCode(mol, a);
  keys.push(hash(`0:${code}`));
  const seen = new Set([a]);
  let layer = [a];
  for (let k = 1; k <= maxSphere; k++) {
    const parts = [];
    const next = [];
    for (const u of layer) {
      for (let i = 0; i < mol.getConnAtoms(u); i++) {
        const v = mol.getConnAtom(u, i);
        if (seen.has(v)) continue;
        seen.add(v);
        const b = mol.getConnBond(u, i);
        const bo = mol.isAromaticBond(b) ? 'a' : mol.getBondOrder(b);
        parts.push(`${bo}${atomCode(mol, v)}`);
        next.push(v);
      }
    }
    if (!parts.length) break;
    parts.sort();
    code += `|${parts.join(',')}`;
    keys.push(hash(`${k}:${code}`));
    layer = next;
  }
  return keys;
}

export function prepareMolecule(mol) {
  const m = mol.getCompactCopy();
  m.removeExplicitHydrogens();
  m.ensureHelperArrays(OCL.Molecule.cHelperRings);
  return m;
}

// ---------- Construcción de la tabla (desde nmrshiftdb.db) ----------

export function buildHoseTable({ dataDir = config.dataDir, onProgress } = {}) {
  const src = new DatabaseSync(path.join(dataDir, 'nmrshiftdb.db'), { readOnly: true });
  const acc = { '13C': new Map(), '1H': new Map() };
  const add = (nuc, key, sphere, value) => {
    const k = `${key}`;
    let s = acc[nuc].get(k);
    if (!s) acc[nuc].set(k, (s = [sphere, 0, 0, 0]));
    s[1]++;
    s[2] += value;
    s[3] += value * value;
  };
  let count = 0;
  for (const r of src.prepare('SELECT molfile, spectra FROM entries').iterate()) {
    count++;
    if (count % 5000 === 0) onProgress?.({ phase: 'hose', records: count });
    let mol;
    try {
      mol = prepareMolecule(OCL.Molecule.fromMolfile(r.molfile));
    } catch {
      continue;
    }
    const n = mol.getAtoms();
    for (const sp of JSON.parse(r.spectra)) {
      if (sp.calculated || (sp.nucleus !== '13C' && sp.nucleus !== '1H')) continue;
      // Valor por átomo pesado (media de los H del mismo átomo en 1H)
      const byAtom = new Map();
      for (const s of sp.signals) {
        if (s.atom === null || s.atom === undefined || s.atom >= n) continue;
        if (!byAtom.has(s.atom)) byAtom.set(s.atom, []);
        byAtom.get(s.atom).push(s.delta);
      }
      for (const [atom, vals] of byAtom) {
        if (sp.nucleus === '13C' && mol.getAtomicNo(atom) !== 6) continue;
        if (sp.nucleus === '1H' && mol.getAllHydrogens(atom) === 0) continue;
        const v = vals.reduce((x, y) => x + y, 0) / vals.length;
        const keys = hoseKeys(mol, atom);
        keys.forEach((key, sphere) => add(sp.nucleus, key, sphere, v));
      }
    }
  }
  src.close();
  const tmp = path.join(dataDir, 'hose.db.tmp');
  fs.rmSync(tmp, { force: true });
  const db = new DatabaseSync(tmp);
  db.exec(`PRAGMA journal_mode = OFF; PRAGMA synchronous = OFF;
    CREATE TABLE hose (nucleus TEXT, key INTEGER, sphere INTEGER, n INTEGER, mean REAL, sd REAL, PRIMARY KEY (nucleus, key)) WITHOUT ROWID;
    CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT);`);
  const ins = db.prepare('INSERT OR REPLACE INTO hose VALUES (?, ?, ?, ?, ?, ?)');
  db.exec('BEGIN');
  let rows = 0;
  for (const nuc of ['13C', '1H']) {
    for (const [k, [sphere, nn, sum, sumsq]] of acc[nuc]) {
      const mean = sum / nn;
      const sd = nn > 1 ? Math.sqrt(Math.max(0, sumsq / nn - mean * mean)) : null;
      ins.run(nuc, Number(k), sphere, nn, mean, sd);
      if (++rows % 200000 === 0) db.exec('COMMIT; BEGIN');
    }
  }
  db.exec('COMMIT');
  db.prepare('INSERT INTO meta VALUES (?, ?)').run('built_at', new Date().toISOString());
  db.prepare('INSERT INTO meta VALUES (?, ?)').run('rows', String(rows));
  db.close();
  fs.renameSync(tmp, path.join(dataDir, 'hose.db'));
  onProgress?.({ phase: 'hose-listo', rows });
  return { rows, molecules: count };
}

// ---------- Predicción local ----------

let hdb = null;
let lookup = null;
function openHose() {
  if (hdb) return true;
  const f = path.join(config.dataDir, 'hose.db');
  if (!fs.existsSync(f)) return false;
  hdb = new DatabaseSync(f, { readOnly: true });
  lookup = hdb.prepare('SELECT sphere, n, mean, sd FROM hose WHERE nucleus = ? AND key = ?');
  return true;
}

export function hoseReady() {
  return openHose();
}

/** Desplazamiento predicho para un átomo: la esfera más específica con datos. */
function predictAtom(mol, a, nucleus) {
  const keys = hoseKeys(mol, a);
  for (let s = keys.length - 1; s >= 1; s--) {
    const row = lookup.get(nucleus, keys[s]);
    if (row && (row.n >= 2 || s >= 3)) return { delta: row.mean, sphere: s, n: row.n, sd: row.sd };
  }
  return null;
}

/** H vecinos (a 3 enlaces) de los H unidos a `a`, sin OH/NH/SH. */
export function vicinalH(mol, a) {
  if ([7, 8, 16].includes(mol.getAtomicNo(a))) return 0;
  let n = 0;
  for (let i = 0; i < mol.getConnAtoms(a); i++) {
    const b = mol.getConnAtom(a, i);
    if (![7, 8, 16].includes(mol.getAtomicNo(b))) n += mol.getAllHydrogens(b);
  }
  return n;
}

/**
 * Predicción local de 1H o 13C. Agrupa átomos equivalentes por simetría
 * (IDs diastereotópicos de OCL) y devuelve [{delta, n, nb, atoms, sphere}].
 */
export function predictLocal(molInput, nucleus) {
  if (!openHose()) return null;
  const mol = prepareMolecule(molInput);
  const ids = mol.getDiastereotopicAtomIDs();
  const groups = new Map();
  for (let a = 0; a < mol.getAtoms(); a++) {
    const z = mol.getAtomicNo(a);
    const nH = mol.getAllHydrogens(a);
    if (nucleus === '13C' && z !== 6) continue;
    if (nucleus === '1H' && (nH === 0 || [7, 8, 16].includes(z))) continue; // OH/NH se tratan aparte
    const p = predictAtom(mol, a, nucleus);
    if (!p) continue;
    const key = ids[a];
    if (!groups.has(key)) groups.set(key, { deltaSum: 0, k: 0, n: 0, atoms: [], nb: vicinalH(mol, a), sphere: p.sphere });
    const g = groups.get(key);
    g.deltaSum += p.delta;
    g.k++;
    g.n += nucleus === '1H' ? nH : 1;
    g.atoms.push(a);
    g.sphere = Math.min(g.sphere, p.sphere);
  }
  return [...groups.values()].map((g) => ({ delta: g.deltaSum / g.k, n: g.n, nb: nucleus === '1H' ? g.nb : null, atoms: g.atoms, sphere: g.sphere }));
}
