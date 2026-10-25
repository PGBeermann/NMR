// Descarga el archivo SD de NMRShiftDB2 y lo convierte en una base SQLite
// indexada por estructura (idcode de OpenChemLib), InChIKey y fórmula.
import fs from 'node:fs';
import path from 'node:path';
import readline from 'node:readline';
import { DatabaseSync } from 'node:sqlite';
import OCL from 'openchemlib';

export async function download(url, dest, onProgress) {
  const res = await fetch(url, { redirect: 'follow', headers: { 'User-Agent': 'nmr-analyzer/1.0' } });
  if (!res.ok) throw new Error(`Descarga fallida (${res.status})`);
  const total = Number(res.headers.get('content-length')) || 0;
  let done = 0;
  let last = 0;
  const tmp = `${dest}.part`;
  const fd = fs.openSync(tmp, 'w');
  try {
    const reader = res.body.getReader();
    for (;;) {
      const { done: end, value } = await reader.read();
      if (end) break;
      fs.writeSync(fd, value);
      done += value.length;
      if (done - last > 2_000_000) {
        last = done;
        onProgress?.({ phase: 'descargando', bytes: done, total });
      }
    }
  } finally {
    fs.closeSync(fd);
  }
  if (total && done !== total) throw new Error(`Descarga incompleta (${done} de ${total} bytes)`);
  fs.renameSync(tmp, dest);
}

// "0:298 " o "2:Quantum ESPRESSO 3:Quantum ESPRESSO " → {0:'298', 2:'Quantum ESPRESSO', 3:'Quantum ESPRESSO'}
function parseIndexed(value) {
  const out = {};
  const parts = (' ' + value).split(/\s(\d+):/);
  for (let i = 1; i < parts.length; i += 2) out[parts[i]] = (parts[i + 1] || '').trim();
  return out;
}

const CALC_FIELDS = ['Program', 'NMRProgram', 'NMRMethod', 'NMRBasisSet', 'GeomMethod'];

function parseRecord(lines) {
  const endIdx = lines.findIndex((l) => l.startsWith('M  END'));
  if (endIdx < 0) return null;
  const molfile = lines.slice(0, endIdx + 1).join('\n');
  const props = {};
  for (let i = endIdx + 1; i < lines.length; i++) {
    const m = lines[i].match(/^>\s*<(.+?)>/);
    if (m) {
      const vals = [];
      while (i + 1 < lines.length && lines[i + 1].trim() !== '' && !lines[i + 1].startsWith('>')) vals.push(lines[++i]);
      props[m[1]] = vals.join(' ');
    }
  }
  const calculated = new Set();
  for (const f of CALC_FIELDS) if (props[f]) for (const k of Object.keys(parseIndexed(props[f]))) calculated.add(k);
  const solvents = parseIndexed(props.Solvent || '');
  const fields = parseIndexed(props['Field Strength [MHz]'] || '');
  const temps = parseIndexed(props['Temperature [K]'] || '');
  const programs = parseIndexed(props.Program || props.NMRProgram || '');
  const assign = parseIndexed(props['Assignment Method'] || '');

  const spectra = [];
  for (const [key, value] of Object.entries(props)) {
    const m = key.match(/^Spectrum (\S+) (\d+)$/);
    if (!m) continue;
    const [, nucleus, num] = m;
    const signals = [];
    for (const part of value.split('|')) {
      const [shift, multField, atom] = part.split(';');
      if (shift === undefined || shift.trim() === '' || Number.isNaN(Number(shift))) continue;
      const mult = (multField || '').replace(/[\d.]/g, '').trim();
      signals.push({ delta: Number(shift), mult: mult || null, atom: atom !== undefined && atom !== '' ? Number(atom) : null });
    }
    spectra.push({
      nucleus,
      number: Number(num),
      calculated: calculated.has(num),
      program: programs[num] || null,
      solvent: solvents[num] || null,
      field: Number(fields[num]) || null,
      temperature: temps[num] || null,
      assignment: assign[num] || null,
      signals,
    });
  }
  return {
    id: Number(props['nmrshiftdb2 ID']),
    inchikey: props['INChI key']?.trim() || null,
    inchi: props.INChI?.trim() || null,
    name: (props.CHEMNAME || props.Name || '').trim() || null,
    molfile,
    spectra,
  };
}

export async function buildIndex({ url, dataDir, onProgress } = {}) {
  fs.mkdirSync(dataDir, { recursive: true });
  const sdFile = path.join(dataDir, 'nmrshiftdb2withsignals.sd');
  if (!fs.existsSync(sdFile)) {
    onProgress?.({ phase: 'descargando', bytes: 0, total: 0 });
    await download(url, sdFile, onProgress);
  }
  const tmpDb = path.join(dataDir, 'nmrshiftdb.db.tmp');
  for (const f of [tmpDb, `${tmpDb}-wal`, `${tmpDb}-shm`]) fs.rmSync(f, { force: true });
  const db = new DatabaseSync(tmpDb);
  db.exec(`
    PRAGMA journal_mode = OFF; PRAGMA synchronous = OFF;
    CREATE TABLE entries (
      id INTEGER PRIMARY KEY, inchikey TEXT, ikey14 TEXT, formula TEXT,
      idcode TEXT, idcode_ns TEXT, smiles TEXT, name TEXT, molfile TEXT, spectra TEXT);
    CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT);
  `);
  const ins = db.prepare(
    'INSERT OR REPLACE INTO entries (id, inchikey, ikey14, formula, idcode, idcode_ns, smiles, name, molfile, spectra) VALUES (?,?,?,?,?,?,?,?,?,?)',
  );

  const rl = readline.createInterface({ input: fs.createReadStream(sdFile, { encoding: 'latin1' }), crlfDelay: Infinity });
  let lines = [];
  let count = 0;
  let errors = 0;
  db.exec('BEGIN');
  for await (const line of rl) {
    if (line.startsWith('$$$$')) {
      try {
        const rec = parseRecord(lines);
        if (rec && rec.id && rec.spectra.length) {
          const mol = OCL.Molecule.fromMolfile(rec.molfile);
          const ns = mol.getCompactCopy();
          ns.stripStereoInformation();
          const formula = mol.getMolecularFormula().formula;
          ins.run(
            rec.id,
            rec.inchikey,
            rec.inchikey ? rec.inchikey.slice(0, 14) : null,
            formula,
            mol.getIDCode(),
            ns.getIDCode(),
            mol.toIsomericSmiles(),
            rec.name,
            rec.molfile,
            JSON.stringify(rec.spectra),
          );
          count++;
          if (count % 2000 === 0) {
            db.exec('COMMIT; BEGIN');
            onProgress?.({ phase: 'indexando', records: count });
          }
        }
      } catch {
        errors++;
      }
      lines = [];
    } else {
      lines.push(line);
    }
  }
  db.exec('COMMIT');
  db.exec(`
    CREATE INDEX idx_formula ON entries(formula);
    CREATE INDEX idx_idcode ON entries(idcode);
    CREATE INDEX idx_idcode_ns ON entries(idcode_ns);
    CREATE INDEX idx_ikey ON entries(inchikey);
    CREATE INDEX idx_ikey14 ON entries(ikey14);
  `);
  const setMeta = db.prepare('INSERT OR REPLACE INTO meta (key, value) VALUES (?, ?)');
  setMeta.run('built_at', new Date().toISOString());
  setMeta.run('records', String(count));
  setMeta.run('errors', String(errors));
  setMeta.run('source', url);
  db.close();
  const finalDb = path.join(dataDir, 'nmrshiftdb.db');
  fs.renameSync(tmpDb, finalDb);
  // El archivo SD original ya no es necesario (ahorra ~160 MB en el VPS)
  if (process.env.NMRSHIFTDB_KEEP_SD !== '1') fs.rmSync(sdFile, { force: true });
  onProgress?.({ phase: 'listo', records: count, errors });
  return { records: count, errors, file: finalDb };
}
