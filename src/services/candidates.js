// Propuesta de estructuras a partir de un espectro experimental (lista de picos).
// 1) Busca en todos los espectros de NMRShiftDB2 (experimentales y calculados).
// 2) Si se indica la fórmula molecular, predice los espectros de los compuestos de
//    PubChem con esa fórmula y los compara.
// La concordancia (%) combina: cercanía de desplazamientos, cobertura de las señales
// del candidato e integración (sólo 1H).
import { getAllSpectra, getMolfile } from './nmrshiftdb.js';
import * as pubchem from './pubchem.js';
import { OCL, moleculeSvg, parseMolecule } from './molecule.js';
import { runPrediction } from './analysis.js';
import { parseMultiplicity } from './chemistry.js';
import { getCachedPrediction, setCachedPrediction } from '../db.js';

const r = (v, d = 1) => Math.round(v * 10 ** d) / 10 ** d;
const TOL = { '1H': 0.12, '13C': 2.5 };
const MARGIN = { '1H': 0.3, '13C': 5 };

/** Nº de H vecinos implicado por una multiplicidad (regla n+1): t → 2, dd → 2, q → 3… ('m' → desconocido). */
export function neighborsFromMultiplicity(m) {
  if (!m) return null;
  const parts = parseMultiplicity(m).filter((p) => p.token !== 'br');
  if (!parts.length || parts.some((p) => p.n === null)) return null;
  return parts.reduce((a, p) => a + p.n, 0);
}

/** H vecinos (a 3 enlaces) de los H unidos al átomo pesado `a`, sin contar OH/NH/SH. */
function vicinalH(mol, a) {
  if ([7, 8, 16].includes(mol.getAtomicNo(a))) return 0;
  let n = 0;
  for (let i = 0; i < mol.getConnAtoms(a); i++) {
    const b = mol.getConnAtom(a, i);
    if (![7, 8, 16].includes(mol.getAtomicNo(b))) n += mol.getAllHydrogens(b);
  }
  return n;
}

/**
 * @param exp [{delta, w}] señales experimentales (w = fracción de la integral)
 * @param cand [{delta, n}] señales del candidato (n = nº de H o de C)
 */
export function concordance(exp, cand, nucleus, { rangeOnly = false, totalH = null, expTotalH = null } = {}) {
  const tol = TOL[nucleus];
  let c = cand;
  if (rangeOnly) {
    const lo = Math.min(...exp.map((e) => e.delta)) - MARGIN[nucleus];
    const hi = Math.max(...exp.map((e) => e.delta)) + MARGIN[nucleus];
    c = cand.filter((x) => x.delta >= lo && x.delta <= hi);
  }
  if (!c.length || !exp.length) return null;
  const sumN = c.reduce((a, x) => a + x.n, 0) || 1;
  const f = c.map((x) => x.n / sumN);
  const sim = (d) => Math.exp(-((d / tol) ** 2));

  // Asignación uno a uno (voraz, por cercanía): cada señal experimental a una señal distinta del candidato
  const pairs = [];
  for (let i = 0; i < exp.length; i++) for (let j = 0; j < c.length; j++) pairs.push([i, j, Math.abs(exp[i].delta - c[j].delta)]);
  pairs.sort((a, b) => a[2] - b[2]);
  const ei = new Array(exp.length).fill(-1);
  const cj = new Array(c.length).fill(-1);
  for (const [i, j] of pairs) {
    if (ei[i] < 0 && cj[j] < 0) {
      ei[i] = j;
      cj[j] = i;
    }
  }
  const nearest = (d0, list) => list.reduce((m, x) => Math.min(m, Math.abs(x.delta - d0)), Infinity);
  let shift = 0;
  let integErr = 0;
  let multSum = 0;
  let multN = 0;
  const matches = [];
  exp.forEach((e, i) => {
    if (ei[i] >= 0) {
      const x = c[ei[i]];
      const d = Math.abs(x.delta - e.delta);
      shift += e.w * sim(d);
      integErr += Math.abs(e.w - f[ei[i]]);
      let multOk = null;
      if (e.nb !== null && e.nb !== undefined && x.nb !== null && x.nb !== undefined) {
        const dn = Math.abs(e.nb - x.nb);
        multOk = dn === 0 ? 1 : dn === 1 ? 0.5 : 0;
        multSum += e.w * multOk;
        multN += e.w;
      }
      matches.push({
        experimental: e.delta,
        candidate: r(x.delta, 2),
        diff: r(x.delta - e.delta, 2),
        n: x.n,
        expectedNeighbors: x.nb ?? null,
        multiplicityOk: multOk,
        note: x.note || null,
      });
    } else {
      // Señal sobrante: crédito parcial si cae sobre una señal ya asignada (multiplete desdoblado)
      shift += 0.4 * e.w * sim(nearest(e.delta, c));
      integErr += e.w;
      matches.push({ experimental: e.delta, candidate: null, diff: null, n: null, note: 'sin señal propia en el candidato' });
    }
  });
  let coverage = 0;
  c.forEach((x, j) => {
    if (cj[j] >= 0) coverage += f[j] * sim(Math.abs(x.delta - exp[cj[j]].delta));
    else {
      coverage += 0.4 * f[j] * sim(nearest(x.delta, exp));
      integErr += f[j];
    }
  });
  const count = Math.min(exp.length, c.length) / Math.max(exp.length, c.length);
  let total;
  let integration = null;
  let hAgreement = null;
  let multiplicity = null;
  if (nucleus === '1H') {
    integration = Math.max(0, 1 - 0.5 * integErr);
    // nº total de H: sólo se usa si las integrales del espectro son enteras fiables (dadas por el usuario)
    let agree = count;
    if (!rangeOnly && totalH && expTotalH) {
      hAgreement = Math.min(totalH, expTotalH) / Math.max(totalH, expTotalH);
      agree = (count + hAgreement) / 2;
    }
    // Acuerdo de multiplicidades (regla n+1) cuando se puede calcular
    multiplicity = multN > 0.3 ? multSum / multN : null;
    total =
      multiplicity === null
        ? 0.4 * shift + 0.25 * coverage + 0.2 * integration + 0.15 * agree
        : 0.35 * shift + 0.2 * coverage + 0.15 * integration + 0.1 * agree + 0.2 * multiplicity;
  } else {
    total = (0.5 * shift + 0.5 * coverage) * (0.7 + 0.3 * count);
  }
  return {
    score: r(total * 100, 1),
    shift: r(shift * 100, 1),
    coverage: r(coverage * 100, 1),
    integration: integration === null ? null : r(integration * 100, 1),
    signalCount: r(count * 100, 1),
    hAgreement: hAgreement === null ? null : r(hAgreement * 100, 1),
    multiplicity: multiplicity === null ? null : r(multiplicity * 100, 1),
    unexplained: c.filter((x, j) => cj[j] < 0).map((x) => r(x.delta, 2)),
    matches,
  };
}

function expSignals(analysis) {
  if (analysis.nucleus === '13C') {
    const n = analysis.signals.length;
    return analysis.signals.map((s) => ({ delta: s.delta, w: 1 / n }));
  }
  return analysis.signals.map((s) => ({ delta: s.delta, w: s.fraction, mult: s.multiplicity, nb: neighborsFromMultiplicity(s.multiplicity) }));
}

const hCount = (formula) => Number(String(formula || '').match(/H(\d*)(?![a-z])/)?.[1] || (/H(?![a-z])/.test(formula || '') ? 1 : 0));

/**
 * Los H de OH/NH/SH no suelen aparecer en las bases de datos ni en la predicción.
 * Si el candidato los tiene y el espectro muestra un singlete sin asignar, se
 * considera una posible señal intercambiable (con crédito parcial).
 */
function withExchangeable(exp, candSignals, smiles, nucleus, totalH) {
  if (nucleus !== '1H' || !smiles) return null;
  let exch = 0;
  try {
    const mol = parseMolecule({ smiles });
    for (let a = 0; a < mol.getAtoms(); a++) if ([7, 8, 16].includes(mol.getAtomicNo(a))) exch += mol.getAllHydrogens(a);
  } catch {
    return null;
  }
  if (!exch) return null;
  const present = candSignals.reduce((a, x) => a + x.n, 0);
  if (totalH && present > totalH - exch + 0.5) return null;
  const tol = TOL['1H'];
  const free = exp.filter((e) => e.mult === 's' && e.delta > 0.5 && e.delta < 13.5 && candSignals.every((x) => Math.abs(x.delta - e.delta) > 2 * tol));
  if (!free.length) return null;
  const target = totalH ? exch / totalH : 0.1;
  const best = free.reduce((a, e) => (Math.abs(e.w - target) < Math.abs(a.w - target) ? e : a));
  return [...candSignals, { delta: best.delta + 0.4 * tol, n: exch, nb: 0, note: 'OH/NH intercambiable (posible)' }];
}

async function attachNames(items) {
  const keys = [...new Set(items.filter((i) => i.inchikey && !i.name).map((i) => i.inchikey))];
  if (!keys.length) return;
  try {
    const props = await pubchem.propertiesByInchikeys(keys);
    const byKey = new Map();
    for (const p of props) if (!byKey.has(p.inchikey)) byKey.set(p.inchikey, p);
    for (const i of items) {
      const p = byKey.get(i.inchikey);
      if (p) Object.assign(i, { name: i.name || p.title, iupac: i.iupac || p.iupac, cid: i.cid || p.cid });
    }
  } catch {
    /* los nombres son opcionales */
  }
}

async function mapLimit(items, limit, fn) {
  const out = new Array(items.length);
  let k = 0;
  await Promise.all(
    Array.from({ length: Math.min(limit, items.length) }, async () => {
      while (k < items.length) {
        const i = k++;
        out[i] = await fn(items[i]).catch(() => null);
      }
    }),
  );
  return out;
}

export async function findCandidates(analysis, { formula = null, rangeOnly = false, maxResults = 15, maxPredict = 12 } = {}) {
  const nucleus = analysis.nucleus;
  const exp = expSignals(analysis);
  // Sólo cuando el usuario fijó los nH de cada multiplete se compara el nº total de H
  const expTotalH = nucleus === '1H' && analysis.userIntegrals ? analysis.signals.reduce((a, s) => a + (s.nH || 0), 0) : null;
  const warnings = [];
  const results = [];
  const score = (signals, f, smiles) => {
    const totalH = hCount(f);
    const base = concordance(exp, signals, nucleus, { rangeOnly, totalH, expTotalH });
    return { base, totalH, smiles, signals };
  };
  // 2.ª pasada: multiplicidades esperadas (desde la estructura) y OH/NH intercambiables
  const refine = (x) => {
    let signals = x.signals;
    let base = x.base;
    if (nucleus === '1H' && x.id && signals.some((s) => s.atoms?.length)) {
      try {
        const mol = OCL.Molecule.fromMolfile(getMolfile(x.id));
        mol.ensureHelperArrays(OCL.Molecule.cHelperNeighbours);
        signals = signals.map((s) => ({ ...s, nb: s.atoms?.length ? vicinalH(mol, s.atoms[0]) : null }));
        base = concordance(exp, signals, nucleus, { rangeOnly, totalH: x.totalH, expTotalH }) || base;
      } catch {
        /* sin multiplicidades */
      }
    }
    const ext = withExchangeable(exp, signals, x.smiles, nucleus, x.totalH);
    if (!ext) return base;
    const c2 = concordance(exp, ext, nucleus, { rangeOnly, totalH: x.totalH, expTotalH });
    return c2 && c2.score > base.score ? { ...c2, exchangeable: true } : base;
  };

  // 1) NMRShiftDB2
  const all = getAllSpectra(nucleus);
  if (!all) warnings.push('NMRShiftDB2 aún no está disponible (se está descargando o indexando).');
  const dbBest = new Map();
  for (const entry of all || []) {
    if (formula && entry.formula !== formula) continue;
    for (const sp of entry.spectra) {
      const s = score(sp.signals, entry.formula, entry.smiles);
      if (!s.base) continue;
      const value = sp.calculated ? s.base.score * 0.97 : s.base.score;
      const key = entry.inchikey || `id${entry.id}`;
      const prev = dbBest.get(key);
      if (!prev || value > prev.score) {
        dbBest.set(key, {
          ...s.base,
          _s: { ...s, id: entry.id },
          calcFactor: sp.calculated ? 0.97 : 1,
          score: r(value, 1),
          source: sp.calculated ? 'NMRShiftDB2 (espectro calculado)' : 'NMRShiftDB2 (espectro experimental)',
          kind: sp.calculated ? 'calculado' : 'experimental',
          solvent: /unreported|unknown/i.test(sp.solvent || '') ? null : sp.solvent,
          field: sp.field,
          nmrshiftdbId: entry.id,
          url: `https://nmrshiftdb.nmr.uni-koeln.de/portal/js_pane/P-Results/nmrshiftdbaction/showDetailsFromHome/molNumber/${entry.id}`,
          smiles: entry.smiles,
          inchikey: entry.inchikey,
          formula: entry.formula,
          candidateSignals: sp.signals.map((x) => ({ delta: r(x.delta, 2), n: x.n })),
        });
      }
    }
  }
  const firstPass = [...dbBest.values()].sort((a, b) => b.score - a.score).slice(0, 300);
  for (const x of firstPass) {
    const c = refine(x._s);
    Object.assign(x, c, { score: r(c.score * x.calcFactor, 1) });
    delete x._s;
    delete x.calcFactor;
  }
  results.push(...firstPass.sort((a, b) => b.score - a.score).slice(0, maxResults * 2));

  // 2) Predicción para compuestos de PubChem con la fórmula indicada
  let predicted = 0;
  if (formula) {
    try {
      const cids = await pubchem.cidsByFormula(formula, 60);
      // Se excluyen sales/mezclas (".") e isótopos marcados
      const props = (await pubchem.properties(cids)).filter(
        (p) => p.formula === formula && p.smiles && !p.smiles.includes('.') && !/\[2H\]|\[13C\]|\[H\+\]/.test(p.smiles),
      );
      const todo = props.slice(0, maxPredict);
      let limited = false;
      const preds = await mapLimit(todo, 1, async (p) => {
        const mol = parseMolecule({ smiles: p.smiles });
        const key = `v2:${mol.getIDCode()}`;
        let cached = getCachedPrediction(key);
        if (!cached) {
          if (limited) return null;
          let pr;
          try {
            pr = await runPrediction(mol, { retries: 1 });
          } catch (e) {
            if (/Too Many Requests|429/i.test(e.message)) limited = true;
            throw e;
          }
          // nb: nº de H vecinos con J ≥ 3 Hz (los acoplamientos menores no se resuelven en la práctica)
          const nbOf = (s) =>
            (s.js || []).filter((jj) => jj.coupling >= 3).reduce((a, jj) => a + (neighborsFromMultiplicity(jj.multiplicity) ?? 1), 0);
          const simp = (sp) =>
            sp ? sp.joinedSignals.map((s) => ({ delta: s.delta, n: s.nbAtoms, m: s.multiplicity || null, nb: nbOf(s) })) : [];
          cached = { proton: simp(pr.spectra.proton), carbon: simp(pr.spectra.carbon) };
          setCachedPrediction(key, cached.proton, cached.carbon);
        }
        const signals = (nucleus === '1H' ? cached.proton : cached.carbon).map((s) => ({
          delta: s.delta,
          n: s.n,
          nb: nucleus === '1H' ? (s.nb ?? neighborsFromMultiplicity(s.m)) : null,
        }));
        const s = score(signals, p.formula, p.smiles);
        if (!s.base) return null;
        return { ...refine(s), p, signals };
      });
      for (const x of preds.filter(Boolean)) {
        predicted++;
        results.push({
          ...x,
          source: 'Predicción (compuesto de PubChem)',
          kind: 'predicho',
          smiles: x.p.smiles,
          inchikey: x.p.inchikey,
          formula: x.p.formula,
          name: x.p.title,
          iupac: x.p.iupac,
          cid: x.p.cid,
          candidateSignals: x.signals.map((s) => ({ delta: r(s.delta, 2), n: s.n })),
        });
        delete results[results.length - 1].p;
        delete results[results.length - 1].signals;
      }
      if (limited) {
        warnings.push(
          `El servicio de predicción limitó la frecuencia de consultas: se evaluaron ${predicted} de ${todo.length} candidatos de PubChem. Repita la búsqueda en 1–2 minutos (las predicciones ya hechas quedan guardadas).`,
        );
      }
    } catch (e) {
      warnings.push(`PubChem/predicción: ${e.message}`);
    }
  }

  // Unir por compuesto (InChIKey) quedándose con la mejor concordancia
  const merged = new Map();
  for (const x of results) {
    const key = x.inchikey || x.smiles;
    const prev = merged.get(key);
    if (!prev) merged.set(key, { ...x, sources: [x.source] });
    else {
      if (!prev.sources.includes(x.source)) prev.sources.push(x.source);
      if (x.score > prev.score) merged.set(key, { ...prev, ...x, sources: prev.sources });
    }
  }
  const top = [...merged.values()].sort((a, b) => b.score - a.score).slice(0, maxResults);
  await attachNames(top);
  for (const t of top) {
    try {
      t.svg = moleculeSvg(parseMolecule({ smiles: t.smiles }), { width: 240, height: 180, numbers: false, id: `c${Math.random().toString(36).slice(2, 8)}` });
    } catch {
      t.svg = null;
    }
  }
  return {
    candidates: top,
    searched: { nmrshiftdb: dbBest.size, predicted },
    warnings,
    method:
      nucleus === '1H'
        ? 'Concordancia = 45 % cercanía de desplazamientos (tolerancia ±0.12 ppm) + 30 % cobertura de las señales del candidato + 25 % acuerdo de integrales.'
        : 'Concordancia = 50 % cercanía de desplazamientos (±2.5 ppm) + 50 % cobertura de los carbonos del candidato, corregida por el nº de señales.',
  };
}
