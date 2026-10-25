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
import { describeHydrogenCarrier, parseMultiplicity } from './chemistry.js';
import { getCachedPrediction, setCachedPrediction } from '../db.js';
import { hoseReady, predictLocal } from './hose.js';

const r = (v, d = 1) => Math.round(v * 10 ** d) / 10 ** d;
// Criterio de aceptación de desplazamientos: desviación estándar (DE) de los errores
// δ(candidato) − δ(experimental) de las señales asignadas ≤ SD_MAX (ppm).
// La DE no castiga un desplazamiento sistemático común a todas las señales (p. ej. por el disolvente);
// ese sesgo (error medio) se penaliza aparte de forma suave con OFFSET.
export const SD_MAX = { '1H': 0.1, '13C': 1.5 };
export const MIN_SCORE = 90;
const OFFSET = { '1H': 0.25, '13C': 3 };
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

/** Media y desviación estándar muestral (n − 1) de una lista de errores. */
export function errorStats(errs) {
  const n = errs.length;
  if (!n) return { n: 0, mean: null, sd: null, rmse: null };
  const mean = errs.reduce((a, b) => a + b, 0) / n;
  // Con una sola señal la DE no está definida: se usa el valor absoluto del error
  const sd = n > 1 ? Math.sqrt(errs.reduce((a, e) => a + (e - mean) ** 2, 0) / (n - 1)) : Math.abs(errs[0]);
  const rmse = Math.sqrt(errs.reduce((a, e) => a + e * e, 0) / n);
  return { n, mean, sd, rmse };
}

/**
 * Señales del candidato superpuestas: si varias caen dentro del intervalo de un mismo
 * multiplete experimental (p. ej. los 4 H aromáticos de un anillo a 60–90 MHz), se
 * agrupan en una sola señal (suma de H, δ medio ponderado). Además se calcula qué
 * fracción de la intensidad de las líneas del multiplete queda cerca de alguna señal
 * del candidato (`lineCoverage`), para distinguir isómeros con el mismo intervalo.
 */
/**
 * Semejanza de forma (0–1) entre las líneas de un multiplete experimental y las señales
 * del candidato que caen en él: cada señal se dibuja como multiplete de primer orden
 * (regla n+1, J ≈ 7.5 Hz), se alinea por el centro y se mide el solapamiento de las envolventes.
 */
function envelopeOverlap(e, list, freq) {
  if (!e.lines?.length || !freq) return 1;
  const J = 7.5 / freq;
  const sig = Math.max(0.006, 3 / freq);
  const n = list.reduce((a, x) => a + x.n, 0);
  const center = list.reduce((a, x) => a + x.delta * x.n, 0) / n;
  const shift = e.delta - center;
  const sticks = [];
  for (const x of list) {
    const k = Math.min(6, Math.max(0, x.nb ?? 1));
    let c = 1;
    const coef = [];
    for (let i = 0; i <= k; i++) {
      coef.push(c);
      c = (c * (k - i)) / (i + 1);
    }
    const tot = coef.reduce((a, b) => a + b, 0);
    coef.forEach((w, i) => sticks.push({ ppm: x.delta + shift + (i - k / 2) * J, w: (x.n * w) / tot }));
  }
  const lo = e.from - 0.15;
  const hi = e.to + 0.15;
  const step = 0.002;
  let sp = 0;
  let sq = 0;
  const P = [];
  const Q = [];
  for (let v = lo; v <= hi; v += step) {
    let p = 0;
    let q = 0;
    for (const l of e.lines) p += l.intensity * Math.exp(-0.5 * ((v - l.ppm) / sig) ** 2);
    for (const s of sticks) q += s.w * Math.exp(-0.5 * ((v - s.ppm) / sig) ** 2);
    P.push(p);
    Q.push(q);
    sp += p;
    sq += q;
  }
  let overlap = 0;
  for (let i = 0; i < P.length; i++) overlap += Math.min(P[i] / (sp || 1), Q[i] / (sq || 1));
  return overlap;
}

function mergeIntoMultiplets(exp, cand, freq, margin = 0.05) {
  const sumN = cand.reduce((a, x) => a + x.n, 0) || 1;
  const groups = new Map();
  const out = [];
  const pending = [];
  for (const x of cand) {
    if (x.credit) {
      out.push(x);
      continue;
    }
    let best = -1;
    let bd = Infinity;
    exp.forEach((e, i) => {
      if (e.from === undefined || e.to === undefined || e.to - e.from < 0.005) return;
      if (x.delta >= e.from - margin && x.delta <= e.to + margin) {
        const d = Math.abs(x.delta - e.delta);
        if (d < bd) {
          bd = d;
          best = i;
        }
      }
    });
    if (best < 0) out.push(x);
    else pending.push({ x, i: best, d: bd });
  }
  // Sólo se agrupan tantas señales como admita la integral del multiplete (las más cercanas primero)
  pending.sort((a, b) => a.d - b.d);
  const used = new Map();
  for (const p of pending) {
    const cap = exp[p.i].w * 1.35 + 0.02;
    const now = used.get(p.i) || 0;
    if (now + p.x.n / sumN <= cap || !groups.has(p.i)) {
      if (!groups.has(p.i)) groups.set(p.i, []);
      groups.get(p.i).push(p.x);
      used.set(p.i, now + p.x.n / sumN);
    } else out.push(p.x);
  }
  for (const [i, list] of groups) {
    if (list.length === 1) {
      // Una sola señal frente a un multiplete ancho (> ~25 Hz): también se compara la forma
      const e = exp[i];
      const widthHz = (e.to - e.from) * (freq || 400);
      if (widthHz > 25 && e.lines?.length > 3) {
        const shape = envelopeOverlap(e, list, freq);
        out.push({ ...list[0], shape, note: `${list[0].note ? `${list[0].note}; ` : ''}semejanza de forma con el multiplete ${Math.round(shape * 100)} %` });
      } else out.push(list[0]);
      continue;
    }
    const e = exp[i];
    const n = list.reduce((a, x) => a + x.n, 0);
    const delta = list.reduce((a, x) => a + x.delta * x.n, 0) / n;
    const shape = envelopeOverlap(e, list, freq);
    out.push({
      delta,
      n,
      nb: null,
      shape,
      merged: list.map((x) => x.delta),
      note: `${list.length} señales del candidato (${list.map((x) => x.delta.toFixed(2)).join(', ')}) superpuestas en el multiplete; semejanza de forma ${Math.round(shape * 100)} %`,
    });
  }
  return out;
}

/**
 * @param exp [{delta, w}] señales experimentales (w = fracción de la integral)
 * @param cand [{delta, n}] señales del candidato (n = nº de H o de C)
 * @param sdMax desviación estándar máxima aceptada para los errores (ppm)
 */
export function concordance(exp, cand, nucleus, { rangeOnly = false, totalH = null, expTotalH = null, sdMax = SD_MAX[nucleus] } = {}) {
  let c = cand;
  if (rangeOnly) {
    const lo = Math.min(...exp.map((e) => e.delta)) - MARGIN[nucleus];
    const hi = Math.max(...exp.map((e) => e.delta)) + MARGIN[nucleus];
    c = cand.filter((x) => x.delta >= lo && x.delta <= hi);
  }
  if (!c.length || !exp.length) return null;
  if (nucleus === '1H') c = mergeIntoMultiplets(exp, c, exp[0]?.freq);
  const sumN = c.reduce((a, x) => a + x.n, 0) || 1;
  const f = c.map((x) => x.n / sumN);
  const off = OFFSET[nucleus];
  const gauss = (d, s) => Math.exp(-0.5 * (d / s) ** 2);

  // Asignación uno a uno (voraz): cada señal experimental a una señal distinta del candidato.
  // El coste combina la diferencia de desplazamiento con el desacuerdo de multiplicidad e integración,
  // para no intercambiar señales cercanas (p. ej. H orto/para de un anillo).
  const pairs = [];
  for (let i = 0; i < exp.length; i++) {
    for (let j = 0; j < c.length; j++) {
      // Un OH/NH intercambiable sólo puede corresponder a un singlete con integral compatible
      if (c[j].credit && nucleus === '1H' && (exp[i].mult !== 's' || exp[i].w > 2 * f[j] + 0.05)) continue;
      let cost = Math.abs(exp[i].delta - c[j].delta) / off;
      if (exp[i].nb !== null && exp[i].nb !== undefined && c[j].nb !== null && c[j].nb !== undefined) cost += 0.6 * Math.min(2, Math.abs(exp[i].nb - c[j].nb));
      if (nucleus === '1H') cost += 2 * Math.abs(exp[i].w - f[j]);
      pairs.push([i, j, cost]);
    }
  }
  pairs.sort((a, b) => a[2] - b[2]);
  const ei = new Array(exp.length).fill(-1);
  const cj = new Array(c.length).fill(-1);
  for (const [i, j] of pairs) {
    if (ei[i] < 0 && cj[j] < 0) {
      ei[i] = j;
      cj[j] = i;
    }
  }

  // Errores de las señales asignadas (los OH/NH intercambiables supuestos no tienen δ predicho: se excluyen)
  const errs = [];
  exp.forEach((e, i) => {
    if (ei[i] >= 0 && !c[ei[i]].credit) errs.push(c[ei[i]].delta - e.delta);
  });
  const st = errorStats(errs);
  const sdOk = st.sd !== null && st.sd <= sdMax + 1e-9;
  // Similitud de cada señal: residuo respecto del error medio. La aceptación estricta es DE ≤ sdMax;
  // dentro del % se usa una escala más suave (2·sdMax) para que, entre candidatos aceptados,
  // pese más la forma de los multipletes superpuestos (x.shape)
  const simPair = (x, e) => (x.credit ? x.credit : gauss(x.delta - e.delta - (st.mean ?? 0), 2 * sdMax) * (x.shape ?? 1));
  const offsetFactor = st.mean === null ? 1 : gauss(st.mean, off);
  const nearest = (d0, list) => list.reduce((m, x) => Math.min(m, Math.abs(x.delta - d0)), Infinity);

  let shift = 0;
  let integErr = 0;
  let multSum = 0;
  let multN = 0;
  const matches = [];
  exp.forEach((e, i) => {
    if (ei[i] >= 0) {
      const x = c[ei[i]];
      shift += e.w * simPair(x, e);
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
        diff: x.credit ? null : r(x.delta - e.delta, 3),
        n: x.n,
        expectedNeighbors: x.nb ?? null,
        multiplicityOk: multOk,
        note: x.note || null,
      });
    } else {
      // Señal sobrante: crédito parcial si cae cerca de una señal ya asignada (multiplete desdoblado)
      shift += 0.4 * e.w * gauss(nearest(e.delta, c), off);
      integErr += e.w;
      matches.push({ experimental: e.delta, candidate: null, diff: null, n: null, note: 'sin señal propia en el candidato' });
    }
  });
  shift *= offsetFactor;
  let coverage = 0;
  c.forEach((x, j) => {
    if (cj[j] >= 0) coverage += f[j] * simPair(x, exp[cj[j]]);
    else {
      coverage += 0.4 * f[j] * gauss(nearest(x.delta, exp), off);
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
  // Toda señal experimental debe tener una señal propia en el candidato: las que quedan sin
  // explicar (p. ej. un OH en una estructura sin OH) reducen la concordancia
  const expExplained = ei.filter((j) => j >= 0).length / exp.length;
  total *= 0.85 + 0.15 * expExplained;
  return {
    score: r(total * 100, 1),
    expExplained: r(expExplained * 100, 1),
    shift: r(shift * 100, 1),
    coverage: r(coverage * 100, 1),
    integration: integration === null ? null : r(integration * 100, 1),
    signalCount: r(count * 100, 1),
    hAgreement: hAgreement === null ? null : r(hAgreement * 100, 1),
    multiplicity: multiplicity === null ? null : r(multiplicity * 100, 1),
    errorSd: st.sd === null ? null : r(st.sd, 3),
    errorMean: st.mean === null ? null : r(st.mean, 3),
    errorRmse: st.rmse === null ? null : r(st.rmse, 3),
    nErrors: st.n,
    sdOk,
    unexplained: c.filter((x, j) => cj[j] < 0).map((x) => r(x.delta, 2)),
    matches,
  };
}

function expSignals(analysis) {
  if (analysis.nucleus === '13C') {
    const n = analysis.signals.length;
    return analysis.signals.map((s) => ({ delta: s.delta, w: 1 / n }));
  }
  return analysis.signals.map((s) => ({
    delta: s.delta,
    w: s.fraction,
    mult: s.multiplicity,
    nb: neighborsFromMultiplicity(s.multiplicity),
    from: s.from,
    to: s.to,
    lines: s.lines,
    freq: analysis.frequency,
  }));
}

const hCount = (formula) => Number(String(formula || '').match(/H(\d*)(?![a-z])/)?.[1] || (/H(?![a-z])/.test(formula || '') ? 1 : 0));

/**
 * Los H de OH/NH/SH no suelen aparecer en las bases de datos ni en la predicción.
 * Si el candidato los tiene y el espectro muestra un singlete sin asignar, se
 * considera una posible señal intercambiable (con crédito parcial).
 */
function withExchangeable(exp, candSignals, smiles, nucleus, totalH) {
  const tol = OFFSET['1H'];
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
  // Intervalos plausibles según el tipo de OH/NH del candidato (alcohol, fenol, ácido, amina, amida…)
  const ranges = [];
  try {
    const mol = parseMolecule({ smiles });
    for (let a = 0; a < mol.getAtoms(); a++) {
      if ([7, 8, 16].includes(mol.getAtomicNo(a)) && mol.getAllHydrogens(a)) {
        const rg = describeHydrogenCarrier(mol, a).range;
        if (rg) ranges.push(rg);
      }
    }
  } catch {
    /* sin intervalos */
  }
  const plausible = (d) => !ranges.length || ranges.some(([lo, hi]) => d >= lo - 0.5 && d <= hi + 0.5);
  const free = exp.filter(
    (e) => e.mult === 's' && e.delta > 0.5 && e.delta < 13.5 && plausible(e.delta) && candSignals.every((x) => Math.abs(x.delta - e.delta) > tol),
  );
  if (!free.length) return null;
  const target = totalH ? exch / totalH : 0.1;
  // El singlete elegido debe tener una integral compatible con los H intercambiables (no un CH₃, p. ej.)
  const compatible = free.filter((e) => e.w <= 2 * target + 0.05);
  if (!compatible.length) return null;
  const best = compatible.reduce((a, e) => (Math.abs(e.w - target) < Math.abs(a.w - target) ? e : a));
  // Crédito parcial (85 %): su posición depende de concentración y disolvente
  return [...candSignals, { delta: best.delta, n: exch, nb: 0, credit: 0.85, note: 'OH/NH intercambiable (posible)' }];
}

async function attachNames(items) {
  const keys = [...new Set(items.filter((i) => i.inchikey && !i.name).map((i) => i.inchikey))];
  if (!keys.length) return;
  try {
    const props = await pubchem.propertiesByInchikeys(keys);
    // Para un mismo InChIKey PubChem puede devolver varios registros: el de menor CID es el compuesto principal
    const byKey = new Map();
    for (const p of [...props].sort((a, b) => a.cid - b.cid)) if (!byKey.has(p.inchikey)) byKey.set(p.inchikey, p);
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

export async function findCandidates(analysis, { formula = null, rangeOnly = false, maxResults = 15, maxPredict = 6, maxPubchem = 150, sdMax = null, minScore = MIN_SCORE } = {}) {
  const nucleus = analysis.nucleus;
  const sdLimit = Number(sdMax) > 0 ? Number(sdMax) : SD_MAX[nucleus];
  const exp = expSignals(analysis);
  // Sólo cuando el usuario fijó los nH de cada multiplete se compara el nº total de H
  const expTotalH = nucleus === '1H' && analysis.userIntegrals ? analysis.signals.reduce((a, s) => a + (s.nH || 0), 0) : null;
  const warnings = [];
  const results = [];
  const score = (signals, f, smiles) => {
    const totalH = hCount(f);
    const base = concordance(exp, signals, nucleus, { rangeOnly, totalH, expTotalH, sdMax: sdLimit });
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
        signals = signals.map((s) => {
          if (!s.atoms?.length) return { ...s, nb: null };
          const atoms = s.atoms.filter((a) => a < mol.getAtoms());
          if (!atoms.length) return { ...s, nb: null };
          // Nº real de H según la estructura (algunos espectros listan un solo H por CH₃/CH₂)
          const nH = atoms.reduce((acc, a) => acc + mol.getAllHydrogens(a), 0) || s.n;
          // OH/NH/SH medidos: posición variable → intercambiables (crédito parcial, fuera de la DE)
          if ([7, 8, 16].includes(mol.getAtomicNo(atoms[0]))) return { ...s, n: nH, nb: 0, credit: 0.85, note: 'OH/NH intercambiable (dato de la base)' };
          return { ...s, n: nH, nb: vicinalH(mol, atoms[0]) };
        });
        base = concordance(exp, signals, nucleus, { rangeOnly, totalH: x.totalH, expTotalH, sdMax: sdLimit }) || base;
      } catch {
        /* sin multiplicidades */
      }
    }
    const ext = withExchangeable(exp, signals, x.smiles, nucleus, x.totalH);
    if (!ext) return base;
    const c2 = concordance(exp, ext, nucleus, { rangeOnly, totalH: x.totalH, expTotalH, sdMax: sdLimit });
    return c2 && c2.score > base.score ? { ...c2, exchangeable: true } : base;
  };

  // Predicción completa NMRium (con J) para un candidato, con caché persistente.
  // Si el servicio limita la frecuencia de consultas se deja de llamar durante esta búsqueda.
  const nmriumState = { limited: false, used: 0 };
  const nmriumScore = async (p) => {
    let mol;
    try {
      mol = parseMolecule({ smiles: p.smiles });
    } catch {
      return null;
    }
    const key = `v2:${mol.getIDCode()}`;
    let cached = getCachedPrediction(key);
    if (!cached) {
      if (nmriumState.limited) return null;
      let pr;
      try {
        pr = await runPrediction(mol, { retries: 0 });
      } catch (e) {
        if (/Too Many Requests|429/i.test(e.message)) nmriumState.limited = true;
        return null;
      }
      // nb: nº de H vecinos con J ≥ 3 Hz (los acoplamientos menores no se resuelven en la práctica)
      const nbOf = (s) => (s.js || []).filter((jj) => jj.coupling >= 3).reduce((a, jj) => a + (neighborsFromMultiplicity(jj.multiplicity) ?? 1), 0);
      const simp = (sp) => (sp ? sp.joinedSignals.map((s) => ({ delta: s.delta, n: s.nbAtoms, m: s.multiplicity || null, nb: nbOf(s) })) : []);
      cached = { proton: simp(pr.spectra.proton), carbon: simp(pr.spectra.carbon) };
      setCachedPrediction(key, cached.proton, cached.carbon);
    }
    nmriumState.used++;
    const signals = (nucleus === '1H' ? cached.proton : cached.carbon).map((s) => ({
      delta: s.delta,
      n: s.n,
      nb: nucleus === '1H' ? (s.nb ?? neighborsFromMultiplicity(s.m)) : null,
    }));
    const s = score(signals, p.formula, p.smiles);
    return s.base ? { ...refine(s), p, signals, method: 'nmrium' } : null;
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

  // Sin fórmula: los mejores candidatos de NMRShiftDB2 también se evalúan con el predictor local
  // y con NMRium (en CDCl₃ y con J), porque la base puede tener el espectro en otro disolvente
  if (!formula && hoseReady()) {
    for (const x of firstPass.slice(0, 40)) {
      let sig;
      try {
        sig = predictLocal(parseMolecule({ smiles: x.smiles }), nucleus);
      } catch {
        continue;
      }
      if (!sig?.length) continue;
      const s = score(sig, x.formula, x.smiles);
      if (!s.base) continue;
      results.push({
        ...refine(s),
        source: 'Predicción local (entornos HOSE de NMRShiftDB2)',
        kind: 'predicho',
        smiles: x.smiles,
        inchikey: x.inchikey,
        formula: x.formula,
        candidateSignals: sig.map((q) => ({ delta: r(q.delta, 2), n: q.n })),
      });
    }
  }
  if (!formula && nucleus === '1H') {
    for (const x of firstPass.slice(0, 8)) {
      const y = await nmriumScore({ smiles: x.smiles, formula: x.formula, inchikey: x.inchikey });
      if (!y) continue;
      const item = { ...y, source: 'Predicción NMRium (con constantes J)', kind: 'predicho', smiles: x.smiles, inchikey: x.inchikey, formula: x.formula };
      item.candidateSignals = y.signals.map((s) => ({ delta: r(s.delta, 2), n: s.n }));
      delete item.p;
      delete item.signals;
      delete item.method;
      results.push(item);
    }
  }

  // 2) Compuestos de PubChem con la fórmula indicada:
  //    predictor local (HOSE/NMRShiftDB2) para todos los isómeros y NMRium para refinar los mejores
  let predicted = 0;
  if (formula) {
    try {
      const cids = await pubchem.cidsByFormula(formula, maxPubchem);
      // Se excluyen sales/mezclas (".") e isótopos marcados
      const props = (await pubchem.properties(cids)).filter(
        (p) => p.formula === formula && p.smiles && !p.smiles.includes('.') && !/\[\d+[A-Z]|\[H\+\]/.test(p.smiles),
      );
      const local = hoseReady();
      if (!local) warnings.push('El predictor local se está construyendo (unos minutos tras el primer arranque); por ahora sólo se usa el servicio NMRium, que limita el nº de consultas.');
      const scored = [];
      if (local) {
        for (const p of props) {
          let sig;
          try {
            sig = predictLocal(parseMolecule({ smiles: p.smiles }), nucleus);
          } catch {
            continue;
          }
          if (!sig?.length) continue;
          const s = score(sig, p.formula, p.smiles);
          if (s.base) scored.push({ ...refine(s), p, signals: sig, method: 'local' });
        }
        scored.sort((a, b) => b.score - a.score);
      }
      // NMRium: predicción completa (con J) para los mejores candidatos, si el servicio responde
      const remoteTodo = local ? scored.slice(0, 3).map((x) => x.p) : props.slice(0, maxPredict);
      const remote = [];
      for (const p of remoteTodo) {
        const x = await nmriumScore(p);
        if (x) remote.push(x);
      }
      const limited = nmriumState.limited;
      predicted = scored.length + remote.length;
      for (const x of [...scored.slice(0, maxResults * 2), ...remote]) {
        const item = {
          ...x,
          source: x.method === 'local' ? 'Predicción local (entornos HOSE de NMRShiftDB2)' : 'Predicción NMRium (con constantes J)',
          kind: 'predicho',
          smiles: x.p.smiles,
          inchikey: x.p.inchikey,
          formula: x.p.formula,
          name: x.p.title,
          iupac: x.p.iupac,
          cid: x.p.cid,
          candidateSignals: x.signals.map((s) => ({ delta: r(s.delta, 2), n: s.n })),
        };
        delete item.p;
        delete item.signals;
        delete item.method;
        results.push(item);
      }
      if (limited && !local) {
        warnings.push(
          `El servicio de predicción NMRium limitó la frecuencia de consultas: se evaluaron ${remote.length} candidatos. Repita la búsqueda en 1–2 minutos (las predicciones ya hechas quedan guardadas).`,
        );
      }
      if (props.length >= maxPubchem) warnings.push(`Se evaluaron los primeros ${maxPubchem} isómeros de PubChem con esta fórmula.`);
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
  const ranked = [...merged.values()].sort((a, b) => b.score - a.score);
  // Sólo se muestran los candidatos con concordancia ≥ minScore
  // Se muestran sólo los candidatos con DE de los errores ≤ sdLimit y concordancia ≥ minScore
  const passes = (x) => x.sdOk && x.score >= minScore;
  const top = ranked.filter(passes).slice(0, maxResults);
  // Si ninguno supera el umbral se informa del mejor candidato encontrado
  const bestBelow = !top.length && ranked[0]
    ? { score: ranked[0].score, errorSd: ranked[0].errorSd, name: ranked[0].name || null, smiles: ranked[0].smiles, inchikey: ranked[0].inchikey }
    : null;
  await attachNames(bestBelow ? [...top, bestBelow] : top);
  for (const t of top) {
    try {
      t.svg = moleculeSvg(parseMolecule({ smiles: t.smiles }), { width: 240, height: 180, numbers: false, id: `c${Math.random().toString(36).slice(2, 8)}` });
    } catch {
      t.svg = null;
    }
  }
  return {
    candidates: top,
    params: { sdMax: sdLimit, minScore },
    evaluated: ranked.length,
    belowThreshold: ranked.length - ranked.filter(passes).length,
    bestBelow,
    searched: { nmrshiftdb: dbBest.size, predicted },
    warnings,
    method:
      nucleus === '1H'
        ? `Criterio de desplazamientos: la desviación estándar de los errores δ(candidato) − δ(experimental) de las señales asignadas debe ser ≤ ${sdLimit} ppm. Concordancia = 35 % desplazamientos + 20 % cobertura de las señales del candidato + 15 % integración + 10 % nº de señales + 20 % multiplicidad (regla n+1). Se muestran los candidatos que cumplen la DE y tienen concordancia ≥ ${minScore} %.`
        : `Criterio de desplazamientos: la desviación estándar de los errores δ(candidato) − δ(experimental) debe ser ≤ ${sdLimit} ppm. Concordancia = 50 % desplazamientos + 50 % cobertura de los carbonos del candidato, corregida por el nº de señales. Se muestran los candidatos que cumplen la DE y tienen concordancia ≥ ${minScore} %.`,
  };
}
