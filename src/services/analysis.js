// Núcleo del análisis: predicción, interpretación, búsqueda en bases de datos,
// comparación predicho/experimental y guardado de los espectros.
import { predict } from 'nmr-processing';
import { OCL, atomMapping, moleculeInfo, moleculeSvg, parseMolecule } from './molecule.js';
import {
  SOLVENTS,
  carbonRegion,
  deptPhase,
  describeCarbon,
  describeHydrogenCarrier,
  multiplicityName,
  parseMultiplicity,
  protonRegion,
} from './chemistry.js';
import { carbonFrequency, pickLines, protonRange, simulateProton } from './spectrum.js';
import * as pubchem from './pubchem.js';
import * as nmrshiftdb from './nmrshiftdb.js';
import { getCompoundSpectra, upsertCompound, upsertSpectrum } from '../db.js';

const r2 = (v, d = 2) => Math.round(v * 10 ** d) / 10 ** d;
const predictionCache = new Map();

async function runPrediction(mol) {
  const key = mol.getIDCode();
  if (predictionCache.has(key)) return predictionCache.get(key);
  const p = Promise.race([
    predict(mol, { predictOptions: { H: {}, C: {} } }),
    new Promise((_, rej) => setTimeout(() => rej(new Error('La predicción tardó demasiado')), 90000)),
  ]);
  predictionCache.set(key, p);
  if (predictionCache.size > 100) predictionCache.delete(predictionCache.keys().next().value);
  try {
    return await p;
  } catch (e) {
    predictionCache.delete(key);
    throw e;
  }
}

// Desplazamientos aproximados para H intercambiables (no los da el predictor)
const EXCHANGEABLE_DELTA = { 'H-acido': 11.5, 'H-fenol': 5.5, 'H-alcohol': 2.0, 'H-amina': 1.5, 'H-amida': 6.5, 'H-tiol': 1.6, 'H-otro': 2.0 };

/** Sencilla relación de números enteros para las integrales. */
function integrationRatio(values) {
  const gcd = (a, b) => (b ? gcd(b, a % b) : a);
  const g = values.reduce((a, v) => gcd(a, v), 0) || 1;
  return values.map((v) => v / g);
}

function protonAnalysis(mol, prediction, frequency) {
  const mh = prediction.moleculeWithHydrogens;
  const pred = prediction.spectra.proton;
  const heavyOf = (h) => mh.getConnAtom(h, 0);
  const joined = pred.joinedSignals;
  const byDia = new Map();
  for (const s of joined) for (const d of s.diaIDs || []) byDia.set(d, s);

  const signals = joined.map((s) => {
    const parents = [...new Set(s.atoms.map(heavyOf))];
    const env = describeHydrogenCarrier(mol, parents[0]);
    const js = (s.js || []).map((j) => {
      const partner = (j.diaIDs || []).map((d) => byDia.get(d)).find(Boolean);
      const n = parseMultiplicity(j.multiplicity)[0]?.n;
      return {
        J: r2(j.coupling, 2),
        multiplicity: j.multiplicity,
        bonds: j.pathLength,
        nPartners: n,
        partner: partner
          ? { delta: r2(partner.delta, 2), atoms: [...new Set(partner.atoms.map(heavyOf))].map((a) => a + 1) }
          : null,
      };
    });
    return {
      delta: r2(s.delta, 3),
      multiplicity: s.multiplicity || 's',
      multiplicityName: multiplicityName(s.multiplicity || 's'),
      nH: s.nbAtoms,
      atoms: parents.map((a) => a + 1),
      hAtoms: s.atoms,
      js,
      environment: env.text,
      envCode: env.code,
      typicalRange: env.range,
      inTypicalRange: env.range ? s.delta >= env.range[0] - 0.3 && s.delta <= env.range[1] + 0.3 : null,
      region: protonRegion(s.delta),
      exchangeable: false,
      estimated: false,
    };
  });

  // Añadir H intercambiables (OH/NH/SH) no predichos
  const predictedH = new Set(joined.flatMap((s) => s.atoms));
  const extraRaw = [];
  const groups = new Map();
  for (let h = mh.getAtoms(); h < mh.getAllAtoms(); h++) {
    if (mh.getAtomicNo(h) !== 1 || predictedH.has(h)) continue;
    const parent = heavyOf(h);
    if (mh.getAtomicNo(parent) === 6) continue;
    if (!groups.has(parent)) groups.set(parent, []);
    groups.get(parent).push(h);
  }
  for (const [parent, hs] of groups) {
    const env = describeHydrogenCarrier(mol, parent);
    const delta = EXCHANGEABLE_DELTA[env.code] ?? 2;
    signals.push({
      delta,
      multiplicity: 'br s',
      multiplicityName: 'singlete ancho',
      nH: hs.length,
      atoms: [parent + 1],
      hAtoms: hs,
      js: [],
      environment: env.text,
      envCode: env.code,
      typicalRange: env.range,
      inTypicalRange: true,
      region: protonRegion(delta),
      exchangeable: true,
      estimated: true,
    });
    for (const h of hs) extraRaw.push({ delta, atoms: [h], nbAtoms: 1, js: [] });
  }

  signals.sort((a, b) => b.delta - a.delta);
  const ratio = integrationRatio(signals.map((s) => s.nH));
  signals.forEach((s, i) => {
    s.relativeIntegration = ratio[i];
    s.explanation = explainProton(s);
  });

  const rawSignals = [...pred.signals, ...extraRaw].map((s) => ({
    delta: s.delta,
    atoms: s.atoms,
    nbAtoms: s.nbAtoms ?? 1,
    js: (s.js || []).map((j) => ({ coupling: j.coupling, atoms: j.atoms, multiplicity: j.multiplicity })),
  }));
  const range = protonRange(signals.map((s) => s.delta));
  const xy = simulateProton(rawSignals, { frequency, ...range });
  const lines = pickLines(xy, frequency);
  for (const l of lines) {
    let best = null;
    for (const s of signals) {
      const d = Math.abs(s.delta - l.ppm);
      if (!best || d < best.d) best = { d, s };
    }
    l.signal = best ? r2(best.s.delta, 2) : null;
  }
  const integrals = signals.map((s) => {
    const width = (s.js.reduce((a, j) => a + j.J * (j.nPartners || 1), 0) / 2 + 2.5) / frequency;
    return { from: r2(s.delta - width, 4), to: r2(s.delta + width, 4), nH: s.nH };
  });
  return { signals, lines, rawSignals, range, integrals };
}

function explainProton(s) {
  const parts = [];
  if (s.exchangeable) {
    parts.push(
      `${s.environment}: señal ancha y de posición variable (depende de concentración, temperatura y disolvente); desaparece al agitar con D₂O. Desplazamiento estimado.`,
    );
    return parts.join(' ');
  }
  if (!s.js.length) {
    parts.push(`Singlete: los ${s.nH} H equivalentes no tienen H vecinos acoplados a 2–3 enlaces (o los vecinos no portan H).`);
  } else {
    const desc = s.js.map((j) => {
      const who = j.partner ? `H de C${j.partner.atoms.join('/C')} (δ ${j.partner.delta})` : 'H vecinos';
      const n = j.nPartners;
      return `${n ? `${n} ` : ''}${who} con ${j.bonds ? `${j.bonds}J` : 'J'} = ${j.J} Hz`;
    });
    parts.push(`${cap(s.multiplicityName)} por acoplamiento con ${desc.join('; ')}.`);
    const n = s.js.reduce((a, j) => a + (j.nPartners || 0), 0);
    if (s.js.length === 1 && n) parts.push(`Regla n+1: ${n} H vecinos → ${n + 1} líneas.`);
  }
  parts.push(`Entorno: ${s.environment}.`);
  if (s.typicalRange)
    parts.push(
      `Rango típico ${s.typicalRange[0]}–${s.typicalRange[1]} ppm${s.inTypicalRange ? ' (coincide)' : ' (fuera de rango: influencia de sustituyentes vecinos)'}.`,
    );
  return parts.join(' ');
}

const cap = (s) => s.charAt(0).toUpperCase() + s.slice(1);

function carbonAnalysis(mol, prediction) {
  const signals = prediction.spectra.carbon.joinedSignals.map((s) => {
    const a = s.atoms[0];
    const env = describeCarbon(mol, a);
    const nH = mol.getAllHydrogens(a);
    const dept = deptPhase(nH);
    return {
      delta: r2(s.delta, 2),
      nC: s.nbAtoms,
      atoms: s.atoms.map((x) => x + 1),
      type: env.label,
      nH,
      environment: env.text,
      envCode: env.code,
      typicalRange: env.range13C,
      inTypicalRange: env.range13C ? s.delta >= env.range13C[0] - 5 && s.delta <= env.range13C[1] + 5 : null,
      region: carbonRegion(s.delta),
      dept135: dept.dept135,
      dept90: dept.dept90,
      // Intensidad aproximada: proporcional al nº de C equivalentes; los C cuaternarios
      // aparecen más débiles (sin NOE y con T1 largos).
      rawIntensity: s.nbAtoms * (nH === 0 ? 0.45 : 1),
      stats: s.statistics ? { n: s.statistics.nb, min: s.statistics.min, max: s.statistics.max } : null,
    };
  });
  const max = Math.max(...signals.map((s) => s.rawIntensity), 1);
  for (const s of signals) {
    s.intensity = r2((s.rawIntensity / max) * 100, 1);
    delete s.rawIntensity;
    s.explanation =
      `${s.type} (${s.dept135} en DEPT-135; ${s.dept90} en DEPT-90). ${s.environment}.` +
      (s.nC > 1 ? ` Representa ${s.nC} carbonos equivalentes por simetría.` : '') +
      (s.typicalRange ? ` Rango típico ${s.typicalRange[0]}–${s.typicalRange[1]} ppm.` : '');
  }
  signals.sort((a, b) => b.delta - a.delta);
  return { signals };
}

function literatureString(nucleus, freq, solvent, signals) {
  const solv = solvent && solvent !== 'ninguno' ? solvent : 'CDCl3';
  if (nucleus === '1H') {
    const items = signals.map((s) => {
      const J = s.js.length ? `, J = ${[...new Set(s.js.map((j) => j.J.toFixed(1)))].join(', ')} Hz` : '';
      return `${s.delta.toFixed(2)} (${s.multiplicity}${J}, ${s.nH}H)`;
    });
    return `¹H RMN (${freq} MHz, ${solv}) δ ${items.join(', ')}.`;
  }
  return `¹³C RMN (${Math.round(freq)} MHz, ${solv}) δ ${signals.map((s) => s.delta.toFixed(1)).join(', ')}.`;
}

function functionalGroups(carbonSignals, protonSignals) {
  const codes = new Set([...carbonSignals.map((s) => s.envCode), ...protonSignals.map((s) => s.envCode)]);
  const names = {
    aromatico: 'Anillo aromático', cetona: 'Cetona', aldehido: 'Aldehído', acido: 'Ácido carboxílico',
    ester: 'Éster', amida: 'Amida', haluro_acilo: 'Haluro de acilo', anhidrido: 'Anhídrido', carbonato: 'Carbonato',
    nitrilo: 'Nitrilo', alquino: 'Alquino', alqueno: 'Alqueno', imina: 'Imina', 'C-X': 'Haluro de alquilo',
    'H-alcohol': 'Alcohol', 'H-fenol': 'Fenol', 'H-amina': 'Amina', 'C-S': 'Compuesto de azufre (tioéter/tiol)',
  };
  const out = [];
  for (const [k, v] of Object.entries(names)) if (codes.has(k)) out.push(v);
  if (codes.has('C-O') && !codes.has('H-alcohol') && !codes.has('ester') && !codes.has('acido')) out.push('Éter');
  return out;
}

// ---------- Datos experimentales ----------

function normalizePeaks(peaks) {
  const max = Math.max(...peaks.map((p) => p.intensity), 1);
  return peaks.map((p) => ({ ppm: p.ppm, intensity: r2((p.intensity / max) * 100, 2) })).sort((a, b) => b.ppm - a.ppm);
}

/** Emparejamiento uno a uno por cercanía (greedy) entre desplazamientos predichos y experimentales. */
function matchShifts(pred, exp) {
  const pairs = [];
  for (const p of pred) for (const e of exp) pairs.push({ p, e, d: Math.abs(p.delta - e.delta) });
  pairs.sort((a, b) => a.d - b.d);
  const usedP = new Set();
  const usedE = new Set();
  const out = [];
  for (const pr of pairs) {
    if (usedP.has(pr.p) || usedE.has(pr.e)) continue;
    usedP.add(pr.p);
    usedE.add(pr.e);
    out.push(pr);
  }
  return out;
}

function comparisonSummary(rows) {
  if (!rows.length) return null;
  const devs = rows.map((r) => Math.abs(r.diff));
  return {
    n: rows.length,
    mae: r2(devs.reduce((a, b) => a + b, 0) / devs.length, 3),
    max: r2(Math.max(...devs), 3),
    rows: rows.sort((a, b) => b.experimental - a.experimental),
  };
}

function compareToPrediction(nucleus, expSignals, predicted, mapping) {
  if (!predicted) return null;
  const rows = [];
  if (mapping && expSignals.some((s) => s.atom !== null && s.atom !== undefined)) {
    // Comparación por átomo asignado (NMRShiftDB)
    const predByAtom = new Map();
    for (const s of predicted.signals) for (const a of s.atoms) predByAtom.set(a, s.delta);
    const seen = new Set();
    for (const s of expSignals) {
      if (s.atom === null || s.atom === undefined) continue;
      const key = `${s.atom}:${s.delta}`;
      if (seen.has(key)) continue;
      seen.add(key);
      const ours = mapping[s.atom];
      const pd = ours !== undefined ? predByAtom.get(ours + 1) : undefined;
      if (pd === undefined) continue;
      rows.push({ atom: ours + 1, experimental: s.delta, predicted: pd, diff: r2(pd - s.delta, 3) });
    }
    if (rows.length) return comparisonSummary(rows);
  }
  const exp = expSignals.map((s) => ({ delta: s.delta ?? s.ppm }));
  const pred = predicted.signals.map((s) => ({ delta: s.delta, atoms: s.atoms }));
  for (const m of matchShifts(pred, exp)) {
    if (nucleus === '1H' && m.d > 0.6) continue;
    if (nucleus === '13C' && m.d > 12) continue;
    rows.push({ atom: m.p.atoms.join(','), experimental: m.e.delta, predicted: m.p.delta, diff: r2(m.p.delta - m.e.delta, 3) });
  }
  return comparisonSummary(rows);
}

function nmrshiftdbSpectra(entries, mol) {
  const out = [];
  for (const e of entries) {
    let mapping = null;
    try {
      mapping = atomMapping(OCL.Molecule.fromMolfile(e.molfile), mol);
    } catch {
      mapping = null;
    }
    for (const s of e.spectra) {
      if (s.nucleus !== '1H' && s.nucleus !== '13C') continue;
      // Agrupar señales con el mismo desplazamiento para estimar intensidades
      const grouped = new Map();
      for (const sig of s.signals) {
        const k = sig.delta.toFixed(s.nucleus === '1H' ? 2 : 1);
        if (!grouped.has(k)) grouped.set(k, { delta: sig.delta, mult: sig.mult, atoms: [], count: 0 });
        const g = grouped.get(k);
        g.count++;
        if (!g.mult && sig.mult) g.mult = sig.mult;
        if (sig.atom !== null) g.atoms.push(mapping ? mapping[sig.atom] + 1 : sig.atom + 1);
      }
      const multMap = { Q: 'CH₃', T: 'CH₂', D: 'CH', S: 'C' };
      const signals = [...grouped.values()]
        .map((g) => ({
          delta: g.delta,
          type: s.nucleus === '13C' ? multMap[g.mult] || null : null,
          multiplicity: s.nucleus === '1H' ? g.mult || null : null,
          count: g.count,
          atoms: g.atoms,
        }))
        .sort((a, b) => b.delta - a.delta);
      const peaks = normalizePeaks(signals.map((g) => ({ ppm: g.delta, intensity: g.count })));
      out.push({
        nucleus: s.nucleus,
        kind: s.calculated ? 'calculado' : 'experimental',
        source: 'NMRShiftDB2',
        sourceId: `${e.id}-${s.nucleus}-${s.number}`,
        frequency: s.field || 0,
        solvent: s.solvent,
        signals,
        peaks,
        meta: {
          url: e.url,
          nmrshiftdbId: e.id,
          temperature: s.temperature,
          program: s.program,
          assignment: s.assignment,
          atomMapping: Boolean(mapping),
          rawSignals: s.signals,
          mapping,
          note: s.calculated
            ? 'Espectro calculado (predicción/cálculo cuántico) depositado en NMRShiftDB2'
            : 'Espectro experimental asignado (NMRShiftDB2). Las intensidades se estiman por nº de átomos equivalentes.',
        },
      });
    }
  }
  return out;
}

function pubchemSpectra(nmr) {
  return nmr.spectra.map((s) => ({
    nucleus: s.nucleus,
    kind: 'experimental',
    source: `PubChem · ${s.source}`,
    sourceId: s.sourceId,
    frequency: s.frequency,
    solvent: s.solvent,
    signals: s.peaks.map((p) => ({ delta: p.ppm, intensity: p.intensity })),
    peaks: s.peaks.length ? normalizePeaks(s.peaks) : [],
    meta: { url: s.url, license: s.license, fields: s.fields, thumbnail: s.thumbnail, hasPeaks: s.peaks.length > 0 },
  }));
}

// ---------- Flujo principal ----------

/**
 * Analiza un compuesto completo y guarda los espectros.
 * @param {object} input {smiles|molfile|cid, name, frequency, solvent, mode, query}
 */
export async function analyze(input) {
  const frequency = Math.min(1200, Math.max(40, Number(input.frequency) || 400));
  const solvent = SOLVENTS[input.solvent] ? input.solvent : 'CDCl3';
  const warnings = [];
  let cid = input.cid ? Number(input.cid) : null;
  let props = null;

  if (cid && !input.smiles && !input.molfile) {
    props = (await pubchem.properties([cid]))[0];
    if (!props) throw Object.assign(new Error(`CID ${cid} no encontrado en PubChem`), { status: 404 });
  }
  const mol = parseMolecule({ smiles: input.smiles || props?.smiles, molfile: input.molfile, idcode: input.idcode });
  const info = moleculeInfo(mol);

  // PubChem: identificar el compuesto (nombres, CID) a partir de la estructura
  const pubchemTask = (async () => {
    if (!cid) {
      const cids = await pubchem.cidsBySmiles(info.smiles);
      cid = cids[0] || null;
    }
    if (!cid) return { props: null, synonyms: [], nmr: { spectra: [], links: [] } };
    const [pp, syn, nmr] = await Promise.all([
      props ? Promise.resolve([props]) : pubchem.properties([cid]),
      pubchem.synonyms(cid).catch(() => []),
      pubchem.nmrSpectra(cid).catch((e) => {
        warnings.push(`PubChem (espectros): ${e.message}`);
        return { spectra: [], links: [] };
      }),
    ]);
    return { props: pp[0] || null, synonyms: syn, nmr };
  })().catch((e) => {
    warnings.push(`PubChem: ${e.message}`);
    return { props: null, synonyms: [], nmr: { spectra: [], links: [] } };
  });

  const predictionTask = runPrediction(mol).catch((e) => {
    warnings.push(`Predicción no disponible: ${e.message}`);
    return null;
  });

  const [pc, prediction] = await Promise.all([pubchemTask, predictionTask]);
  const inchikey = pc.props?.inchikey || input.inchikey || null;

  let dbEntries = [];
  if (nmrshiftdb.isReady()) {
    dbEntries = nmrshiftdb.findByStructure({ idcode: info.idcode, idcodeNoStereo: info.idcodeNoStereo, inchikey });
  } else {
    warnings.push(`NMRShiftDB2 local: ${nmrshiftdb.status.state}${nmrshiftdb.status.detail ? ` (${nmrshiftdb.status.detail})` : ''}`);
  }

  const name = input.name || pc.props?.title || null;
  const compound = {
    idcode: info.idcode,
    name,
    iupac: pc.props?.iupac || null,
    formula: info.formula,
    mw: info.mw,
    exactMass: info.exactMass,
    smiles: info.smiles,
    inchikey,
    cid,
    molfile: info.molfile,
    info: {
      unsaturation: info.unsaturation,
      nC: info.nC,
      nH: info.nH,
      elements: info.elements,
      synonyms: pc.synonyms,
      pubchemLinks: pc.nmr.links,
      nmrshiftdbIds: dbEntries.map((e) => e.id),
    },
  };
  const compoundId = upsertCompound(compound);

  const spectra = [];
  const c13freq = carbonFrequency(frequency);
  let protonPred = null;
  let carbonPred = null;
  if (prediction) {
    try {
      protonPred = protonAnalysis(mol, prediction, frequency);
      carbonPred = carbonAnalysis(mol, prediction);
    } catch (e) {
      warnings.push(`Error interpretando la predicción: ${e.message}`);
    }
  }
  if (protonPred) {
    const hTypes = protonPred.signals.length;
    spectra.push({
      nucleus: '1H',
      kind: 'predicho',
      source: 'Predicción (nmr-processing / NMRium)',
      sourceId: '',
      frequency,
      solvent,
      signals: protonPred.signals,
      peaks: protonPred.lines,
      analysis: {
        literature: literatureString('1H', frequency, solvent, protonPred.signals),
        environments: hTypes,
        totalH: protonPred.signals.reduce((a, s) => a + s.nH, 0),
        integrationRatio: protonPred.signals.map((s) => s.relativeIntegration).join(' : '),
      },
      meta: {
        rawSignals: protonPred.rawSignals,
        range: protonPred.range,
        integrals: protonPred.integrals,
        lineWidth: 1,
        method:
          'Desplazamientos por códigos HOSE/red neuronal (servicio de predicción de NMRium); espectro simulado como sistema de espín completo (incluye efectos de segundo orden).',
      },
    });
  }
  if (carbonPred) {
    spectra.push({
      nucleus: '13C',
      kind: 'predicho',
      source: 'Predicción (nmr-processing / NMRium)',
      sourceId: '',
      frequency: c13freq,
      solvent,
      signals: carbonPred.signals,
      peaks: carbonPred.signals.map((s) => ({ ppm: s.delta, intensity: s.intensity, type: s.type, atoms: s.atoms })),
      analysis: {
        literature: literatureString('13C', c13freq, solvent, carbonPred.signals),
        environments: carbonPred.signals.length,
        totalC: carbonPred.signals.reduce((a, s) => a + s.nC, 0),
        dept: {
          CH3: carbonPred.signals.filter((s) => s.nH === 3).length,
          CH2: carbonPred.signals.filter((s) => s.nH === 2).length,
          CH: carbonPred.signals.filter((s) => s.nH === 1).length,
          C: carbonPred.signals.filter((s) => s.nH === 0).length,
        },
      },
      meta: {
        method: 'Espectro desacoplado de protones {¹H}; intensidades aproximadas (proporcionales al nº de C equivalentes, C cuaternarios atenuados).',
      },
    });
  }
  const predBy = { '1H': protonPred, '13C': carbonPred };
  const experimental = [...pubchemSpectra(pc.nmr), ...nmrshiftdbSpectra(dbEntries, mol)];
  for (const e of experimental) {
    const mapping = e.meta?.mapping || null;
    if (e.signals.length) {
      e.analysis = {
        comparison: compareToPrediction(e.nucleus, e.signals.map((s) => ({ delta: s.delta })), predBy[e.nucleus], null),
      };
      if (mapping && e.meta?.rawSignals) {
        const byAtom = compareToPrediction(e.nucleus, e.meta.rawSignals, predBy[e.nucleus], mapping);
        if (byAtom) e.analysis.comparison = { ...byAtom, method: 'por átomo asignado' };
      }
      if (e.analysis.comparison && !e.analysis.comparison.method) e.analysis.comparison.method = 'por cercanía de desplazamientos';
    }
    if (e.meta) {
      delete e.meta.mapping;
    }
    spectra.push(e);
  }

  for (const s of spectra) s.id = upsertSpectrum(compoundId, s);

  const fg = protonPred && carbonPred ? functionalGroups(carbonPred.signals, protonPred.signals) : [];
  const summary = {
    formula: info.formula,
    unsaturation: info.unsaturation,
    protonEnvironments: protonPred?.signals.length ?? null,
    carbonEnvironments: carbonPred?.signals.length ?? null,
    nC: info.nC,
    nH: info.nH,
    symmetry: carbonPred && carbonPred.signals.length < info.nC ? `Hay simetría: ${info.nC} carbonos dan ${carbonPred.signals.length} señales de ¹³C.` : null,
    functionalGroups: fg,
    experimentalCount: experimental.length,
    notes: [
      `Índice de deficiencia de hidrógeno (insaturaciones): ${info.unsaturation} (anillos + enlaces π).`,
      info.unsaturation >= 4 ? 'IDH ≥ 4: compatible con un anillo aromático.' : null,
    ].filter(Boolean),
  };

  return {
    compound: { id: compoundId, ...compound, svg: info.svg, pubchemUrl: cid ? pubchem.pubchemUrl(cid) : null },
    summary,
    spectra: spectra.map(publicSpectrum),
    warnings,
  };
}

export function publicSpectrum(s) {
  const meta = { ...(s.meta || {}) };
  delete meta.rawSignals;
  return {
    id: s.id,
    nucleus: s.nucleus,
    kind: s.kind,
    source: s.source,
    sourceId: s.sourceId,
    frequency: s.frequency,
    solvent: s.solvent,
    signals: s.signals,
    peaks: s.peaks,
    analysis: s.analysis,
    meta,
    image: `/api/spectra/${s.id}/image.svg`,
  };
}

export function savedAnalysis(compound) {
  const mol = parseMolecule({ molfile: compound.molfile, smiles: compound.smiles });
  return {
    compound: {
      ...compound,
      svg: moleculeSvg(mol),
      pubchemUrl: compound.cid ? pubchem.pubchemUrl(compound.cid) : null,
    },
    spectra: getCompoundSpectra(compound.id).map(publicSpectrum),
  };
}
