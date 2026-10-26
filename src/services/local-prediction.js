// Predicción local de respaldo (cuando el servicio NMRium no responde).
// Desplazamientos: entornos tipo HOSE construidos a partir de NMRShiftDB2 (hose.js).
// Acoplamientos: regla n+1 con constantes ³J típicas según el tipo de enlace.
// Devuelve un objeto con la misma forma que la predicción de nmr-processing, para que
// el resto del análisis (tablas, simulación, imágenes, guardado) funcione igual.
import OCL from 'openchemlib';
import { hoseReady, predictLocal } from './hose.js';

const LETTERS = ['s', 'd', 't', 'q', 'quint', 'sext', 'sept'];
const HETERO = [7, 8, 16];

/** ³J típica (Hz) entre los H de dos átomos vecinos a y b. */
function typicalJ(mol, a, b) {
  const bond = mol.getBond(a, b);
  const isCarbonyl = (x) => {
    for (let i = 0; i < mol.getConnAtoms(x); i++) {
      if (mol.getAtomicNo(mol.getConnAtom(x, i)) === 8 && mol.getConnBondOrder(x, i) === 2) return true;
    }
    return false;
  };
  if (mol.isAromaticAtom(a) && mol.isAromaticAtom(b)) return 7.8; // orto aromático
  if (isCarbonyl(a) || isCarbonyl(b)) return 2.5; // H de aldehído con CH vecino
  if (bond >= 0 && mol.getBondOrder(bond) === 2) return 10.0; // vinílico (promedio cis/trans)
  return 7.0; // alifático
}

export function localPredictionAvailable() {
  return hoseReady();
}

/**
 * @param {OCL.Molecule} mol molécula con hidrógenos implícitos (como la de parseMolecule)
 */
export function localPrediction(mol) {
  if (!hoseReady()) return null;
  const m = mol.getCompactCopy();
  m.removeExplicitHydrogens();
  m.ensureHelperArrays(OCL.Molecule.cHelperRings);
  const hGroups = predictLocal(m, '1H') || [];
  const cGroups = predictLocal(m, '13C') || [];
  if (!hGroups.length && !cGroups.length) return null;

  // Molécula con H explícitos (los índices de los átomos pesados no cambian)
  const mh = m.getCompactCopy();
  mh.addImplicitHydrogens();
  mh.ensureHelperArrays(OCL.Molecule.cHelperNeighbours);
  const hOf = (a) => {
    const out = [];
    for (let i = 0; i < mh.getAllConnAtoms(a); i++) {
      const b = mh.getConnAtom(a, i);
      if (mh.getAtomicNo(b) === 1) out.push(b);
    }
    return out;
  };

  const groupOf = new Map();
  hGroups.forEach((g, i) => g.atoms.forEach((a) => groupOf.set(a, i)));
  const dia = (i) => `local-${i}`;

  const joinedSignals = [];
  const signals = [];
  hGroups.forEach((g, gi) => {
    const a0 = g.atoms[0];
    // Vecinos con H (sin OH/NH/SH ni H equivalentes del mismo grupo), agrupados por J
    const buckets = new Map();
    for (let i = 0; i < m.getConnAtoms(a0); i++) {
      const b = m.getConnAtom(a0, i);
      if (HETERO.includes(m.getAtomicNo(b)) || !m.getAllHydrogens(b)) continue;
      const bg = groupOf.get(b);
      if (bg === gi) continue;
      const J = typicalJ(m, a0, b);
      if (!buckets.has(J)) buckets.set(J, { count: 0, groups: new Set() });
      const k = buckets.get(J);
      k.count += m.getAllHydrogens(b);
      if (bg !== undefined) k.groups.add(dia(bg));
    }
    const js = [...buckets.entries()]
      .sort((x, y) => y[0] - x[0])
      .map(([J, k]) => ({ coupling: J, multiplicity: LETTERS[k.count] || 'm', diaIDs: [...k.groups], pathLength: 3 }));
    let multiplicity = js.length ? js.map((j) => j.multiplicity).join('') : 's';
    if (js.length > 3 || js.some((j) => j.multiplicity === 'm')) multiplicity = 'm';

    const hAtoms = g.atoms.flatMap(hOf);
    joinedSignals.push({ delta: g.delta, nbAtoms: hAtoms.length, atoms: hAtoms, diaIDs: [dia(gi)], js, multiplicity });

    // Señales por átomo de H (para simular el espectro como sistema de espín)
    for (const a of g.atoms) {
      const perHjs = [];
      for (let i = 0; i < m.getConnAtoms(a); i++) {
        const b = m.getConnAtom(a, i);
        if (HETERO.includes(m.getAtomicNo(b)) || groupOf.get(b) === gi) continue;
        const J = typicalJ(m, a, b);
        for (const hb of hOf(b)) perHjs.push({ coupling: J, atoms: [hb], multiplicity: 'd' });
      }
      for (const h of hOf(a)) signals.push({ delta: g.delta, atoms: [h], nbAtoms: 1, js: perHjs });
    }
  });

  const carbon = cGroups.map((g) => ({ delta: g.delta, atoms: g.atoms, nbAtoms: g.n, diaIDs: [], statistics: null }));

  return {
    source: 'local',
    moleculeWithHydrogens: mh,
    spectra: {
      proton: { joinedSignals, signals },
      carbon: { joinedSignals: carbon },
    },
  };
}
