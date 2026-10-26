// Utilidades de estructura química basadas en OpenChemLib.
import OCL from 'openchemlib';

const { Molecule } = OCL;

export function parseMolecule({ smiles, molfile, idcode } = {}) {
  let mol;
  try {
    if (molfile && /M\s+END/.test(molfile)) mol = Molecule.fromMolfile(molfile);
    else if (smiles && smiles.trim()) mol = Molecule.fromSmiles(smiles.trim());
    else if (idcode) mol = Molecule.fromIDCode(idcode);
  } catch (e) {
    throw Object.assign(new Error(`Estructura no válida: ${e.message}`), { status: 400 });
  }
  if (!mol || mol.getAllAtoms() === 0) {
    throw Object.assign(new Error('La estructura está vacía. Dibuje una molécula.'), { status: 400 });
  }
  // Se trabaja siempre con hidrógenos implícitos para numerar sólo átomos pesados
  mol.removeExplicitHydrogens();
  mol.ensureHelperArrays(Molecule.cHelperRings);
  // Si viene sin coordenadas (todas cero) se generan coordenadas 2D
  let allZero = true;
  for (let i = 0; i < mol.getAllAtoms() && allZero; i++) {
    if (mol.getAtomX(i) !== 0 || mol.getAtomY(i) !== 0) allZero = false;
  }
  if (allZero && mol.getAllAtoms() > 1) mol.inventCoordinates();
  return mol;
}

export function elementCounts(formula) {
  const counts = {};
  for (const m of String(formula).matchAll(/([A-Z][a-z]?)(\d*)/g)) {
    counts[m[1]] = (counts[m[1]] || 0) + (m[2] ? Number(m[2]) : 1);
  }
  return counts;
}

/** Normaliza una fórmula escrita por el usuario al orden de Hill (C, H, resto alfabético). */
export function normalizeFormula(input) {
  const clean = String(input || '').replace(/\s+/g, '').replace(/[₀-₉]/g, (d) => String(d.charCodeAt(0) - 8320));
  if (!/^([A-Z][a-z]?\d*)+$/.test(clean)) {
    throw Object.assign(new Error('Fórmula molecular no válida. Ejemplo: C4H8O2'), { status: 400 });
  }
  const c = elementCounts(clean);
  const order = [];
  if (c.C) {
    order.push('C');
    if (c.H) order.push('H');
  }
  for (const el of Object.keys(c).sort()) if (!order.includes(el)) order.push(el);
  return order.map((el) => el + (c[el] > 1 ? c[el] : '')).join('');
}

/** Índice de deficiencia de hidrógeno (insaturaciones = anillos + enlaces π). */
export function unsaturation(formula) {
  const c = elementCounts(formula);
  const C = (c.C || 0) + (c.Si || 0);
  const H = c.H || 0;
  const N = (c.N || 0) + (c.P || 0);
  const X = (c.F || 0) + (c.Cl || 0) + (c.Br || 0) + (c.I || 0);
  return (2 * C + 2 + N - H - X) / 2;
}

export function moleculeSvg(mol, { width = 420, height = 320, numbers = true, id = 'mol' } = {}) {
  const m = mol.getCompactCopy();
  if (numbers) {
    // Una etiqueta que empieza por "]" se dibuja junto al símbolo del átomo
    for (let i = 0; i < m.getAllAtoms(); i++) m.setAtomCustomLabel(i, `]${i + 1}`);
  }
  return m.toSVG(width, height, id, {
    autoCrop: true,
    autoCropMargin: 25,
    suppressChiralText: true,
    noCarbonLabelWithCustomLabel: true,
  });
}

export const SIGNAL_LETTERS = 'abcdefghijklmnopqrstuvwxyz'.split('');
export const signalLetter = (i) => (i < 26 ? SIGNAL_LETTERS[i] : SIGNAL_LETTERS[Math.floor(i / 26) - 1] + SIGNAL_LETTERS[i % 26]);

/**
 * Estructura con las letras de las señales de ¹H junto al átomo que porta esos H
 * (p. ej. «a» en el CH₃ cuyo pico es la señal a). `labels`: Map átomo → [letras].
 */
export const LETTER_COLORS = { '1H': 'rgb(31,78,156)', '13C': 'rgb(176,90,0)' };

export function moleculeSvgLabeled(mol, labels, { width = 440, height = 330, id = 'molL', color = LETTER_COLORS['1H'] } = {}) {
  const m = mol.getCompactCopy();
  for (const [a, ls] of labels) {
    if (a >= 0 && a < m.getAllAtoms() && ls?.length) m.setAtomCustomLabel(a, `]${[...new Set(ls)].join(',')}`);
  }
  const svg = m.toSVG(width, height, id, { autoCrop: true, autoCropMargin: 28, suppressChiralText: true, noCarbonLabelWithCustomLabel: true });
  // Las letras (minúsculas) se resaltan: negrita, mayor tamaño y color propio
  return svg.replace(/<text([^>]*)>([a-z]{1,2}(?:,[a-z]{1,2})*)<\/text>/g, (all, attrs, txt) => {
    const clean = attrs.replace(/\s(fill|font-size)="[^"]*"/g, '');
    return `<text${clean} font-size="14" font-weight="bold" fill="${color}">${txt}</text>`;
  });
}

export function moleculeInfo(mol) {
  const f = mol.getMolecularFormula();
  const formula = f.formula;
  const counts = elementCounts(formula);
  const noStereo = mol.getCompactCopy();
  noStereo.stripStereoInformation();
  return {
    formula,
    mw: Number(f.relativeWeight.toFixed(4)),
    exactMass: Number(f.absoluteWeight.toFixed(5)),
    smiles: mol.toIsomericSmiles(),
    idcode: mol.getIDCode(),
    idcodeNoStereo: noStereo.getIDCode(),
    molfile: mol.toMolfile(),
    nC: counts.C || 0,
    nH: counts.H || 0,
    elements: counts,
    unsaturation: unsaturation(formula),
    heavyAtoms: mol.getAtoms(),
    svg: moleculeSvg(mol),
  };
}

/** Correspondencia de átomos entre dos estructuras idénticas (la de la BD y la del usuario). */
export function atomMapping(fromMol, toMol) {
  try {
    const frag = fromMol.getCompactCopy();
    frag.setFragment(true);
    const s = new OCL.SSSearcher();
    s.setMol(frag, toMol);
    if (!s.findFragmentInMolecule()) return null;
    const match = s.getMatchList()[0];
    return match ? Array.from(match) : null;
  } catch {
    return null;
  }
}

export { OCL };
