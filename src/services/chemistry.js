// Interpretación química: entorno de cada átomo, rangos típicos de desplazamiento,
// nombres de multiplicidades y regiones del espectro.

const HALOGENS = { 9: 'F', 17: 'Cl', 35: 'Br', 53: 'I' };

function neighbors(mol, a) {
  const out = [];
  for (let i = 0; i < mol.getAllConnAtoms(a); i++) {
    out.push({ atom: mol.getConnAtom(a, i), order: mol.getConnBondOrder(a, i), z: mol.getAtomicNo(mol.getConnAtom(a, i)) });
  }
  return out;
}

function isCarbonyl(mol, c) {
  return mol.getAtomicNo(c) === 6 && neighbors(mol, c).some((n) => n.z === 8 && n.order === 2);
}

function carbonylType(mol, c) {
  const nb = neighbors(mol, c);
  const others = nb.filter((n) => !(n.z === 8 && n.order === 2));
  const oSingle = others.filter((n) => n.z === 8 && n.order === 1);
  if (oSingle.length) {
    if (oSingle.length === 2) return { code: 'carbonato', text: 'C=O de carbonato' };
    const o = oSingle[0].atom;
    if (mol.getAllHydrogens(o) > 0) return { code: 'acido', text: 'C=O de ácido carboxílico' };
    if (neighbors(mol, o).some((n) => n.atom !== c && isCarbonyl(mol, n.atom))) return { code: 'anhidrido', text: 'C=O de anhídrido' };
    return { code: 'ester', text: 'C=O de éster' };
  }
  if (others.some((n) => n.z === 7)) return { code: 'amida', text: 'C=O de amida' };
  if (others.some((n) => n.z === 17 || n.z === 35)) return { code: 'haluro_acilo', text: 'C=O de haluro de acilo' };
  if (mol.getAllHydrogens(c) > 0) return { code: 'aldehido', text: 'C=O de aldehído' };
  return { code: 'cetona', text: 'C=O de cetona' };
}

function oxygenType(mol, o, from) {
  const others = neighbors(mol, o).filter((n) => n.atom !== from);
  if (mol.getAllHydrogens(o) > 0) return 'alcohol';
  if (others.some((n) => isCarbonyl(mol, n.atom))) return 'éster';
  if (others.some((n) => mol.isAromaticAtom(n.atom))) return 'éter arílico';
  return 'éter';
}

const H_LABEL = ['C', 'CH', 'CH₂', 'CH₃', 'CH₄'];

/** Describe el entorno químico de un carbono. Devuelve código, texto y rango típico 13C. */
export function describeCarbon(mol, a) {
  const nH = mol.getAllHydrogens(a);
  const label = H_LABEL[nH] || 'C';
  const nb = neighbors(mol, a);
  const hetero = nb.filter((n) => n.z !== 6 && n.z !== 1);

  if (mol.isAromaticAtom(a)) {
    const subs = hetero.filter((n) => !mol.isAromaticAtom(n.atom) || n.z !== 6);
    let t = nH ? 'CH aromático' : 'C aromático cuaternario (ipso)';
    if (nb.some((n) => n.z === 7 && mol.isAromaticAtom(n.atom))) t += ' adyacente a N del heterociclo';
    else if (subs.some((n) => n.z === 8)) t += ' unido a O';
    else if (subs.some((n) => n.z === 7)) t += ' unido a N';
    else if (subs.some((n) => HALOGENS[n.z])) t += ` unido a ${HALOGENS[subs.find((n) => HALOGENS[n.z]).z]}`;
    return { code: 'aromatico', label, nH, text: t, range13C: [110, 160], range1H: [6.5, 8.5] };
  }
  if (isCarbonyl(mol, a)) {
    const ct = carbonylType(mol, a);
    const ranges = {
      cetona: [195, 220], aldehido: [190, 205], acido: [170, 185], ester: [165, 175],
      amida: [160, 175], haluro_acilo: [165, 180], anhidrido: [160, 175], carbonato: [150, 160],
    };
    return { code: ct.code, label, nH, text: ct.text, range13C: ranges[ct.code], range1H: ct.code === 'aldehido' ? [9.3, 10.5] : null };
  }
  if (nb.some((n) => n.z === 7 && n.order === 3)) return { code: 'nitrilo', label, nH, text: 'C de nitrilo (C≡N)', range13C: [110, 125], range1H: null };
  if (nb.some((n) => n.z === 6 && n.order === 3)) return { code: 'alquino', label, nH, text: nH ? 'CH alquínico terminal' : 'C alquínico', range13C: [65, 90], range1H: [1.8, 3.2] };
  if (nb.some((n) => n.z === 7 && n.order === 2)) return { code: 'imina', label, nH, text: 'C de imina (C=N)', range13C: [150, 170], range1H: [7.5, 8.5] };
  if (nb.some((n) => n.z === 6 && n.order === 2)) {
    const conj = nb.some((n) => n.order === 1 && (isCarbonyl(mol, n.atom) || mol.isAromaticAtom(n.atom)));
    return { code: 'alqueno', label, nH, text: `${label} alquénico (vinílico)${conj ? ' conjugado' : ''}`, range13C: [100, 150], range1H: [4.5, 7.0] };
  }
  // Carbono sp3
  if (hetero.length) {
    const h = hetero[0];
    if (h.z === 8) {
      const t = oxygenType(mol, h.atom, a);
      return { code: 'C-O', label, nH, text: `${label} unido a O (${t})`, range13C: [50, 90], range1H: t === 'éster' ? [3.6, 4.6] : [3.3, 4.2] };
    }
    if (h.z === 7) {
      const amide = neighbors(mol, h.atom).some((n) => isCarbonyl(mol, n.atom));
      return { code: 'C-N', label, nH, text: `${label} unido a N (${amide ? 'amida' : 'amina'})`, range13C: [25, 65], range1H: amide ? [2.8, 3.6] : [2.2, 3.0] };
    }
    if (HALOGENS[h.z]) {
      const x = HALOGENS[h.z];
      const r13 = { F: [70, 100], Cl: [30, 70], Br: [20, 65], I: [-10, 40] }[x];
      return { code: 'C-X', label, nH, text: `${label} unido a ${x}`, range13C: r13, range1H: x === 'F' ? [4.2, 4.8] : [2.9, 4.2] };
    }
    if (h.z === 16) return { code: 'C-S', label, nH, text: `${label} unido a S`, range13C: [10, 45], range1H: [2.0, 3.0] };
    if (h.z === 14) return { code: 'C-Si', label, nH, text: `${label} unido a Si`, range13C: [-5, 20], range1H: [-0.2, 0.5] };
  }
  const cNb = nb.filter((n) => n.z === 6);
  if (cNb.some((n) => isCarbonyl(mol, n.atom))) return { code: 'alfa-CO', label, nH, text: `${label} en α a C=O`, range13C: [20, 55], range1H: [2.0, 2.7] };
  if (cNb.some((n) => mol.isAromaticAtom(n.atom))) return { code: 'bencilico', label, nH, text: `${label} bencílico`, range13C: [15, 45], range1H: [2.2, 3.0] };
  if (cNb.some((n) => mol.getAtomPi(n.atom) > 0)) return { code: 'alilico', label, nH, text: `${label} alílico / propargílico`, range13C: [15, 40], range1H: [1.6, 2.4] };
  return { code: 'alquilo', label, nH, text: `${label} alifático (alquílico)`, range13C: [0, 55], range1H: nH === 3 ? [0.7, 1.3] : nH === 2 ? [1.1, 1.8] : [1.3, 2.0] };
}

/** Describe el entorno del átomo pesado que porta un H. */
export function describeHydrogenCarrier(mol, a) {
  const z = mol.getAtomicNo(a);
  if (z === 6) {
    const d = describeCarbon(mol, a);
    let range = d.range1H;
    if (d.code === 'aldehido') range = [9.3, 10.5];
    return { code: `H-${d.code}`, text: `H de ${d.text}`, range, exchangeable: false };
  }
  if (z === 8) {
    const nb = neighbors(mol, a);
    if (nb.some((n) => isCarbonyl(mol, n.atom))) return { code: 'H-acido', text: 'OH de ácido carboxílico (intercambiable)', range: [10, 13], exchangeable: true };
    if (nb.some((n) => mol.isAromaticAtom(n.atom))) return { code: 'H-fenol', text: 'OH fenólico (intercambiable)', range: [4, 8], exchangeable: true };
    return { code: 'H-alcohol', text: 'OH de alcohol (intercambiable)', range: [0.5, 5.5], exchangeable: true };
  }
  if (z === 7) {
    const amide = neighbors(mol, a).some((n) => isCarbonyl(mol, n.atom));
    return amide
      ? { code: 'H-amida', text: 'NH de amida (intercambiable)', range: [5, 9], exchangeable: true }
      : { code: 'H-amina', text: 'NH de amina (intercambiable)', range: [0.5, 5], exchangeable: true };
  }
  if (z === 16) return { code: 'H-tiol', text: 'SH de tiol', range: [1, 4], exchangeable: true };
  return { code: 'H-otro', text: `H unido a ${mol.getAtomLabel(a)}`, range: null, exchangeable: true };
}

const MULT_WORDS = {
  s: ['singlete', 0], d: ['doblete', 1], t: ['triplete', 2], q: ['cuarteto', 3],
  p: ['quintuplete', 4], quint: ['quintuplete', 4], h: ['sextuplete', 5], sext: ['sextuplete', 5],
  hept: ['septuplete', 6], sept: ['septuplete', 6], m: ['multiplete', null], br: ['ancho', null],
};

export function parseMultiplicity(mult) {
  const tokens = String(mult || 's').match(/hept|sept|quint|sext|br|[sdtqphm]/g) || ['s'];
  return tokens.map((t) => ({ token: t, word: MULT_WORDS[t]?.[0] || t, n: MULT_WORDS[t]?.[1] ?? null }));
}

export function multiplicityName(mult) {
  const parts = parseMultiplicity(mult);
  if (parts.length === 1) return parts[0].word;
  // dd → doblete de dobletes; dt → doblete de tripletes; ddd → doblete de doblete de dobletes
  const plural = (w) => (w.endsWith('e') || w.endsWith('o') ? `${w}s` : w);
  return parts
    .map((p, i) => (i === parts.length - 1 ? plural(p.word) : p.word))
    .join(' de ');
}

export function protonRegion(ppm) {
  if (ppm < 0) return 'Zona apantallada (H unidos a Si/metales, anillos aromáticos internos)';
  if (ppm < 1.9) return 'H alifáticos (CH₃, CH₂, CH saturados)';
  if (ppm < 3.0) return 'H alílicos, bencílicos, en α a C=O o a N';
  if (ppm < 4.6) return 'H sobre carbono unido a O, N o halógeno';
  if (ppm < 6.5) return 'H vinílicos (alquenos); un singlete aislado puede ser OH/NH intercambiable';
  if (ppm < 8.6) return 'H aromáticos / heteroaromáticos';
  if (ppm < 10.5) return 'H de aldehído o amida/NH desapantallado';
  return 'H ácidos (COOH, fenoles con puente de H, enoles)';
}

export function carbonRegion(ppm) {
  if (ppm < 50) return 'C sp³ alifáticos (CH₃, CH₂, CH, C)';
  if (ppm < 65) return 'C sp³ unidos a N/O (aminas, metoxilos) o cuaternarios';
  if (ppm < 95) return 'C sp³ unidos a O (alcoholes, éteres, ésteres) / alquinos';
  if (ppm < 160) return 'C sp² de alquenos y aromáticos (y nitrilos 110–125)';
  if (ppm < 185) return 'C=O de ácidos, ésteres, amidas, anhídridos';
  return 'C=O de aldehídos y cetonas';
}

export function deptPhase(nH) {
  if (nH === 0) return { dept135: 'ausente', dept90: 'ausente' };
  if (nH === 2) return { dept135: 'negativa (↓)', dept90: 'ausente' };
  if (nH === 1) return { dept135: 'positiva (↑)', dept90: 'positiva (↑)' };
  return { dept135: 'positiva (↑)', dept90: 'ausente' };
}

/** Disolventes deuterados comunes: señales residuales (Fulmer et al., Organometallics 2010). */
export const SOLVENTS = {
  CDCl3: { name: 'Cloroformo-d (CDCl₃)', H: [7.26], C: [77.16], water: 1.56 },
  'DMSO-d6': { name: 'DMSO-d₆', H: [2.5], C: [39.52], water: 3.33 },
  D2O: { name: 'Agua deuterada (D₂O)', H: [4.79], C: [], water: null },
  CD3OD: { name: 'Metanol-d₄ (CD₃OD)', H: [3.31], C: [49.0], water: 4.87 },
  'Acetona-d6': { name: 'Acetona-d₆', H: [2.05], C: [29.84, 206.26], water: 2.84 },
  C6D6: { name: 'Benceno-d₆', H: [7.16], C: [128.06], water: 0.4 },
  CD3CN: { name: 'Acetonitrilo-d₃', H: [1.94], C: [1.32, 118.26], water: 2.13 },
  CD2Cl2: { name: 'Diclorometano-d₂ (CD₂Cl₂)', H: [5.32], C: [53.84], water: 1.52 },
  ninguno: { name: 'Sin disolvente (solo TMS)', H: [], C: [], water: null },
};
