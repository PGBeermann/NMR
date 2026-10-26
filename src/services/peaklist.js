// Lectura de listas de picos (CSV) y análisis del espectro experimental:
// agrupación en multipletes, centro, constantes J aparentes e integración relativa.
import { carbonRegion, protonRegion, SOLVENTS, multiplicityName } from './chemistry.js';
import { signalLetter } from './molecule.js';

const r = (v, d = 3) => Math.round(v * 10 ** d) / 10 ** d;

/**
 * Acepta separadores coma, punto y coma, tabulador o espacios; encabezado opcional;
 * coma decimal si el separador es ";" o tabulador. Columnas: ppm [, intensidad].
 */
export function parsePeakCsv(text) {
  const lines = String(text || '').replace(/\r/g, '').split('\n').map((l) => l.trim()).filter(Boolean);
  if (!lines.length) throw Object.assign(new Error('La tabla está vacía.'), { status: 400 });
  const sample = lines.slice(0, 10).join('\n');
  // Tabulador o ";" → la coma, si aparece, es decimal. "7.262,581" → coma separadora.
  const sep = /\t/.test(sample) ? '\t' : /;/.test(sample) ? ';' : /-?\d+(\.\d+)?\s*,\s*-?\d/.test(sample) && /\d\.\d/.test(sample) ? ',' : null;
  const peaks = [];
  let ppmCol = 0;
  let intCol = 1;
  let hzCol = -1;
  let areaCol = -1;
  let widthCol = -1;
  for (const line of lines) {
    const cells = (sep ? line.split(sep) : line.split(/\s+/)).map((c) => c.trim().replace(/^"|"$/g, ''));
    const nums = cells.map((c) => {
      const v = sep !== ',' ? c.replace(',', '.') : c;
      return /^[-+]?\d*\.?\d+(e[-+]?\d+)?$/i.test(v) ? Number(v) : NaN;
    });
    if (nums.every((n) => Number.isNaN(n))) {
      // Encabezado: anchura y área se reconocen antes que la altura ("Integral" es área, no altura)
      const kind = (c) => {
        if (/w1\/2|w½|width|ancho|anchura|fwhm|^lw\b/.test(c)) return 'width';
        if (/área|area|integral|integr/.test(c)) return 'area';
        if (/height|altura|amp|intens|^int\.?$|^int\b|^i$/.test(c)) return 'height';
        if (/ppm|shift|δ|delta|despl/.test(c)) return 'ppm';
        if (/^hz$|^\(hz\)$|frec|freq|^ν/.test(c)) return 'hz';
        return null;
      };
      cells.map((c) => c.toLowerCase()).forEach((c, i) => {
        const k = kind(c);
        if (k === 'ppm') ppmCol = i;
        if (k === 'height') intCol = i;
        if (k === 'hz') hzCol = i;
        if (k === 'area' && areaCol < 0) areaCol = i;
        if (k === 'width' && widthCol < 0) widthCol = i;
      });
      // Si sólo hay columna de área (sin altura), se usa también como intensidad
      if (areaCol >= 0 && intCol === areaCol) intCol = -1;
      continue;
    }
    const ppm = nums[ppmCol];
    if (!Number.isFinite(ppm)) continue;
    const area = areaCol >= 0 && Number.isFinite(nums[areaCol]) ? Math.abs(nums[areaCol]) : null;
    const width = widthCol >= 0 && Number.isFinite(nums[widthCol]) ? Math.abs(nums[widthCol]) : null;
    let height = intCol >= 0 && Number.isFinite(nums[intCol]) ? Math.abs(nums[intCol]) : null;
    if (height === null) height = area ?? 1;
    peaks.push({ ppm, height, area, width, hz: hzCol >= 0 ? nums[hzCol] : null });
  }
  if (!peaks.length) throw Object.assign(new Error('No se encontraron picos válidos (columnas: ppm, intensidad).'), { status: 400 });
  peaks.sort((a, b) => b.ppm - a.ppm);
  const max = Math.max(...peaks.map((p) => p.height)) || 1;
  const hasArea = peaks.every((p) => p.area !== null);
  const hasWidth = !hasArea && peaks.every((p) => p.width !== null && p.width > 0);
  const out = peaks.map((p) => ({
    ppm: p.ppm,
    intensity: r((p.height / max) * 100, 2),
    // Área bruta de la línea si el archivo la trae (área o altura × anchura)
    area: hasArea ? p.area : hasWidth ? p.height * p.width : null,
    width: p.width,
  }));
  // Si la tabla trae una columna en Hz, la frecuencia del espectrómetro es Hz/ppm
  const ratios = peaks.filter((p) => Number.isFinite(p.hz) && Math.abs(p.ppm) > 0.5).map((p) => p.hz / p.ppm).sort((a, b) => a - b);
  out.frequency = ratios.length ? Math.round(ratios[Math.floor(ratios.length / 2)]) : null;
  out.areaSource = hasArea ? 'columna de área del archivo' : hasWidth ? 'altura × anchura (W½) del archivo' : null;
  return out;
}

/**
 * Áreas de las líneas cuando el archivo sólo trae alturas: cada línea se modela como una
 * lorentziana de anchura común y se resuelve el solapamiento entre líneas cercanas
 * (deconvolución no negativa), de modo que el área de cada línea es proporcional a su
 * amplitud real y no a la altura aparente.
 */
export function estimateAreas(lines, freq, widthHz) {
  const n = lines.length;
  if (n === 1) return [lines[0].intensity * widthHz];
  const x = lines.map((l) => l.ppm * freq);
  const h = lines.map((l) => l.intensity);
  const L = (d) => 1 / (1 + ((2 * d) / widthHz) ** 2);
  const a = h.slice();
  for (let it = 0; it < 300; it++) {
    let change = 0;
    for (let i = 0; i < n; i++) {
      let s = h[i];
      for (let j = 0; j < n; j++) if (j !== i) s -= a[j] * L(x[i] - x[j]);
      const v = Math.max(0, s);
      change += Math.abs(v - a[i]);
      a[i] = v;
    }
    if (change < 1e-6) break;
  }
  // Área de una lorentziana = π/2 · amplitud · anchura (la constante se cancela al normalizar)
  return a.map((v) => v * widthHz);
}

/** Anchura de línea común estimada (Hz): algo menor que la menor separación entre líneas resueltas. */
export function estimateLineWidth(peaks, freq) {
  const hz = peaks.map((p) => p.ppm * freq).sort((a, b) => a - b);
  let min = Infinity;
  for (let i = 1; i < hz.length; i++) {
    const d = hz[i] - hz[i - 1];
    if (d > 0.05 && d < min) min = d;
  }
  if (!Number.isFinite(min)) return 1;
  return Math.min(1.5, Math.max(0.3, 0.6 * min));
}

export function detectNucleus(peaks) {
  const max = Math.max(...peaks.map((p) => p.ppm));
  return max > 20 ? '13C' : '1H';
}

function clusterValues(values, tol) {
  const sorted = [...values].sort((a, b) => a - b);
  const groups = [];
  for (const v of sorted) {
    const g = groups[groups.length - 1];
    if (g && v - g[g.length - 1] <= tol) g.push(v);
    else groups.push([v]);
  }
  return groups.map((g) => ({ value: g.reduce((a, b) => a + b, 0) / g.length, n: g.length }));
}

/**
 * Fusiona líneas separadas menos de `resHz` (estructura fina, acoplamientos a larga
 * distancia o efectos de segundo orden) en "sub-picos" con su centro ponderado.
 */
function mergeFine(lines, freq, resHz) {
  const out = [];
  for (const l of lines) {
    const g = out[out.length - 1];
    if (g && (g.last - l.ppm) * freq < resHz) {
      g.items.push(l);
      g.last = l.ppm;
    } else out.push({ items: [l], last: l.ppm });
  }
  return out.map((g) => {
    const s = g.items.reduce((a, x) => a + x.intensity, 0) || 1;
    return {
      ppm: g.items.reduce((a, x) => a + x.ppm * x.intensity, 0) / s,
      intensity: s,
      max: Math.max(...g.items.map((x) => x.intensity)),
      n: g.items.length,
    };
  });
}

/** Multiplicidad de primer orden a partir de las posiciones de las líneas (Hz). */
function firstOrder(rawLines, freq, resHz = 2.5) {
  const lines = mergeFine(rawLines, freq, resHz);
  const fine = lines.length < rawLines.length;
  const res = firstOrderPattern(lines, freq);
  return { ...res, fine };
}

function firstOrderPattern(lines, freq) {
  const hz = lines.map((l) => l.ppm * freq);
  const n = hz.length;
  if (n === 1) return { multiplicity: 's', J: [] };
  const sp = [];
  for (let i = 0; i < n - 1; i++) sp.push(hz[i] - hz[i + 1]);
  const mean = sp.reduce((a, b) => a + b, 0) / sp.length;
  const tol = Math.max(0.8, 0.18 * mean);
  const equal = sp.every((s) => Math.abs(s - mean) < tol);
  if (n === 2) return { multiplicity: 'd', J: [r(mean, 1)] };
  if (equal && n <= 7) {
    const names = { 3: 't', 4: 'q', 5: 'quint', 6: 'sext', 7: 'sept' };
    return { multiplicity: names[n], J: [r(mean, 1)] };
  }
  if (n === 4 && Math.abs(sp[0] - sp[2]) < tol) return { multiplicity: 'dd', J: [r(sp[0] + sp[1], 1), r(sp[0], 1)] };
  if (n === 6 && Math.abs(sp[0] - sp[2]) < tol && Math.abs(sp[2] - sp[4]) < tol && Math.abs(sp[1] - sp[3]) < tol) {
    return sp[1] > sp[0]
      ? { multiplicity: 'td', J: [r(sp[0] + sp[1], 1), r(sp[0], 1)] }
      : { multiplicity: 'dt', J: [r(sp[0] * 2 + sp[1], 1), r(sp[0], 1)] };
  }
  if (n === 8 && Math.abs(sp[0] - sp[2]) < tol && Math.abs(sp[0] - sp[4]) < tol && Math.abs(sp[0] - sp[6]) < tol) {
    return { multiplicity: 'ddd', J: [r(sp[0] + sp[1] + sp[2] + sp[3], 1), r(sp[0] + sp[1], 1), r(sp[0], 1)] };
  }
  return { multiplicity: 'm', J: [] };
}

/** Constantes J aparentes: separaciones entre líneas que se repiten. */
function apparentJ(lines, freq) {
  const hz = lines.map((l) => l.ppm * freq);
  const diffs = [];
  for (let i = 0; i < hz.length; i++) for (let k = i + 1; k < hz.length; k++) diffs.push(hz[i] - hz[k]);
  const groups = clusterValues(diffs.filter((d) => d > 0.4 && d < 20), 0.35)
    .filter((g) => g.n >= 2)
    .sort((a, b) => b.n - a.n)
    .slice(0, 4)
    .map((g) => r(g.value, 1))
    .sort((a, b) => b - a);
  return groups;
}

/**
 * Agrupa las líneas de un espectro de 1H en multipletes.
 * Dos líneas pertenecen al mismo multiplete si las separa menos de `gapHz`.
 */
export function analyzePeaks(
  peaks,
  { nucleus, frequency = 400, solvent = 'CDCl3', excludeSolvent = false, gapHz = 15, resolutionHz = 2.5, totalH = null, integrals = null } = {},
) {
  const solv = SOLVENTS[solvent] || SOLVENTS.CDCl3;
  const solventShifts = nucleus === '1H' ? [...solv.H, ...(solv.water ? [solv.water] : [])] : solv.C;
  const notes = [];

  // Origen de las áreas: columna del archivo, altura × anchura, o estimación lorentziana
  const fileAreas = peaks.length && peaks.every((p) => p.area !== null && p.area !== undefined);
  const widthHz = fileAreas ? null : estimateLineWidth(peaks, frequency);
  const areaSource = fileAreas
    ? peaks.areaSource || 'columna de área del archivo'
    : `estimada a partir de las alturas (líneas lorentzianas de anchura común ≈ ${widthHz.toFixed(2)} Hz, con deconvolución de las líneas solapadas)`;

  if (nucleus === '13C') {
    const areas = fileAreas ? peaks.map((p) => p.area) : estimateAreas(peaks, frequency, widthHz);
    const maxA = Math.max(...areas) || 1;
    const signals = peaks
      .map((p, i) => ({
        delta: r(p.ppm, 2),
        height: p.intensity,
        intensity: r((areas[i] / maxA) * 100, 2),
        region: carbonRegion(p.ppm),
        solvent: solventShifts.some((s) => Math.abs(s - p.ppm) < 0.3),
      }))
      .filter((s) => !(excludeSolvent && s.solvent));
    // Letra de cada pico de ¹³C (de mayor a menor δ)
    signals.forEach((s, i) => (s.label = signalLetter(i)));
    if (signals.some((s) => s.solvent)) notes.push(`Hay señales que coinciden con el disolvente (${solvent}).`);
    notes.push(`Intensidades por área (${areaSource}).`);
    return { nucleus, frequency, signals, notes, areaSource };
  }

  const gap = gapHz / frequency;
  const groups = [];
  for (const p of peaks) {
    const g = groups[groups.length - 1];
    if (g && g[g.length - 1].ppm - p.ppm <= gap) g.push(p);
    else groups.push([p]);
  }
  let signals = groups.map((lines) => {
    // Área de cada línea: la del archivo o la estimada (lorentzianas con deconvolución dentro del multiplete)
    const lineAreas = fileAreas ? lines.map((l) => l.area) : estimateAreas(lines, frequency, widthHz);
    const sum = lineAreas.reduce((a, v) => a + v, 0) || 1;
    const center = lines.reduce((a, l, i) => a + l.ppm * lineAreas[i], 0) / sum;
    const fo = firstOrder(lines, frequency, resolutionHz);
    const width = (lines[0].ppm - lines[lines.length - 1].ppm) * frequency;
    return {
      delta: r(center, 3),
      from: r(lines[lines.length - 1].ppm, 4),
      to: r(lines[0].ppm, 4),
      widthHz: r(width, 1),
      nLines: lines.length,
      multiplicity: fo.multiplicity,
      multiplicityName: multiplicityName(fo.multiplicity) + (fo.fine && fo.multiplicity !== 'm' ? ' (aparente, con estructura fina)' : ''),
      J: fo.J.length ? fo.J : [],
      apparentJ: fo.multiplicity === 'm' ? apparentJ(lines, frequency) : fo.J,
      area: sum,
      region: protonRegion(center),
      solvent: lines.length <= 2 && solventShifts.some((s) => Math.abs(s - center) < 0.01),
      lines: lines.map((l, i) => ({ ppm: l.ppm, hz: r(l.ppm * frequency, 2), intensity: l.intensity, area: lineAreas[i] })),
    };
  });
  const solventSignals = signals.filter((s) => s.solvent);
  if (solventSignals.length) {
    notes.push(
      `La señal en δ ${solventSignals.map((s) => s.delta.toFixed(2)).join(', ')} coincide con el disolvente/agua residual (${solvent})` +
        (excludeSolvent ? ' y se excluyó del análisis.' : '; puede excluirse marcando la opción correspondiente.'),
    );
  }
  // Disolvente superpuesto a un multiplete (p. ej. CHCl3 7.26 dentro de un multiplete aromático)
  for (const s of signals) {
    if (!s.solvent && solventShifts.some((x) => x >= s.from - 0.005 && x <= s.to + 0.005)) {
      notes.push(`El multiplete en δ ${s.delta.toFixed(2)} contiene la posición del disolvente residual: su integral puede estar sobreestimada.`);
    }
  }
  if (excludeSolvent) signals = signals.filter((s) => !s.solvent);

  // Integración relativa: suma de las áreas de las líneas de cada multiplete
  const minArea = Math.min(...signals.map((s) => s.area));
  const totalArea = signals.reduce((a, s) => a + s.area, 0);
  for (const s of signals) {
    s.relativeIntegral = r(s.area / minArea, 2);
    s.fraction = r(s.area / totalArea, 4);
  }
  // Con la fórmula: reparto de los H totales (redondeo de mayor resto, mínimo 1 H por señal)
  if (totalH && totalH >= signals.length) {
    const raw = signals.map((s) => s.fraction * totalH);
    const ints = raw.map((v) => Math.max(1, Math.floor(v)));
    let left = totalH - ints.reduce((a, b) => a + b, 0);
    const order = raw.map((v, i) => [v - Math.floor(v), i]).sort((a, b) => b[0] - a[0]);
    for (let k = 0; left > 0 && k < order.length * 3; k++, left--) ints[order[k % order.length][1]]++;
    signals.forEach((s, i) => (s.estimatedH = ints[i]));
  }
  // Relación entera aproximada (probando múltiplos 1–4 del menor)
  // Se elige el múltiplo más pequeño cuyo error de redondeo medio sea aceptable (≤ 0.2)
  let best = null;
  for (let k = 1; k <= 4; k++) {
    const ints = signals.map((s) => s.relativeIntegral * k);
    const err = ints.reduce((a, v) => a + Math.abs(v - Math.round(v)), 0) / ints.length;
    const cand = { k, err, ints: ints.map((v) => Math.max(1, Math.round(v))) };
    if (err <= 0.2) {
      best = cand;
      break;
    }
    if (!best || err < best.err) best = cand;
  }
  signals.forEach((s, i) => {
    s.integerIntegral = best.ints[i];
    // Letra que identifica el multiplete en el espectro, en la tabla y en las estructuras candidatas
    s.label = signalLetter(i);
    s.nH = s.estimatedH ?? s.integerIntegral;
  });
  // Integrales fijadas por el usuario (nº de H de cada multiplete)
  let userIntegrals = false;
  if (Array.isArray(integrals) && integrals.length === signals.length && integrals.every((v) => Number(v) > 0)) {
    const sumU = integrals.reduce((a, v) => a + Number(v), 0);
    signals.forEach((s, i) => {
      s.nH = Number(integrals[i]);
      s.fraction = r(s.nH / sumU, 4);
    });
    userIntegrals = true;
  } else {
    notes.push(
      `Integrales por área: ${areaSource}. Relación aproximada ${best.ints.join(' : ')}${
        totalH ? `; repartiendo ${totalH} H: ${signals.map((s) => s.estimatedH).join(' : ')}` : ''
      }.${
        fileAreas
          ? ''
          : ' Sin columna de área o de anchura (W½) en el archivo, las señales anchas (OH, NH) pueden quedar subestimadas: exporte la tabla con el área de cada pico o corrija el nº de H y recalcule.'
      }`,
    );
  }
  for (const s of signals) delete s.area;
  return { nucleus, frequency, signals, notes, ratio: best.ints.join(' : '), userIntegrals, areaSource };
}

/**
 * Interpretación estructural a partir de los multipletes (reglas empíricas):
 * fragmentos probables con la evidencia que los sustenta.
 */
export function interpret(analysis) {
  const out = [];
  const S = analysis.signals;
  const add = (fragment, evidence) => out.push({ fragment, evidence });
  if (analysis.nucleus === '13C') {
    const n = (a, b) => S.filter((s) => !s.solvent && s.delta >= a && s.delta < b).length;
    if (n(190, 230)) add('Cetona o aldehído', `${n(190, 230)} señal(es) entre 190 y 220 ppm`);
    if (n(160, 190)) add('Ácido carboxílico, éster, amida o anhídrido', `${n(160, 190)} señal(es) entre 160 y 185 ppm`);
    if (n(100, 160)) add('Carbonos aromáticos o alquénicos', `${n(100, 160)} señal(es) entre 100 y 160 ppm${n(100, 160) === 4 ? ' (4 señales: típico de un benceno mono- o para-disustituido)' : ''}`);
    if (n(110, 125) && n(100, 160) === 1) add('Posible nitrilo', 'señal aislada entre 110 y 125 ppm');
    if (n(65, 90)) add('C–O (alcohol, éter, éster) o alquino', `${n(65, 90)} señal(es) entre 65 y 90 ppm`);
    if (n(50, 65)) add('C–O / C–N (p. ej. OCH₃ ~55, NCH)', `${n(50, 65)} señal(es) entre 50 y 65 ppm`);
    if (n(-5, 50)) add('Carbonos alifáticos sp³', `${n(-5, 50)} señal(es) por debajo de 50 ppm`);
    return out;
  }
  const nh = (s) => s.nH;
  const near = (s, a, b) => s.delta >= a && s.delta < b;
  const arom = S.filter((s) => near(s, 6.4, 8.6));
  const arH = arom.reduce((a, s) => a + nh(s), 0);
  if (arH) {
    const dd2 = arom.filter((s) => s.multiplicity === 'd' && Math.abs(nh(s) - 2) < 0.6);
    let pattern = '';
    if (Math.round(arH) === 5) pattern = 'Benceno monosustituido (C₆H₅–)';
    else if (Math.round(arH) === 4 && dd2.length === 2) pattern = 'Benceno para-disustituido (sistema AA′BB′: dos dobletes de 2H)';
    else if (Math.round(arH) === 4) pattern = 'Benceno disustituido (orto o meta)';
    else if (Math.round(arH) === 3) pattern = 'Benceno trisustituido (o heteroaromático)';
    else pattern = 'Anillo(s) aromático(s)';
    add(pattern, `${r(arH, 1)} H aromáticos entre 6.4 y 8.6 ppm (${arom.map((s) => `${s.delta.toFixed(2)} ${s.multiplicity}, ${nh(s)}H`).join('; ')})`);
    const minAr = Math.min(...arom.map((s) => s.delta));
    const maxAr = Math.max(...arom.map((s) => s.delta));
    if (minAr < 7.0) add('Sustituyente donador de electrones en el anillo (–OH, –OR, –NH₂/–NR₂)', `H aromáticos apantallados (δ ${minAr.toFixed(2)} < 7.0 ppm, típicos de posiciones orto/para a un donador)`);
    if (maxAr > 7.9) add('Sustituyente aceptor de electrones (C=O, NO₂, CN, COOR)', `H aromáticos desapantallados (δ ${maxAr.toFixed(2)} > 7.9 ppm)`);
    if (Math.round(arH) === 5 && arom.length === 3 && arom.some((s) => s.multiplicity === 't')) {
      add('Patrón orto (d, 2H) / meta (t, 2H) / para (t, 1H) resuelto', 'tres multipletes 2H : 2H : 1H típicos de C₆H₅–X con X que modifica la densidad electrónica');
    }
  }
  for (const s of S) {
    const n = nh(s);
    if (s.multiplicity === 's' && Math.abs(n - 3) < 0.6 && near(s, 3.55, 4.05)) add('Metoxilo (–OCH₃)', `singlete de 3H en δ ${s.delta.toFixed(2)} (éster metílico ~3.7; anisol ~3.8)`);
    if (s.multiplicity === 's' && Math.abs(n - 3) < 0.6 && near(s, 1.95, 2.6)) {
      // Ar–CH₃ (2.2–2.5) y acetilo (2.0–2.2) se solapan: se decide por la presencia de H aromáticos
      if (arH && s.delta >= 2.15) add('Metilo sobre anillo aromático (Ar–CH₃)', `singlete de 3H en δ ${s.delta.toFixed(2)} y presencia de H aromáticos (un acetilo daría ~2.0–2.2 ppm)`);
      else if (s.delta < 2.2) add('Acetilo (CH₃–C=O)', `singlete de 3H en δ ${s.delta.toFixed(2)}`);
      else add('Metilcetona (CH₃–C=O) o CH₃ alílico', `singlete de 3H en δ ${s.delta.toFixed(2)}`);
    }
    if (s.multiplicity === 's' && Math.abs(n - 9) < 1 && near(s, 0.8, 1.6)) add('terc-Butilo (–C(CH₃)₃)', `singlete de 9H en δ ${s.delta.toFixed(2)}`);
    if (near(s, 9.3, 10.6) && n <= 1.5) add('Aldehído (–CHO)', `señal de 1H en δ ${s.delta.toFixed(2)}`);
    if (s.delta >= 10.5 && n <= 1.5) add('Ácido carboxílico (–COOH) o fenol con puente de H', `señal de 1H en δ ${s.delta.toFixed(2)}`);
    if (s.multiplicity === 's' && n <= 1.5 && near(s, 4.0, 7.0) && !near(s, 6.4, 7.0))
      add('Posible OH (fenol/alcohol) o NH intercambiable', `singlete de 1H en δ ${s.delta.toFixed(2)} sin acoplamientos (confirmar con D₂O: desaparece si es OH/NH)`);
  }
  // Etilo: triplete 3H + cuarteto 2H con la misma J
  const t3 = S.filter((s) => s.multiplicity === 't' && Math.abs(nh(s) - 3) < 0.6);
  const q2 = S.filter((s) => s.multiplicity === 'q' && Math.abs(nh(s) - 2) < 0.6);
  for (const t of t3)
    for (const q of q2) {
      if (t.J[0] && q.J[0] && Math.abs(t.J[0] - q.J[0]) < 0.8) {
        const where = q.delta > 3.8 ? 'unido a O (éster etílico / etoxilo, –OCH₂CH₃)' : q.delta > 2.2 ? 'en α a C=O o a un anillo (–CH₂CH₃)' : 'alifático';
        add(`Grupo etilo ${where}`, `t 3H (δ ${t.delta.toFixed(2)}) y q 2H (δ ${q.delta.toFixed(2)}) con la misma J ≈ ${t.J[0]} Hz`);
      }
    }
  // Isopropilo: doblete 6H
  for (const d6 of S.filter((s) => s.multiplicity === 'd' && Math.abs(nh(s) - 6) < 0.8)) add('Grupo isopropilo (–CH(CH₃)₂)', `doblete de 6H en δ ${d6.delta.toFixed(2)}`);
  return out;
}

/** Texto en formato de publicación. */
export function literature(analysis, solvent) {
  const s = analysis.signals;
  if (analysis.nucleus === '13C') return `¹³C RMN (${analysis.frequency} MHz, ${solvent}) δ ${s.map((x) => x.delta.toFixed(1)).join(', ')}.`;
  const items = s.map((x) => {
    const J = x.J.length ? `, J = ${x.J.join(', ')} Hz` : '';
    const range = x.multiplicity === 'm' ? `${x.to.toFixed(2)}–${x.from.toFixed(2)}` : x.delta.toFixed(2);
    return `${range} (${x.multiplicity}${J}, ${x.nH}H)`;
  });
  return `¹H RMN (${analysis.frequency} MHz, ${solvent}) δ ${items.join(', ')}.`;
}
