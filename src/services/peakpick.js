// Detección de picos y áreas sobre un espectro continuo (p. ej. leído de un JCAMP-DX):
// corrección de línea base, nivel de ruido, máximos locales con prominencia mínima,
// posición por interpolación parabólica y área por integración de la curva.
import { SOLVENTS } from './chemistry.js';

const median = (arr) => {
  if (!arr.length) return 0;
  const a = Float64Array.from(arr).sort();
  const m = a.length >> 1;
  return a.length % 2 ? a[m] : (a[m - 1] + a[m]) / 2;
};

/** Línea base por tramos: mediana de cada tramo y mínimo entre tramos vecinos (ignora señales anchas). */
function baseline(y) {
  const n = y.length;
  const seg = Math.max(64, Math.floor(n / 256));
  const nSeg = Math.ceil(n / seg);
  const med = [];
  for (let k = 0; k < nSeg; k++) {
    // Submuestreo para que la mediana sea rápida en tramos grandes
    const a = k * seg;
    const b = Math.min(n, a + seg);
    const step = Math.max(1, Math.floor((b - a) / 256));
    const vals = [];
    for (let i = a; i < b; i += step) vals.push(y[i]);
    med.push(median(vals));
  }
  const anchor = med.map((_, k) => Math.min(...med.slice(Math.max(0, k - 2), k + 3)));
  const out = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const t = i / seg - 0.5;
    const k = Math.max(0, Math.min(nSeg - 2, Math.floor(t)));
    const f = nSeg > 1 ? Math.max(0, Math.min(1, t - k)) : 0;
    out[i] = nSeg > 1 ? anchor[k] * (1 - f) + anchor[k + 1] * f : anchor[0];
  }
  return out;
}

/** Desviación estándar robusta del ruido (MAD de la curva y de sus diferencias; se toma la mayor). */
function noiseLevel(y) {
  const n = y.length;
  const step = Math.max(1, Math.floor(n / 50000));
  const vals = [];
  const diffs = [];
  for (let i = 0; i + step < n; i += step) {
    vals.push(y[i]);
    diffs.push(y[i + 1] - y[i]);
  }
  const m = median(vals);
  const mad = 1.4826 * median(vals.map((v) => Math.abs(v - m)));
  const md = median(diffs);
  const madD = (1.4826 * median(diffs.map((v) => Math.abs(v - md)))) / Math.SQRT2;
  return Math.max(mad, madD) || 1e-12;
}

/**
 * Picos de un espectro de 1D. `x` en ppm (ascendente), `y` intensidades.
 * Devuelve la lista de picos (ppm descendente) con altura, área (integral de la curva entre
 * los mínimos que separan las líneas) y anchura a media altura (Hz), y la curva corregida.
 */
export function pickPeaks(x, y, { nucleus = '1H', frequency = 400, solvent = null } = {}) {
  const n = x.length;
  if (n < 16) throw Object.assign(new Error('El espectro tiene muy pocos puntos.'), { status: 400 });
  const isH = nucleus === '1H';
  const base = baseline(y);
  const yc = new Float64Array(n);
  for (let i = 0; i < n; i++) yc[i] = y[i] - base[i];
  const sigma = noiseLevel(yc);
  const dxPpm = Math.abs(x[n - 1] - x[0]) / (n - 1);
  const hzPerPt = dxPpm * frequency;

  // Rango químico razonable y señal de referencia (TMS) fuera del análisis
  const lo = isH ? -1 : -10;
  const hi = isH ? 16 : 240;
  const tms = isH ? 0.04 : 0.4;
  const solv = SOLVENTS[solvent];
  const solventShifts = solv ? (isH ? [...solv.H, ...(solv.water ? [solv.water] : [])] : solv.C) : [];
  const nearSolvent = (ppm) => solventShifts.some((s) => Math.abs(s - ppm) < (isH ? 0.03 : 0.6));

  // Máximos locales
  const maxima = [];
  for (let i = 1; i < n - 1; i++) {
    if (yc[i] > yc[i - 1] && yc[i] >= yc[i + 1] && yc[i] > 3 * sigma && x[i] >= lo && x[i] <= hi) maxima.push(i);
  }
  // Altura de referencia: el máximo fuera del disolvente y del TMS
  let ref = 0;
  for (const i of maxima) if (!nearSolvent(x[i]) && Math.abs(x[i]) > tms && yc[i] > ref) ref = yc[i];
  if (!ref) for (const i of maxima) ref = Math.max(ref, yc[i]);
  if (!ref) throw Object.assign(new Error('No se detectaron señales por encima del ruido en el espectro.'), { status: 400 });

  // Prominencia: altura sobre el más alto de los dos mínimos que la separan de un pico mayor
  const prominence = (i) => {
    const h = yc[i];
    let minL = h;
    for (let k = i - 1; k >= 0; k--) {
      if (yc[k] > h) break;
      if (yc[k] < minL) minL = yc[k];
      if (yc[k] < 0) break;
    }
    let minR = h;
    for (let k = i + 1; k < n; k++) {
      if (yc[k] > h) break;
      if (yc[k] < minR) minR = yc[k];
      if (yc[k] < 0) break;
    }
    return h - Math.max(minL, minR);
  };

  let rel = isH ? 0.008 : 0.04;
  const kSigma = isH ? 8 : 6;
  let idx = [];
  for (let attempt = 0; attempt < 8; attempt++) {
    const thr = Math.max(kSigma * sigma, rel * ref);
    idx = maxima.filter((i) => yc[i] >= thr && prominence(i) >= Math.max(4 * sigma, 0.02 * yc[i]));
    if (idx.length <= (isH ? 500 : 300)) break;
    rel *= 1.6;
  }
  if (!idx.length) throw Object.assign(new Error('No se detectaron picos en el espectro.'), { status: 400 });

  // Posición por interpolación parabólica
  const pos = (i) => {
    const a = yc[i - 1];
    const b = yc[i];
    const c = yc[i + 1];
    const d = a - 2 * b + c;
    const off = d ? (0.5 * (a - c)) / d : 0;
    return x[i] + Math.max(-0.5, Math.min(0.5, off)) * (x[i + 1] - x[i]);
  };

  // Límites de integración: mínimo entre picos vecinos; hacia afuera, hasta caer al nivel del ruido
  const bounds = idx.map((i, k) => {
    const level = (h) => Math.max(2 * sigma, 0.002 * h);
    let l = i;
    const lv = level(yc[i]);
    const leftStop = k > 0 ? argmin(yc, idx[k - 1], i) : 0;
    while (l > leftStop && yc[l - 1] > lv) l--;
    let r = i;
    const rightStop = k < idx.length - 1 ? argmin(yc, i, idx[k + 1]) : n - 1;
    while (r < rightStop && yc[r + 1] > lv) r++;
    return [l, r];
  });

  let peaks = idx.map((i, k) => {
    const [l, r] = bounds[k];
    let area = 0;
    for (let j = l; j <= r; j++) area += yc[j];
    area *= hzPerPt;
    // Anchura a media altura (Hz)
    const half = yc[i] / 2;
    let a = i;
    while (a > l && yc[a] > half) a--;
    let b = i;
    while (b < r && yc[b] > half) b++;
    return { ppm: pos(i), height: yc[i], area: Math.max(area, 0), width: Math.max(hzPerPt, (b - a) * hzPerPt), solvent: nearSolvent(x[i]) };
  });

  // Máximos casi coincidentes (ruido sobre la cima de una línea): se fusionan
  const minSepPpm = (isH ? 0.4 : 2) / frequency;
  const merged = [];
  for (const p of [...peaks].sort((a, b) => a.ppm - b.ppm)) {
    const q = merged[merged.length - 1];
    if (q && p.ppm - q.ppm < minSepPpm) {
      if (p.height > q.height) Object.assign(q, { ppm: p.ppm, height: p.height, width: p.width });
      q.area += p.area;
    } else merged.push({ ...p });
  }
  peaks = merged;
  // ¹³C: líneas mínimas pegadas a una línea mucho mayor (ondulaciones de truncamiento, satélites)
  if (!isH) {
    const nb = 25 / frequency;
    const keep = peaks.filter((p) => !peaks.some((q) => q !== p && Math.abs(q.ppm - p.ppm) < nb && q.height > 8 * p.height));
    peaks = keep;
  }

  // Señales anchas (OH, NH, COOH): poca altura pero mucha área. Se buscan en la curva suavizada
  // (≈ 12 Hz) fuera de las zonas ya integradas.
  const covered = new Uint8Array(n);
  for (const [l, r] of bounds) covered.fill(1, l, r + 1);
  const win = Math.max(3, Math.round(12 / hzPerPt) | 1);
  if (isH && win < n / 4) {
    const cs = new Float64Array(n + 1);
    for (let i = 0; i < n; i++) cs[i + 1] = cs[i] + yc[i];
    const sm = new Float64Array(n);
    const h = win >> 1;
    for (let i = 0; i < n; i++) {
      const a = Math.max(0, i - h);
      const b = Math.min(n, i + h + 1);
      sm[i] = (cs[b] - cs[a]) / (b - a);
    }
    const sSigma = sigma / Math.sqrt(win);
    const totalArea = peaks.reduce((s, p) => s + p.area, 0) || 1;
    for (let i = 1; i < n - 1; i++) {
      if (covered[i] || !(sm[i] > sm[i - 1] && sm[i] >= sm[i + 1]) || sm[i] < Math.max(10 * sSigma, 0.002 * ref)) continue;
      if (x[i] < lo || x[i] > hi) continue;
      let l = i;
      while (l > 0 && sm[l - 1] > 2 * sSigma && !covered[l - 1]) l--;
      let r = i;
      while (r < n - 1 && sm[r + 1] > 2 * sSigma && !covered[r + 1]) r++;
      let area = 0;
      for (let j = l; j <= r; j++) area += yc[j];
      area *= hzPerPt;
      const widthHz = (r - l) * hzPerPt;
      if (area / totalArea >= 0.02 && widthHz >= 15) {
        peaks.push({ ppm: x[i], height: sm[i], area, width: widthHz / 2, solvent: nearSolvent(x[i]), broad: true });
        covered.fill(1, l, r + 1);
      }
      i = r;
    }
  }

  const notes = [];
  const broad = peaks.filter((p) => p.broad);
  if (broad.length) notes.push(`Se detectaron ${broad.length} señal(es) ancha(s) de poca altura (δ ${broad.map((p) => p.ppm.toFixed(2)).join(', ')}): posibles OH/NH/COOH intercambiables.`);
  const before = peaks.length;
  peaks = peaks.filter((p) => Math.abs(p.ppm) > tms);
  if (peaks.length < before) notes.push('Se omitió la señal de referencia (TMS, δ 0).');

  // ¹H: grupos con área despreciable (satélites de ¹³C, bandas de rotación, impurezas mínimas)
  if (isH) {
    const gap = 15 / frequency;
    const sorted = [...peaks].sort((p, q) => q.ppm - p.ppm);
    const groups = [];
    for (const p of sorted) {
      const g = groups[groups.length - 1];
      if (g && g[g.length - 1].ppm - p.ppm <= gap) g.push(p);
      else groups.push([p]);
    }
    const total = groups.filter((g) => !g.every((p) => p.solvent)).reduce((a, g) => a + g.reduce((s, p) => s + p.area, 0), 0) || 1;
    const keep = groups.filter((g) => g.every((p) => p.solvent) || g.reduce((s, p) => s + p.area, 0) / total >= 0.012);
    const dropped = groups.length - keep.length;
    if (dropped) notes.push(`Se descartaron ${dropped} grupo(s) de líneas con área < 1.2 % del total (satélites de ¹³C, bandas laterales o impurezas).`);
    peaks = keep.flat();
  }
  peaks.sort((p, q) => q.ppm - p.ppm);
  return { peaks, sigma, baseline: base, corrected: yc, snr: ref / sigma, ref, notes };
}

function argmin(arr, a, b) {
  let m = a;
  for (let i = a; i <= b; i++) if (arr[i] < arr[m]) m = i;
  return m;
}

/**
 * Curva reducida para guardar y dibujar: rango químico útil, a lo sumo `maxPoints` puntos
 * (en cada tramo se conserva el valor de mayor |y| para no perder picos), escala 0–100.
 */
export function compactSpectrum(x, y, { nucleus = '1H', ref = null, maxPoints = 32768 } = {}) {
  const isH = nucleus === '1H';
  const lo = Math.max(x[0], isH ? -1 : -10);
  const hi = Math.min(x[x.length - 1], isH ? 14 : 230);
  let a = 0;
  while (a < x.length - 1 && x[a] < lo) a++;
  let b = x.length - 1;
  while (b > a && x[b] > hi) b--;
  const count = b - a + 1;
  const buckets = Math.min(maxPoints, count);
  const scale = 100 / (ref || Math.max(...y.subarray(a, b + 1)) || 1);
  const out = new Array(buckets);
  for (let k = 0; k < buckets; k++) {
    const s = a + Math.floor((k * count) / buckets);
    const e = a + Math.floor(((k + 1) * count) / buckets);
    let v = y[s];
    for (let i = s + 1; i < e; i++) if (Math.abs(y[i]) > Math.abs(v)) v = y[i];
    out[k] = Math.round(v * scale * 1000) / 1000;
  }
  return { from: x[a], to: x[b], y: out };
}

// Constante de acoplamiento ¹J(C–D) del disolvente deuterado (Hz)
const JCD = { CDCl3: 32.0, 'DMSO-d6': 21.0, CD3OD: 21.4, 'Acetona-d6': 19.4, C6D6: 24.3, CD3CN: 21.0, CD2Cl2: 27.2 };

/**
 * ¹³C: busca el multiplete del disolvente deuterado (línea central con sus dos vecinas a ±J(C–D)
 * y de altura comparable) cerca de su δ esperado. Si está desplazado ≥ 0.25 ppm, el espectro está
 * mal referenciado: devuelve la corrección para re-referenciarlo al disolvente.
 */
export function referenceCheck(peaks, { nucleus, solvent, frequency }) {
  const solv = SOLVENTS[solvent];
  if (nucleus !== '13C' || !solv?.C?.length || !JCD[solvent] || !frequency) return null;
  const expected = solv.C[0];
  const j = JCD[solvent] / frequency;
  const tol = 0.25 * j;
  const near = peaks.filter((p) => Math.abs(p.ppm - expected) < 2);
  const centers = near.filter((c) => {
    const side = (s) => near.find((p) => Math.abs(p.ppm - (c.ppm + s * j)) < tol && p.height >= 0.3 * c.height && p.height <= 3 * c.height);
    return side(1) && side(-1);
  });
  if (!centers.length) return null;
  const center = centers.reduce((a, c) => (Math.abs(c.ppm - expected) < Math.abs(a.ppm - expected) ? c : a)).ppm;
  const diff = center - expected;
  if (Math.abs(diff) < 0.25) return null;
  return {
    diff,
    center,
    note: `El multiplete del disolvente (${solvent}) aparecía en δ ${center.toFixed(2)} en lugar de ${expected.toFixed(2)}: el espectro estaba mal referenciado y se re-referenció al disolvente (corrección de ${diff > 0 ? '−' : '+'}${Math.abs(diff).toFixed(2)} ppm en todos los picos).`,
  };
}
