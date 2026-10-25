// Generación de espectros (curvas x/y) a partir de señales o listas de picos.
import { signalsToXY, xyAutoPeaksPicking } from 'nmr-processing';

export const carbonFrequency = (protonMHz) => Math.round(protonMHz * 0.251449 * 100) / 100;

const xyCache = new Map();

function remember(key, value) {
  if (xyCache.size > 60) xyCache.delete(xyCache.keys().next().value);
  xyCache.set(key, value);
  return value;
}

export function protonRange(deltas) {
  const min = Math.min(...deltas, 0);
  const max = Math.max(...deltas, 9);
  return { from: Math.floor(min - 0.5), to: Math.ceil(max + 0.8) };
}

/**
 * Simulación completa del espectro de 1H (incluye efectos de segundo orden)
 * a partir de las señales por átomo con sus constantes de acoplamiento.
 */
export function simulateProton(rawSignals, { frequency = 400, from, to, lineWidth = 1, nbPoints = 65536, cacheKey } = {}) {
  if (cacheKey && xyCache.has(cacheKey)) return xyCache.get(cacheKey);
  const xy = signalsToXY(rawSignals, { frequency, from, to, nbPoints, lineWidth, maxClusterSize: 9 });
  const out = { x: Float64Array.from(xy.x), y: Float64Array.from(xy.y) };
  return cacheKey ? remember(cacheKey, out) : out;
}

/** Lista de líneas (picos) del espectro simulado con intensidades relativas (0–100). */
export function pickLines(xy, frequency) {
  const peaks = xyAutoPeaksPicking(xy, { frequency, thresholdFactor: 1, minMaxRatio: 0.005 });
  const max = Math.max(...peaks.map((p) => p.y), 1);
  return peaks
    .map((p) => ({
      ppm: Number(p.x.toFixed(4)),
      hz: Number((p.x * frequency).toFixed(2)),
      intensity: Number(((p.y / max) * 100).toFixed(2)),
    }))
    .filter((p) => p.intensity >= 0.5)
    .sort((a, b) => b.ppm - a.ppm);
}

/** Suma de lorentzianas para listas de picos (espectros de 13C o experimentales). */
export function sticksToXY(sticks, { from, to, nbPoints = 32768, width = 0.02, cacheKey } = {}) {
  if (cacheKey && xyCache.has(cacheKey)) return xyCache.get(cacheKey);
  const x = new Float64Array(nbPoints);
  const y = new Float64Array(nbPoints);
  const step = (to - from) / (nbPoints - 1);
  for (let i = 0; i < nbPoints; i++) x[i] = from + i * step;
  const hw = width / 2;
  const span = Math.ceil((hw * 60) / step);
  for (const s of sticks) {
    const c = Math.round((s.ppm - from) / step);
    for (let i = Math.max(0, c - span); i <= Math.min(nbPoints - 1, c + span); i++) {
      const d = x[i] - s.ppm;
      y[i] += (s.intensity * hw * hw) / (d * d + hw * hw);
    }
  }
  const out = { x, y };
  return cacheKey ? remember(cacheKey, out) : out;
}

/** Curva integral acumulada (para dibujarla sobre el espectro de 1H). */
export function cumulativeIntegral(xy) {
  const out = new Float64Array(xy.y.length);
  // Se integra desde δ alto (izquierda del espectro) hacia δ bajo
  let s = 0;
  for (let i = xy.y.length - 1; i >= 0; i--) {
    s += xy.y[i];
    out[i] = s;
  }
  return out;
}
