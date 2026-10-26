// Construye la imagen y los datos x/y de un espectro guardado.
import { SOLVENTS } from './chemistry.js';
import { spectrumSvg } from './plot.js';
import { cumulativeIntegral, protonRange, simulateProton, sticksToXY } from './spectrum.js';

const NUC = { '1H': '¹H', '13C': '¹³C' };

export function spectrumXY(s) {
  if (s.kind === 'predicho' && s.nucleus === '1H' && s.meta?.rawSignals) {
    const range = s.meta.range || protonRange(s.signals.map((x) => x.delta));
    return simulateProton(s.meta.rawSignals, {
      frequency: s.frequency,
      ...range,
      lineWidth: s.meta.lineWidth || 1,
      cacheKey: `h-${s.id}-${s.updatedAt || ''}`,
    });
  }
  const peaks = s.peaks || [];
  if (!peaks.length) return null;
  const ppms = peaks.map((p) => p.ppm);
  const isH = s.nucleus === '1H';
  const from = isH ? Math.floor(Math.min(0, ...ppms) - 0.5) : Math.floor(Math.min(0, ...ppms) - 5);
  const to = isH ? Math.ceil(Math.max(9, ...ppms) + 0.8) : Math.ceil(Math.max(200, ...ppms) + 10);
  const width = isH ? 1.2 / (s.frequency || 400) : 0.06;
  return sticksToXY(peaks, { from, to, width, nbPoints: isH ? 65536 : 32768, cacheKey: `s-${s.id}-${s.updatedAt || ''}` });
}

export function defaultRange(s, xy) {
  if (s.kind === 'predicho' && s.nucleus === '1H' && s.meta?.range) return s.meta.range;
  if (xy) return { from: xy.x[0], to: xy.x[xy.x.length - 1] };
  return s.nucleus === '1H' ? { from: -0.5, to: 10 } : { from: -5, to: 220 };
}

/** Espectro cargado por el usuario (lista de picos CSV). */
export function uploadXY(u) {
  const isH = u.nucleus === '1H';
  const ppms = u.peaks.map((p) => p.ppm);
  const from = isH ? Math.floor(Math.min(0, ...ppms) - 0.5) : Math.floor(Math.min(0, ...ppms) - 5);
  const to = isH ? Math.ceil(Math.max(9, ...ppms) + 0.8) : Math.ceil(Math.max(200, ...ppms) + 10);
  const width = isH ? 0.9 / (u.frequency || 400) : 0.06;
  return sticksToXY(u.peaks, { from, to, width, nbPoints: isH ? 131072 : 32768, cacheKey: `u-${u.id}-${u.updatedAt}` });
}

export function renderUpload(u, { from, to, width, height } = {}) {
  const xy = uploadXY(u);
  const isH = u.nucleus === '1H';
  const sig = u.analysis?.signals || [];
  let lo = Number.isFinite(from) ? from : xy.x[0];
  let hi = Number.isFinite(to) ? to : xy.x[xy.x.length - 1];
  if (lo > hi) [lo, hi] = [hi, lo];
  if (hi - lo < 0.01) hi = lo + 0.01;
  const heightAt = (a, b) => {
    let m = 0;
    for (let i = 0; i < xy.x.length; i++) if (xy.x[i] >= a && xy.x[i] <= b && xy.y[i] > m) m = xy.y[i];
    return m;
  };
  const labels = isH
    ? sig.map((g) => ({ ppm: g.delta, tag: g.label || null, text: `${g.delta.toFixed(2)} (${g.multiplicity}, ${g.nH}H)`, height: heightAt(g.from - 0.002, g.to + 0.002) }))
    : sig.map((g) => ({ ppm: g.delta, text: g.delta.toFixed(1), height: heightAt(g.delta - 0.2, g.delta + 0.2) }));
  const integrals = isH ? sig.map((g) => ({ from: g.from - 0.004, to: g.to + 0.004, text: `${g.nH}H` })) : [];
  return spectrumSvg({
    width: width || 1100,
    height: height || 500,
    from: lo,
    to: hi,
    title: `${NUC[u.nucleus] || u.nucleus} RMN — ${u.name || 'Espectro cargado'}`,
    subtitle: `Espectro experimental cargado (CSV) · ${u.peaks.length} picos${u.formula ? ` · ${u.formula}` : ''}`,
    info: [u.frequency ? `${u.frequency} MHz` : null, u.solvent].filter(Boolean).join(' · '),
    xy,
    labels,
    integrals,
    integralCurve: isH ? { x: xy.x, y: cumulativeIntegral(xy) } : null,
    marks: [{ ppm: 0, text: 'TMS' }],
    color: '#b02a37',
    footer: 'Analizador RMN',
  });
}

export function renderSpectrum(s, compound, { from, to, width, height } = {}) {
  const xy = spectrumXY(s);
  const range = defaultRange(s, xy);
  const f = Number.isFinite(from) ? from : range.from;
  const t = Number.isFinite(to) ? to : range.to;
  const lo = Math.min(f, t);
  const hi = Math.max(f, t) === lo ? lo + 1 : Math.max(f, t);
  const isH = s.nucleus === '1H';
  const solv = SOLVENTS[s.solvent];

  const heightAt = (ppm) => {
    if (!xy) return 0;
    // Altura máxima de la curva cerca del desplazamiento
    const dx = isH ? 0.03 : 0.3;
    let m = 0;
    for (let i = 0; i < xy.x.length; i++) if (Math.abs(xy.x[i] - ppm) <= dx && xy.y[i] > m) m = xy.y[i];
    return m;
  };

  let labels = [];
  let integrals = [];
  let integralCurve = null;
  if (s.kind === 'predicho' && isH) {
    labels = s.signals.map((g) => ({ ppm: g.delta, tag: g.label || null, text: `${g.delta.toFixed(2)} (${g.multiplicity}, ${g.nH}H)`, height: heightAt(g.delta) }));
    integrals = (s.meta?.integrals || []).map((i) => ({ from: i.from, to: i.to, text: `${i.nH}H` }));
    if (xy) integralCurve = { x: xy.x, y: cumulativeIntegral(xy) };
  } else if (s.kind === 'predicho') {
    labels = s.signals.map((g) => ({ ppm: g.delta, text: `${g.delta.toFixed(1)} ${g.type}`, height: heightAt(g.delta) }));
  } else {
    labels = (s.peaks || []).map((p) => ({ ppm: p.ppm, text: p.ppm.toFixed(isH ? 2 : 1), height: heightAt(p.ppm) }));
  }

  const marks = [{ ppm: 0, text: 'TMS' }];
  if (solv && s.kind === 'predicho') for (const v of solv[isH ? 'H' : 'C']) marks.push({ ppm: v, text: s.solvent });

  const kindTxt = s.kind === 'predicho' ? 'Espectro predicho' : s.kind === 'calculado' ? 'Espectro calculado' : 'Espectro experimental';
  const svg = spectrumSvg({
    width: width || 1100,
    height: height || 500,
    from: lo,
    to: hi,
    title: `${NUC[s.nucleus] || s.nucleus} RMN — ${compound?.name || compound?.iupac || compound?.formula || ''}`,
    subtitle: `${kindTxt} · ${s.source}`,
    info: [s.frequency ? `${s.frequency} MHz` : null, s.solvent || null].filter(Boolean).join(' · '),
    xy,
    labels,
    integrals,
    integralCurve,
    marks,
    color: s.kind === 'predicho' ? '#1f4e9c' : s.kind === 'calculado' ? '#8a5a00' : '#b02a37',
    footer: `${compound?.formula || ''}  ·  Analizador RMN`,
  });
  return svg;
}
