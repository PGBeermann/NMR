// Dibuja espectros de RMN como imagen SVG (y PNG mediante resvg).
import fs from 'node:fs';
import path from 'node:path';
import { Resvg } from '@resvg/resvg-js';

const esc = (s) =>
  String(s ?? '').replace(/[&<>"']/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' })[c]);

function niceStep(span, target) {
  const raw = span / target;
  const steps = [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.25, 0.5, 1, 2, 5, 10, 20, 25, 50, 100];
  return steps.find((s) => s >= raw) || 100;
}

function decimalsFor(step) {
  if (step >= 1) return 0;
  return Math.min(3, Math.ceil(-Math.log10(step) - 1e-9));
}

/** Reparte etiquetas para que no se solapen (mínimo `gap` px entre ellas). */
function spreadLabels(items, gap, minX, maxX) {
  const s = [...items].sort((a, b) => a.x - b.x);
  for (let i = 0; i < s.length; i++) s[i].lx = s[i].x;
  for (let iter = 0; iter < 50; iter++) {
    let moved = false;
    for (let i = 1; i < s.length; i++) {
      const d = s[i].lx - s[i - 1].lx;
      if (d < gap) {
        const push = (gap - d) / 2;
        s[i - 1].lx -= push;
        s[i].lx += push;
        moved = true;
      }
    }
    for (const it of s) it.lx = Math.min(maxX, Math.max(minX, it.lx));
    if (!moved) break;
  }
  return s;
}

/**
 * @param {object} o
 * @param {{x: Float64Array, y: Float64Array}} [o.xy] curva continua
 * @param {Array<{ppm:number,intensity:number}>} [o.sticks] líneas verticales
 * @param {Array<{ppm:number,text:string,height?:number}>} [o.labels]
 * @param {Array<{from:number,to:number,text:string}>} [o.integrals]
 * @param {Array<{ppm:number,text:string}>} [o.marks] disolvente / TMS
 */
export function spectrumSvg(o) {
  const W = o.width || 1100;
  const H = o.height || 500;
  const from = o.from;
  const to = o.to;
  const color = o.color || '#1f4e9c';
  const labelBand = o.labels?.length ? 92 : 20;
  const integralBand = o.integrals?.length ? 30 : 0;
  const m = { l: 36, r: 30, t: 48 + labelBand, b: 52 + integralBand };
  const pw = W - m.l - m.r;
  const ph = H - m.t - m.b;
  const base = m.t + ph;
  const X = (ppm) => m.l + ((to - ppm) / (to - from)) * pw;

  // Máximo de intensidad en el rango visible
  let ymax = 0;
  if (o.xy) {
    for (let i = 0; i < o.xy.x.length; i++) {
      const xv = o.xy.x[i];
      if (xv >= from && xv <= to && o.xy.y[i] > ymax) ymax = o.xy.y[i];
    }
  }
  if (o.sticks) for (const s of o.sticks) if (s.ppm >= from && s.ppm <= to && s.intensity > ymax) ymax = s.intensity;
  if (!ymax) ymax = 1;
  const Y = (v) => base - (v / ymax) * ph * 0.94;

  const parts = [];
  parts.push(
    `<svg xmlns="http://www.w3.org/2000/svg" width="${W}" height="${H}" viewBox="0 0 ${W} ${H}" font-family="DejaVu Sans, Arial, Helvetica, sans-serif">`,
    `<rect width="${W}" height="${H}" fill="#ffffff"/>`,
    `<text x="${m.l}" y="22" font-size="16" font-weight="bold" fill="#1d1d1f">${esc(o.title)}</text>`,
    `<text x="${m.l}" y="40" font-size="12" fill="#555">${esc(o.subtitle || '')}</text>`,
    `<text x="${W - m.r}" y="22" font-size="12" fill="#555" text-anchor="end">${esc(o.info || '')}</text>`,
  );

  // Marcas del disolvente y TMS
  for (const mk of o.marks || []) {
    if (mk.ppm < from || mk.ppm > to) continue;
    const x = X(mk.ppm).toFixed(1);
    parts.push(
      `<line x1="${x}" y1="${m.t}" x2="${x}" y2="${base}" stroke="#bbb" stroke-dasharray="4 4"/>`,
      `<text x="${x}" y="${m.t + 12}" font-size="10" fill="#999" text-anchor="middle">${esc(mk.text)}</text>`,
    );
  }

  // Curva continua (decimada a mín/máx por píxel para conservar los picos)
  if (o.xy) {
    const cols = new Map();
    for (let i = 0; i < o.xy.x.length; i++) {
      const xv = o.xy.x[i];
      if (xv < from || xv > to) continue;
      const px = Math.round(X(xv) * 2) / 2;
      const yv = o.xy.y[i];
      const c = cols.get(px);
      if (!c) cols.set(px, { min: yv, max: yv });
      else {
        if (yv < c.min) c.min = yv;
        if (yv > c.max) c.max = yv;
      }
    }
    const keys = [...cols.keys()].sort((a, b) => a - b);
    let d = '';
    for (const k of keys) {
      const c = cols.get(k);
      d += `${d ? 'L' : 'M'}${k.toFixed(1)},${Y(c.min).toFixed(1)}`;
      if (c.max !== c.min) d += `L${k.toFixed(1)},${Y(c.max).toFixed(1)}`;
    }
    parts.push(`<path d="${d}" fill="none" stroke="${color}" stroke-width="1.1" stroke-linejoin="round"/>`);
  }
  if (o.sticks && !o.xy) {
    for (const s of o.sticks) {
      if (s.ppm < from || s.ppm > to) continue;
      const x = X(s.ppm).toFixed(1);
      parts.push(`<line x1="${x}" y1="${base}" x2="${x}" y2="${Y(s.intensity).toFixed(1)}" stroke="${color}" stroke-width="1.6"/>`);
    }
  }

  // Curvas de integración (1H)
  if (o.integralCurve) {
    const { x, y } = o.integralCurve;
    let d = '';
    let yMaxI = 0;
    for (let i = 0; i < x.length; i++) if (x[i] >= from && x[i] <= to) yMaxI = Math.max(yMaxI, y[i]);
    let yMinI = Infinity;
    for (let i = 0; i < x.length; i++) if (x[i] >= from && x[i] <= to) yMinI = Math.min(yMinI, y[i]);
    const span = yMaxI - yMinI || 1;
    const step = Math.max(1, Math.floor(x.length / 3000));
    for (let i = 0; i < x.length; i += step) {
      if (x[i] < from || x[i] > to) continue;
      const yy = base - 6 - ((y[i] - yMinI) / span) * ph * 0.55;
      d += `${d ? 'L' : 'M'}${X(x[i]).toFixed(1)},${yy.toFixed(1)}`;
    }
    parts.push(`<path d="${d}" fill="none" stroke="#2e9d57" stroke-width="1" opacity="0.75"/>`);
  }

  // Línea base
  parts.push(`<line x1="${m.l}" y1="${base}" x2="${W - m.r}" y2="${base}" stroke="#ccc"/>`);

  // Etiquetas de picos
  if (o.labels?.length) {
    const visible = o.labels.filter((l) => l.ppm >= from && l.ppm <= to);
    const maxLabels = Math.floor(pw / 13);
    const chosen = visible.length > maxLabels ? [...visible].sort((a, b) => (b.height || 0) - (a.height || 0)).slice(0, maxLabels) : visible;
    const items = spreadLabels(chosen.map((l) => ({ ...l, x: X(l.ppm) })), 13, m.l + 6, W - m.r - 6);
    const bandBottom = m.t - 6;
    for (const it of items) {
      const peakY = it.height !== undefined ? Y(it.height) - 4 : base - 10;
      const ly = bandBottom - 14;
      parts.push(
        `<path d="M${it.x.toFixed(1)},${Math.max(peakY, bandBottom).toFixed(1)} L${it.x.toFixed(1)},${bandBottom} L${it.lx.toFixed(1)},${ly}" fill="none" stroke="#999" stroke-width="0.7"/>`,
        `<text transform="translate(${(it.lx + 3.5).toFixed(1)},${(ly - 3).toFixed(1)}) rotate(-90)" font-size="10.5" fill="#222">${esc(it.text)}</text>`,
      );
    }
  }

  // Integrales bajo la línea base
  if (o.integrals?.length) {
    const yI = base + 10;
    for (const it of o.integrals) {
      if (it.to < from || it.from > to) continue;
      const x1 = X(Math.min(to, it.to));
      const x2 = X(Math.max(from, it.from));
      const xm = (x1 + x2) / 2;
      parts.push(
        `<path d="M${x1.toFixed(1)},${yI - 4} L${x1.toFixed(1)},${yI} L${x2.toFixed(1)},${yI} L${x2.toFixed(1)},${yI - 4}" fill="none" stroke="#2e9d57" stroke-width="1.2"/>`,
        `<text x="${xm.toFixed(1)}" y="${yI + 13}" font-size="10.5" fill="#2e7d4f" text-anchor="middle">${esc(it.text)}</text>`,
      );
    }
  }

  // Eje de desplazamiento químico
  const axisY = base + integralBand + 8;
  parts.push(`<line x1="${m.l}" y1="${axisY}" x2="${W - m.r}" y2="${axisY}" stroke="#333"/>`);
  const step = niceStep(to - from, pw / 70);
  const dec = decimalsFor(step);
  const minor = step / 5;
  for (let v = Math.ceil(from / minor) * minor; v <= to + 1e-9; v += minor) {
    const x = X(v).toFixed(1);
    const isMajor = Math.abs(v / step - Math.round(v / step)) < 1e-6;
    parts.push(`<line x1="${x}" y1="${axisY}" x2="${x}" y2="${axisY + (isMajor ? 6 : 3)}" stroke="#333"/>`);
    if (isMajor) parts.push(`<text x="${x}" y="${axisY + 19}" font-size="11" text-anchor="middle" fill="#222">${(Math.round(v / step) * step).toFixed(dec)}</text>`);
  }
  parts.push(`<text x="${(m.l + pw / 2).toFixed(1)}" y="${axisY + 37}" font-size="12" text-anchor="middle" fill="#222">δ (ppm)</text>`);
  if (o.footer) parts.push(`<text x="${W - m.r}" y="${H - 6}" font-size="9.5" fill="#888" text-anchor="end">${esc(o.footer)}</text>`);
  parts.push('</svg>');
  return parts.join('\n');
}

// Archivos de fuentes conocidos (en Docker: paquete fonts-dejavu-core).
// Se pasan explícitamente porque resvg no siempre detecta las fuentes del sistema.
let fontFiles = null;
function findFonts() {
  if (fontFiles) return fontFiles;
  fontFiles = [];
  const dirs = ['/usr/share/fonts', '/usr/local/share/fonts', 'C:/Windows/Fonts'];
  const walk = (dir, depth = 0) => {
    let entries = [];
    try {
      entries = fs.readdirSync(dir, { withFileTypes: true });
    } catch {
      return;
    }
    for (const e of entries) {
      const p = path.join(dir, e.name);
      if (e.isDirectory() && depth < 4) walk(p, depth + 1);
      else if (/^(DejaVuSans(-Bold)?|arial|arialbd)\.ttf$/i.test(e.name)) fontFiles.push(p);
    }
  };
  for (const d of dirs) walk(d);
  return fontFiles;
}

export function svgToPng(svg, width = 2200) {
  const files = findFonts();
  const family = files.some((f) => /DejaVu/i.test(f)) ? 'DejaVu Sans' : 'Arial';
  const r = new Resvg(svg, {
    fitTo: { mode: 'width', value: width },
    background: 'white',
    font: { fontFiles: files, loadSystemFonts: files.length === 0, defaultFontFamily: family, sansSerifFamily: family },
  });
  return r.render().asPng();
}
