// Exportación de espectros: JCAMP-DX 5.01, CSV y JSON.
import { spectrumXY } from './render.js';

const NUC = { '1H': '^1H', '13C': '^13C' };

function header(s, c, dataType) {
  const lines = [
    `##TITLE=${c?.name || c?.iupac || c?.formula || 'Compuesto'} - ${s.nucleus} NMR (${s.kind})`,
    '##JCAMP-DX=5.01',
    '##DATA TYPE=NMR SPECTRUM',
    `##DATA CLASS=${dataType}`,
    `##ORIGIN=${s.source}`,
    '##OWNER=Analizador RMN',
    `##LONGDATE=${new Date().toISOString()}`,
    `##.OBSERVE NUCLEUS=${NUC[s.nucleus] || s.nucleus}`,
    `##.OBSERVE FREQUENCY=${s.frequency || ''}`,
    `##.SOLVENT NAME=${s.solvent || ''}`,
  ];
  if (c?.formula) lines.push(`##MOLFORM=${c.formula.replace(/([A-Z][a-z]?)(\d*)/g, (m, el, n) => `${el} ${n || 1} `).trim()}`);
  if (c?.smiles) lines.push(`##$SMILES=${c.smiles}`);
  if (c?.inchikey) lines.push(`##$INCHIKEY=${c.inchikey}`);
  if (c?.cid) lines.push(`##$PUBCHEM CID=${c.cid}`);
  // Señales con multiplicidad, J e integración como bloque de usuario
  for (const g of s.signals || []) {
    const J = g.js?.length ? ` J=${g.js.map((j) => j.J).join(',')}Hz` : '';
    const mult = g.multiplicity || g.type || '';
    const n = g.nH ? ` ${g.nH}H` : g.nC ? ` ${g.nC}C` : '';
    const at = g.atoms?.length ? ` atomos=${g.atoms.join(',')}` : '';
    lines.push(`##$SIGNAL=${Number(g.delta).toFixed(4)} ${mult}${J}${n}${at}`.trimEnd());
  }
  return lines;
}

export function toJcamp(s, c) {
  const xy = spectrumXY(s);
  if (!xy) return toJcampPeaks(s, c);
  const n = xy.x.length;
  let ymax = 0;
  for (let i = 0; i < n; i++) ymax = Math.max(ymax, Math.abs(xy.y[i]));
  const yFactor = ymax / 1e6 || 1;
  const deltaX = (xy.x[n - 1] - xy.x[0]) / (n - 1);
  const lines = header(s, c, 'XYDATA');
  lines.push(
    '##XUNITS=PPM',
    '##YUNITS=ARBITRARY UNITS',
    '##XFACTOR=1',
    `##YFACTOR=${yFactor.toExponential(8)}`,
    `##FIRSTX=${xy.x[0]}`,
    `##LASTX=${xy.x[n - 1]}`,
    `##DELTAX=${deltaX}`,
    `##NPOINTS=${n}`,
    `##FIRSTY=${Math.round(xy.y[0] / yFactor) * yFactor}`,
    '##XYDATA=(X++(Y..Y))',
  );
  for (let i = 0; i < n; i += 10) {
    const ys = [];
    for (let k = i; k < Math.min(n, i + 10); k++) ys.push(Math.round(xy.y[k] / yFactor));
    lines.push(`${xy.x[i].toFixed(6)} ${ys.join(' ')}`);
  }
  lines.push('##END=');
  return lines.join('\n') + '\n';
}

export function toJcampPeaks(s, c) {
  const lines = header(s, c, 'PEAK TABLE');
  lines.push('##XUNITS=PPM', '##YUNITS=RELATIVE INTENSITY', `##NPOINTS=${(s.peaks || []).length}`, '##PEAK TABLE=(XY..XY)');
  for (const p of s.peaks || []) lines.push(`${p.ppm.toFixed(4)}, ${p.intensity}`);
  lines.push('##END=');
  return lines.join('\n') + '\n';
}

const csvCell = (v) => {
  const t = Array.isArray(v) ? v.join(' ') : v === null || v === undefined ? '' : String(v);
  return /[",;\n]/.test(t) ? `"${t.replace(/"/g, '""')}"` : t;
};

export function toCsv(s, c) {
  const rows = [];
  rows.push(['# Compuesto', c?.name || c?.iupac || ''], ['# Formula', c?.formula || ''], ['# SMILES', c?.smiles || '']);
  rows.push(['# Nucleo', s.nucleus], ['# Tipo', s.kind], ['# Fuente', s.source], ['# Frecuencia (MHz)', s.frequency], ['# Disolvente', s.solvent || '']);
  rows.push([]);
  rows.push(['SEÑALES']);
  if (s.nucleus === '1H' && s.kind === 'predicho') {
    rows.push(['delta_ppm', 'multiplicidad', 'J_Hz', 'nH', 'integral_relativa', 'atomos', 'entorno']);
    for (const g of s.signals || [])
      rows.push([g.delta, g.multiplicity, (g.js || []).map((j) => j.J).join(' / '), g.nH, g.relativeIntegration, (g.atoms || []).join(' '), g.environment]);
  } else if (s.kind === 'predicho') {
    rows.push(['delta_ppm', 'tipo', 'nC', 'intensidad_relativa', 'DEPT135', 'atomos', 'entorno']);
    for (const g of s.signals || []) rows.push([g.delta, g.type, g.nC, g.intensity, g.dept135, (g.atoms || []).join(' '), g.environment]);
  } else {
    rows.push(['delta_ppm', 'tipo/multiplicidad', 'intensidad', 'atomos']);
    for (const g of s.signals || []) rows.push([g.delta, g.type || g.multiplicity || '', g.intensity ?? g.count ?? '', (g.atoms || []).join(' ')]);
  }
  rows.push([]);
  rows.push(['LÍNEAS / PICOS']);
  rows.push(['delta_ppm', 'Hz', 'intensidad_relativa_%', 'señal']);
  for (const p of s.peaks || []) rows.push([p.ppm, p.hz ?? (s.frequency ? (p.ppm * s.frequency).toFixed(2) : ''), p.intensity, p.signal ?? '']);
  return rows.map((r) => r.map(csvCell).join(',')).join('\n') + '\n';
}

export function toXyCsv(s) {
  const xy = spectrumXY(s);
  if (!xy) return 'ppm,intensidad\n';
  const out = ['ppm,intensidad'];
  const step = Math.max(1, Math.floor(xy.x.length / 32768));
  for (let i = 0; i < xy.x.length; i += step) out.push(`${xy.x[i].toFixed(6)},${xy.y[i].toExponential(6)}`);
  return out.join('\n') + '\n';
}
