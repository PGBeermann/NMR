// Lectura de archivos JCAMP-DX (4.24, 5.x y 6.0) de espectros de RMN de 1D:
// XYDATA=(X++(Y..Y)) comprimido (ASDF: SQZ, DIF, DUP, PAC) o AFFN, NTUPLES de Bruker/MestReNova
// (página de la parte real), XYPOINTS y PEAK TABLE. Convierte el eje a ppm con la referencia
// del archivo y obtiene núcleo, frecuencia y disolvente.

const fail = (msg) => Object.assign(new Error(msg), { status: 400 });

/** ¿El texto parece un archivo JCAMP-DX? */
export function isJcamp(text) {
  return /^\s*##\s*(TITLE|JCAMP)/im.test(String(text || '').slice(0, 4000));
}

const normLabel = (l) => l.replace(/[\s\-/_]/g, '').toUpperCase();

/** Divide el archivo en bloques (cada ##TITLE= abre un bloque; ##END= lo cierra). */
function splitBlocks(text) {
  const lines = String(text).replace(/\r/g, '').split('\n');
  const blocks = [];
  let cur = null;
  let lastLabel = null;
  for (let raw of lines) {
    const c = raw.indexOf('$$');
    if (c >= 0) raw = raw.slice(0, c);
    const line = raw.trimEnd();
    const m = /^\s*##([^=]*)=(.*)$/.exec(line);
    if (m) {
      const label = normLabel(m[1]);
      const value = m[2].trim();
      if (label === 'TITLE') {
        cur = { ldr: [], data: {} };
        blocks.push(cur);
      }
      if (!cur) {
        cur = { ldr: [], data: {} };
        blocks.push(cur);
      }
      if (label === 'END') {
        lastLabel = null;
        continue;
      }
      cur.ldr.push({ label, value, lines: [] });
      lastLabel = cur.ldr[cur.ldr.length - 1];
    } else if (lastLabel && line.trim()) {
      lastLabel.lines.push(line.trim());
    }
  }
  return blocks;
}

// --- ASDF -----------------------------------------------------------------------------
const SQZ = { '@': 0, A: 1, B: 2, C: 3, D: 4, E: 5, F: 6, G: 7, H: 8, I: 9, a: -1, b: -2, c: -3, d: -4, e: -5, f: -6, g: -7, h: -8, i: -9 };
const DIF = { '%': 0, J: 1, K: 2, L: 3, M: 4, N: 5, O: 6, P: 7, Q: 8, R: 9, j: -1, k: -2, l: -3, m: -4, n: -5, o: -6, p: -7, q: -8, r: -9 };
const DUP = { S: 1, T: 2, U: 3, V: 4, W: 5, X: 6, Y: 7, Z: 8, s: 9 };
const AFFN_LINE = /^[\s\d.+\-eE,;?]+$/;
const AFFN_NUM = /^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$/;

/**
 * Decodifica una línea de datos. Devuelve el valor de X (primer número), los valores de Y
 * y si la línea termina en forma DIF (en ese caso la línea siguiente repite el último Y).
 */
export function decodeLine(line) {
  // AFFN: números separados por espacios/comas (con decimales o exponentes). Lo demás
  // (incluido PAC, "123+45-6") se decodifica como ASDF.
  const tokens = line.trim().split(/[\s,;]+/).filter(Boolean);
  if (AFFN_LINE.test(line) && tokens.length >= 2 && tokens.every((t) => t === '?' || AFFN_NUM.test(t))) {
    const nums = tokens.map((t) => (t === '?' ? NaN : Number(t)));
    return { x: nums[0], ys: nums.slice(1), endsDif: false };
  }
  const vals = [];
  let tok = '';
  let kind = null;
  let lastDiff = 0;
  let lastWasDif = false;
  const flush = () => {
    if (!kind) return;
    if (kind === 'abs') {
      vals.push(Number(tok));
      lastWasDif = false;
    } else if (kind === 'dif') {
      lastDiff = Number(tok);
      vals.push(vals[vals.length - 1] + lastDiff);
      lastWasDif = true;
    } else if (kind === 'dup') {
      const n = Number(tok);
      for (let k = 1; k < n; k++) vals.push(lastWasDif ? vals[vals.length - 1] + lastDiff : vals[vals.length - 1]);
    }
    tok = '';
    kind = null;
  };
  for (const ch of line) {
    if (ch >= '0' && ch <= '9') {
      if (!kind) kind = 'abs';
      tok += ch;
    } else if (ch === '.') {
      if (!kind) kind = 'abs';
      tok += ch;
    } else if (ch in SQZ) {
      flush();
      kind = 'abs';
      tok = String(SQZ[ch]);
    } else if (ch in DIF) {
      flush();
      kind = 'dif';
      tok = String(DIF[ch]);
    } else if (ch in DUP) {
      flush();
      kind = 'dup';
      tok = String(DUP[ch]);
    } else if (ch === '+' || ch === '-') {
      flush();
      kind = 'abs';
      tok = ch === '-' ? '-' : '';
    } else if (ch === '?') {
      flush();
      vals.push(NaN);
      lastWasDif = false;
    } else {
      flush(); // espacio, coma, etc.
    }
  }
  flush();
  return { x: vals[0], ys: vals.slice(1), endsDif: lastWasDif };
}

/**
 * Decodifica un bloque (X++(Y..Y)). Con nPoints/first/delta conocidos, cada línea se coloca
 * según su valor de X (así el valor de control repetido de DIF se sobrescribe sin duplicarse).
 */
function decodeXYData(lines, { xFactor = 1, yFactor = 1, firstX, lastX, nPoints }) {
  const decoded = lines.map(decodeLine).filter((d) => Number.isFinite(d.x) && d.ys.length);
  if (!decoded.length) return null;
  const n = Number(nPoints) || 0;
  if (n > 1 && Number.isFinite(firstX) && Number.isFinite(lastX) && firstX !== lastX) {
    const delta = (lastX - firstX) / (n - 1);
    const y = new Float64Array(n).fill(NaN);
    for (const d of decoded) {
      let idx = Math.round((d.x * xFactor - firstX) / delta);
      for (const v of d.ys) {
        if (idx >= 0 && idx < n) y[idx] = v * yFactor;
        idx++;
      }
    }
    // Huecos (no deberían existir): se rellenan con el valor vecino
    for (let i = 0; i < n; i++) if (Number.isNaN(y[i])) y[i] = i ? y[i - 1] : 0;
    const x = new Float64Array(n);
    for (let i = 0; i < n; i++) x[i] = firstX + i * delta;
    return { x, y };
  }
  // Sin metadatos suficientes: concatenación secuencial descartando el valor de control DIF
  const ys = [];
  let prevDif = false;
  for (const d of decoded) {
    const vals = prevDif ? d.ys.slice(1) : d.ys;
    for (const v of vals) ys.push(v * yFactor);
    prevDif = d.endsDif;
  }
  const first = Number.isFinite(firstX) ? firstX : decoded[0].x * xFactor;
  const last = Number.isFinite(lastX) ? lastX : first + ys.length - 1;
  const delta = ys.length > 1 ? (last - first) / (ys.length - 1) : 1;
  return { x: Float64Array.from(ys, (_, i) => first + i * delta), y: Float64Array.from(ys) };
}

/** Pares (x, y[, w]) de XYPOINTS / PEAK TABLE. */
function decodePairs(lines, arity, { xFactor = 1, yFactor = 1 }) {
  const nums = lines
    .join(' ')
    .split(/[\s,;]+/)
    .filter(Boolean)
    .map(Number)
    .filter((v) => Number.isFinite(v));
  const out = [];
  for (let i = 0; i + arity - 1 < nums.length; i += arity) {
    out.push({ x: nums[i] * xFactor, y: nums[i + 1] * yFactor, w: arity > 2 ? nums[i + 2] * xFactor : null });
  }
  return out;
}

const numList = (v) => String(v || '').split(',').map((s) => s.trim());
const num = (v) => {
  const n = parseFloat(String(v ?? '').replace(/[<>]/g, ''));
  return Number.isFinite(n) ? n : null;
};

function nucleusFrom(v) {
  const s = String(v || '').toUpperCase().replace(/[\s^<>]/g, '');
  if (/^1H$|^H1$|PROTON/.test(s)) return '1H';
  if (/^13C$|^C13$|CARBON/.test(s)) return '13C';
  return s || null;
}

/** Disolvente del archivo → clave de SOLVENTS. */
export function solventKey(name) {
  const s = String(name || '').toUpperCase().replace(/[\s<>\-_]/g, '');
  if (!s) return null;
  if (/CDCL3|CHLOROFORM/.test(s)) return 'CDCl3';
  if (/DMSO/.test(s)) return 'DMSO-d6';
  if (/D2O|DEUTERIUMOXIDE|WATER/.test(s)) return 'D2O';
  if (/CD3OD|MEOD|METHANOL/.test(s)) return 'CD3OD';
  if (/ACETONE|CD3COCD3/.test(s)) return 'Acetona-d6';
  if (/C6D6|BENZENE/.test(s)) return 'C6D6';
  if (/CD3CN|ACETONITRIL/.test(s)) return 'CD3CN';
  if (/CD2CL2|DICHLOROMETHANE|METHYLENECHLORIDE/.test(s)) return 'CD2Cl2';
  return null;
}

/** Interpreta un bloque: metadatos y datos (curva o tabla de picos). */
function readBlock(block) {
  const get = (label) => block.ldr.find((l) => l.label === label)?.value ?? null;
  const meta = {
    title: get('TITLE'),
    dataType: (get('DATATYPE') || '').toUpperCase(),
    dataClass: (get('DATACLASS') || '').toUpperCase(),
    nucleus: nucleusFrom(get('.OBSERVENUCLEUS')) || nucleusFrom(get('$NUC1')),
    observeFrequency: num(get('.OBSERVEFREQUENCY')) ?? num(get('$SFO1')) ?? num(get('$BF1')),
    sf: num(get('$SF')),
    solventName: (get('.SOLVENTNAME') || get('$SOLVENT') || '').replace(/[<>]/g, '').trim() || null,
    shiftReference: get('.SHIFTREFERENCE'),
    offset: num(get('$OFFSET')),
    origin: get('ORIGIN'),
    xUnits: (get('XUNITS') || '').toUpperCase(),
  };
  const result = { meta, kind: null };
  const isFid = /FID/.test(meta.dataType) || /^(SECONDS|TIME)/.test(meta.xUnits);
  const is2D = /2D|ND NMR/.test(meta.dataType) || block.ldr.some((l) => l.label === 'PAGE' && /F1\s*=/.test(l.value));

  // NTUPLES (Bruker, MestReNova): variables con FIRST/LAST/FACTOR por columna; página de la parte real
  const ntuples = get('NTUPLES');
  if (ntuples) {
    if (/FID/i.test(ntuples)) return { ...result, kind: 'fid' };
    const symbols = numList(get('SYMBOL'));
    const names = numList(get('VARNAME')).map((s) => s.toUpperCase());
    const units = numList(get('UNITS')).map((s) => s.toUpperCase());
    const first = numList(get('FIRST')).map(num);
    const last = numList(get('LAST')).map(num);
    const factor = numList(get('FACTOR')).map(num);
    const dims = numList(get('VARDIM')).map(num);
    const xi = Math.max(0, symbols.indexOf('X'));
    let yi = names.findIndex((n) => /REAL/.test(n));
    if (yi < 0) yi = symbols.indexOf('R') >= 0 ? symbols.indexOf('R') : symbols.indexOf('Y');
    const ySym = symbols[yi] || 'R';
    if (/TIME|SECOND/.test(units[xi] || '') || /FID/.test(names[yi] || '')) return { ...result, kind: 'fid' };
    // Busca la tabla de datos de la variable real
    const tables = block.ldr.filter((l) => l.label === 'DATATABLE' || l.label === 'XYDATA');
    const table = tables.find((t) => new RegExp(`\\(\\s*${ySym}\\s*\\.\\.`, 'i').test(t.value)) || tables[0];
    if (!table) return result;
    const pageIdx = block.ldr.indexOf(table);
    // NPOINTS de la página (si existe) o VAR_DIM
    let nPoints = dims[xi] || null;
    for (let k = pageIdx - 1; k >= 0 && block.ldr[k].label !== 'PAGE'; k--) if (block.ldr[k].label === 'NPOINTS') nPoints = num(block.ldr[k].value);
    if (/XYPOINTS|PEAK/i.test(table.value) && !/\+\+/.test(table.value)) {
      const pts = decodePairs(table.lines, 2, { xFactor: factor[xi] ?? 1, yFactor: factor[yi] ?? 1 });
      return { ...result, kind: 'peaks', peaks: pts, xUnits: units[xi] };
    }
    const xy = decodeXYData(table.lines, {
      xFactor: factor[xi] ?? 1,
      yFactor: factor[yi] ?? 1,
      firstX: first[xi],
      lastX: last[xi],
      nPoints,
    });
    if (!xy) return result;
    return { ...result, kind: is2D ? '2d' : 'xy', xy, xUnits: units[xi] || meta.xUnits };
  }

  if (isFid) return { ...result, kind: 'fid' };
  const xFactor = num(get('XFACTOR')) ?? 1;
  const yFactor = num(get('YFACTOR')) ?? 1;
  const xyd = block.ldr.find((l) => l.label === 'XYDATA');
  if (xyd && /\+\+/.test(xyd.value)) {
    const xy = decodeXYData(xyd.lines, { xFactor, yFactor, firstX: num(get('FIRSTX')), lastX: num(get('LASTX')), nPoints: num(get('NPOINTS')) });
    if (xy) return { ...result, kind: is2D ? '2d' : 'xy', xy, xUnits: meta.xUnits };
  }
  const xyp = block.ldr.find((l) => l.label === 'XYPOINTS' || (l.label === 'XYDATA' && !/\+\+/.test(l.value)));
  if (xyp) {
    const pts = decodePairs(xyp.lines, 2, { xFactor, yFactor });
    if (pts.length > 1) {
      pts.sort((a, b) => a.x - b.x);
      return { ...result, kind: 'xy', xy: { x: Float64Array.from(pts, (p) => p.x), y: Float64Array.from(pts, (p) => p.y) }, xUnits: meta.xUnits };
    }
  }
  const pt = block.ldr.find((l) => l.label === 'PEAKTABLE');
  if (pt) {
    const arity = /XYW|XYM/i.test(pt.value) ? 3 : 2;
    return { ...result, kind: 'peaks', peaks: decodePairs(pt.lines, arity, { xFactor, yFactor }), xUnits: meta.xUnits };
  }
  return result;
}

/**
 * Conversión del eje X a ppm:
 * - XUNITS = PPM: directo.
 * - XUNITS = HZ con .SHIFT REFERENCE (tipo, compuesto, punto, δ) o $OFFSET (δ del primer punto,
 *   Bruker): δ = δref + (x − x_ref) / SF.
 * - XUNITS = HZ sin referencia: δ = x / frecuencia de observación.
 */
function toPpm(xs, units, meta) {
  const u = String(units || '').toUpperCase();
  if (/PPM/.test(u)) return { x: xs, referencing: 'eje en ppm' };
  const freq = meta.sf || meta.observeFrequency;
  if (!freq) throw fail('El archivo JCAMP-DX tiene el eje en Hz pero no indica la frecuencia del espectrómetro (.OBSERVE FREQUENCY).');
  let refPoint = null;
  let refPpm = null;
  let how = '';
  if (meta.shiftReference) {
    const parts = meta.shiftReference.split(',').map((s) => s.trim());
    const point = parseInt(parts[2], 10);
    const value = parseFloat(parts[3]);
    if (Number.isFinite(point) && Number.isFinite(value) && point >= 1 && point <= xs.length) {
      refPoint = point - 1;
      refPpm = value;
      how = `.SHIFT REFERENCE (${parts.slice(0, 2).join(', ')})`;
    }
  }
  if (refPoint === null && Number.isFinite(meta.offset)) {
    refPoint = 0;
    refPpm = meta.offset;
    how = '$OFFSET de Bruker';
  }
  const out = new Float64Array(xs.length);
  if (refPoint !== null) {
    const x0 = xs[refPoint];
    for (let i = 0; i < xs.length; i++) out[i] = refPpm + (xs[i] - x0) / freq;
    return { x: out, referencing: `Hz → ppm con ${how} y SF = ${freq} MHz` };
  }
  for (let i = 0; i < xs.length; i++) out[i] = xs[i] / freq;
  return { x: out, referencing: `Hz → ppm dividiendo por ${freq} MHz` };
}

/**
 * Lee un archivo JCAMP-DX y devuelve el espectro de 1D (curva en ppm o tabla de picos).
 * `nucleus` ('1H' | '13C') elige el bloque cuando el archivo trae varios.
 */
export function parseJcamp(text, { nucleus = null } = {}) {
  const blocks = splitBlocks(text).map(readBlock);
  const usable = blocks.filter((b) => b.kind === 'xy' || b.kind === 'peaks');
  if (!usable.length) {
    if (blocks.some((b) => b.kind === 'fid'))
      throw fail('El archivo JCAMP-DX contiene sólo el FID (dominio del tiempo). Procese el espectro (transformada de Fourier, fase y referencia) y expórtelo como espectro.');
    if (blocks.some((b) => b.kind === '2d')) throw fail('El archivo es un espectro de 2D; sólo se admiten espectros de ¹H o ¹³C de 1D.');
    throw fail('No se encontraron datos de espectro en el archivo JCAMP-DX (XYDATA, NTUPLES, XYPOINTS o PEAK TABLE).');
  }
  const wanted = nucleus === '1H' || nucleus === '13C' ? nucleus : null;
  const pick =
    usable.find((b) => b.kind === 'xy' && (!wanted || b.meta.nucleus === wanted)) ||
    usable.find((b) => !wanted || b.meta.nucleus === wanted) ||
    usable[0];
  const meta = pick.meta;
  const otherNuclei = [...new Set(usable.map((b) => b.meta.nucleus).filter((n) => n && n !== meta.nucleus))];
  const base = {
    title: meta.title,
    nucleus: meta.nucleus === '1H' || meta.nucleus === '13C' ? meta.nucleus : null,
    rawNucleus: meta.nucleus,
    frequency: meta.observeFrequency ? Math.round(meta.observeFrequency * 100) / 100 : null,
    solvent: solventKey(meta.solventName),
    solventName: meta.solventName,
    origin: meta.origin,
    blocks: usable.length,
    otherNuclei,
  };
  if (base.rawNucleus && !base.nucleus) throw fail(`El archivo JCAMP-DX es de ${base.rawNucleus}; sólo se admiten espectros de ¹H y ¹³C.`);
  if (pick.kind === 'peaks') {
    const conv = toPpm(Float64Array.from(pick.peaks, (p) => p.x), pick.xUnits, meta);
    const freq = meta.sf || meta.observeFrequency;
    const peaks = pick.peaks.map((p, i) => ({
      ppm: conv.x[i],
      height: Math.abs(p.y),
      width: p.w !== null && Number.isFinite(p.w) ? (/PPM/.test(pick.xUnits || '') && freq ? p.w * freq : p.w) : null,
    }));
    return { ...base, kind: 'peaks', peaks, referencing: conv.referencing };
  }
  const conv = toPpm(pick.xy.x, pick.xUnits, meta);
  // Eje ascendente en ppm
  let x = conv.x;
  let y = pick.xy.y;
  if (x.length > 1 && x[0] > x[x.length - 1]) {
    x = Float64Array.from(x).reverse();
    y = Float64Array.from(y).reverse();
  }
  return { ...base, kind: 'xy', x, y, referencing: conv.referencing };
}
