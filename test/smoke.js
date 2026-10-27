// Pruebas rápidas sin conexión: npm test
import assert from 'node:assert/strict';
import { normalizeFormula, unsaturation, parseMolecule, moleculeInfo } from '../src/services/molecule.js';
import { describeCarbon, describeHydrogenCarrier, multiplicityName } from '../src/services/chemistry.js';
import { spectrumSvg, svgToPng } from '../src/services/plot.js';
import { sticksToXY } from '../src/services/spectrum.js';

assert.equal(normalizeFormula('H8C4O2'), 'C4H8O2');
assert.equal(normalizeFormula('C4 H8 O2'), 'C4H8O2');
assert.equal(unsaturation('C6H6'), 4);
assert.equal(unsaturation('C4H8O2'), 1);
assert.equal(multiplicityName('dd'), 'doblete de dobletes');
assert.equal(multiplicityName('q'), 'cuarteto');

const mol = parseMolecule({ smiles: 'CCOC(C)=O' });
const info = moleculeInfo(mol);
assert.equal(info.formula, 'C4H8O2');
assert.equal(describeCarbon(mol, 3).code, 'ester');
assert.equal(describeCarbon(mol, 1).code, 'C-O');
assert.equal(describeHydrogenCarrier(parseMolecule({ smiles: 'CC(=O)O' }), 3).code, 'H-acido');

const sticks = [{ ppm: 170.7, intensity: 45 }, { ppm: 60.4, intensity: 100 }];
const xy = sticksToXY(sticks, { from: -5, to: 220, width: 0.1 });
const svg = spectrumSvg({ from: -5, to: 220, xy, title: 'Prueba', labels: sticks.map((s) => ({ ppm: s.ppm, text: String(s.ppm), height: s.intensity })) });
assert.ok(svg.startsWith('<svg'));
assert.ok(svgToPng(svg, 600).length > 1000);

// Modo 4: lista de picos → multipletes, integrales e interpretación
import fs from 'node:fs';
import { parsePeakCsv, analyzePeaks, interpret } from '../src/services/peaklist.js';
import { concordance, neighborsFromMultiplicity } from '../src/services/candidates.js';

const peaks = parsePeakCsv(fs.readFileSync(new URL('./fixtures/fenol-cdcl3-400MHz.tsv', import.meta.url), 'utf8'));
assert.equal(peaks.frequency, 400);
const an = analyzePeaks(peaks, { nucleus: '1H', frequency: 400, totalH: 6 });
assert.deepEqual(an.signals.map((s) => s.multiplicity), ['t', 't', 'd', 's']);
assert.deepEqual(an.signals.map((s) => s.nH), [2, 1, 2, 1]);
assert.ok(interpret(an).some((x) => /monosustituido/.test(x.fragment)));
assert.equal(parsePeakCsv('ppm;Int\n7,26;10\n1,25;5').length, 2);
assert.equal(neighborsFromMultiplicity('dd'), 2);
assert.equal(neighborsFromMultiplicity('m'), null);
const exp = an.signals.map((s) => ({ delta: s.delta, w: s.fraction, nb: neighborsFromMultiplicity(s.multiplicity) }));
const good = concordance(exp, [{ delta: 7.24, n: 2, nb: 2 }, { delta: 6.93, n: 1, nb: 2 }, { delta: 6.84, n: 2, nb: 1 }, { delta: 5.3, n: 1, nb: 0 }], '1H');
const bad = concordance(exp, [{ delta: 3.7, n: 3, nb: 0 }, { delta: 1.2, n: 3, nb: 2 }], '1H');
assert.ok(good.score > 90 && bad.score < 30, `concordancia ${good.score} / ${bad.score}`);
// Criterio de DE de los errores ≤ 0.1 ppm: un desplazamiento sistemático común no se penaliza en la DE…
const shifted = concordance(exp, [{ delta: 7.32, n: 2, nb: 2 }, { delta: 7.01, n: 1, nb: 2 }, { delta: 6.92, n: 2, nb: 1 }, { delta: 5.43, n: 1, nb: 0 }], '1H');
assert.ok(shifted.sdOk && shifted.errorSd < 0.02, `DE ${shifted.errorSd}`);
// …pero errores dispersos (+0.2 / −0.2) sí superan el límite
const scattered = concordance(exp, [{ delta: 7.44, n: 2, nb: 2 }, { delta: 6.73, n: 1, nb: 2 }, { delta: 7.04, n: 2, nb: 1 }, { delta: 5.15, n: 1, nb: 0 }], '1H');
assert.ok(!scattered.sdOk && scattered.errorSd > 0.1, `DE ${scattered.errorSd}`);

// o-Cresol a 90 MHz: los 4 H aromáticos forman un solo multiplete; las señales del candidato
// se agrupan y la forma del multiplete favorece al isómero orto frente al para
const pc = parsePeakCsv(fs.readFileSync(new URL('./fixtures/o-cresol-90MHz.tsv', import.meta.url), 'utf8'));
assert.equal(pc.frequency, 90);
const ac = analyzePeaks(pc, { nucleus: '1H', frequency: 90, totalH: 8 });
assert.equal(ac.signals.length, 3);
const expC = ac.signals.map((s) => ({ delta: s.delta, w: s.fraction, mult: s.multiplicity, nb: neighborsFromMultiplicity(s.multiplicity), from: s.from, to: s.to, lines: s.lines, freq: 90 }));
const ortho = concordance(expC, [{ delta: 2.2, n: 3, nb: 0 }, { delta: 7.12, n: 1, nb: 1 }, { delta: 7.08, n: 1, nb: 2 }, { delta: 6.86, n: 1, nb: 2 }, { delta: 6.77, n: 1, nb: 1 }], '1H');
const para = concordance(expC, [{ delta: 2.27, n: 3, nb: 0 }, { delta: 7.03, n: 2, nb: 1 }, { delta: 6.74, n: 2, nb: 1 }], '1H');
assert.ok(ortho.sdOk, `DE orto ${ortho.errorSd}`);
assert.ok(ortho.score > para.score, `orto ${ortho.score} vs para ${para.score}`);

// JCAMP-DX: decodificación ASDF, conversión Hz → ppm, detección de picos y áreas
const { decodeLine, parseJcamp, isJcamp } = await import('../src/services/jcamp.js');
const { pickPeaks } = await import('../src/services/peakpick.js');
// Ejemplos de la norma: los tres codifican X = 1000 e Y = 1000 2000 2001 2002 2003 2003 2003
assert.deepEqual(decodeLine('1000 1000 2000 2001 2002 2003 2003 2003').ys, [1000, 2000, 2001, 2002, 2003, 2003, 2003]); // AFFN
assert.deepEqual(decodeLine('1000A000B000B001B002B003B003B003').ys, [1000, 2000, 2001, 2002, 2003, 2003, 2003]); // SQZ
assert.deepEqual(decodeLine('1000A000J000JJJ%%').ys, [1000, 2000, 2001, 2002, 2003, 2003, 2003]); // DIF
assert.deepEqual(decodeLine('1000A000J000JU%T').ys, [1000, 2000, 2001, 2002, 2003, 2003, 2003]); // DIFDUP
assert.deepEqual(decodeLine('1000+1000+2000-5').ys, [1000, 2000, -5]); // PAC
assert.ok(isJcamp('##TITLE=x\n##JCAMP-DX=5.01'));
assert.throws(() => parseJcamp('##TITLE=fid\n##JCAMP-DX=6.0\n##DATA TYPE=NMR FID\n##.OBSERVE NUCLEUS=^1H\n##END='), /FID/);
const pt = parseJcamp('##TITLE=picos\n##JCAMP-DX=5.01\n##DATA TYPE=NMR PEAK TABLE\n##.OBSERVE NUCLEUS=^13C\n##.OBSERVE FREQUENCY=100\n##XUNITS=HZ\n##PEAK TABLE=(XY..XY)\n17000,50 6050,100\n##END=');
assert.equal(pt.kind, 'peaks');
assert.equal(pt.nucleus, '13C');
assert.ok(Math.abs(pt.peaks[0].ppm - 170) < 1e-9 && Math.abs(pt.peaks[1].ppm - 60.5) < 1e-9);

const jh = parseJcamp(fs.readFileSync(new URL('../public/examples/acetato-de-etilo-1H-400MHz.jdx', import.meta.url), 'latin1'));
assert.equal(jh.nucleus, '1H');
assert.equal(jh.frequency, 400.13);
assert.equal(jh.solvent, 'CDCl3');
assert.equal(jh.x.length, 32768);
assert.ok(Math.abs(jh.x[jh.x.length - 1] - 10.5) < 1e-6 && Math.abs(jh.x[0] + 0.5) < 1e-6, 'eje en ppm con .SHIFT REFERENCE');
const pkH = pickPeaks(jh.x, jh.y, { nucleus: '1H', frequency: jh.frequency, solvent: 'CDCl3' });
const maxH = Math.max(...pkH.peaks.map((p) => p.height));
const peaksH = pkH.peaks.map((p) => ({ ppm: p.ppm, intensity: (p.height / maxH) * 100, area: p.area }));
const aj = analyzePeaks(peaksH, { nucleus: '1H', frequency: jh.frequency, solvent: 'CDCl3', excludeSolvent: true, totalH: 8 });
assert.deepEqual(
  aj.signals.map((s) => [s.delta.toFixed(2), s.multiplicity, s.nH]),
  [
    ['4.12', 'q', 2],
    ['2.05', 's', 3],
    ['1.26', 't', 3],
  ],
);
assert.ok(Math.abs(aj.signals[0].J[0] - 7.1) < 0.15, 'J del cuarteto');

const jc = parseJcamp(fs.readFileSync(new URL('../public/examples/acetato-de-etilo-13C-100MHz.jdx', import.meta.url), 'latin1'));
assert.equal(jc.nucleus, '13C');
const pkC = pickPeaks(jc.x, jc.y, { nucleus: '13C', frequency: jc.frequency, solvent: 'CDCl3' });
const maxC = Math.max(...pkC.peaks.map((p) => p.height));
const acj = analyzePeaks(pkC.peaks.map((p) => ({ ppm: p.ppm, intensity: (p.height / maxC) * 100, area: p.area })), {
  nucleus: '13C',
  frequency: jc.frequency,
  solvent: 'CDCl3',
});
assert.deepEqual(acj.signals.filter((s) => !s.solvent).map((s) => Math.round(s.delta * 10) / 10), [171.1, 60.5, 21, 14.2]);
assert.equal(acj.signals.filter((s) => s.solvent).length, 3, 'triplete del CDCl3');

// Predictor local de respaldo (sólo si la tabla de entornos ya está construida en data/)
const { hoseReady } = await import('../src/services/hose.js');
if (hoseReady()) {
  const { localPrediction } = await import('../src/services/local-prediction.js');
  const lp = localPrediction(parseMolecule({ smiles: 'CCOC(C)=O' }));
  const byMult = Object.fromEntries(lp.spectra.proton.joinedSignals.map((s) => [s.multiplicity, s]));
  assert.ok(byMult.q && Math.abs(byMult.q.delta - 4.1) < 0.2 && byMult.q.nbAtoms === 2, 'OCH2 cuarteto');
  assert.ok(byMult.t && Math.abs(byMult.t.delta - 1.25) < 0.2 && byMult.t.nbAtoms === 3, 'CH3 triplete');
  assert.equal(lp.spectra.carbon.joinedSignals.length, 4);
} else {
  console.log('(predictor local no construido: se omite su prueba)');
}

console.log('Todas las pruebas pasaron ✔');
