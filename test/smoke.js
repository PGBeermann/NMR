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

console.log('Todas las pruebas pasaron ✔');
