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

console.log('Todas las pruebas pasaron ✔');
