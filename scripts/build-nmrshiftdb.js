// Uso: npm run index-nmrshiftdb   (descarga NMRShiftDB2 y construye data/nmrshiftdb.db)
import { config } from '../src/config.js';
import { buildIndex } from '../src/services/nmrshiftdb-indexer.js';

const t0 = Date.now();
const r = await buildIndex({
  url: config.nmrshiftdbUrl,
  dataDir: config.dataDir,
  onProgress: (p) => console.log(JSON.stringify(p)),
});
console.log(`Listo: ${r.records} compuestos (${r.errors} errores) en ${((Date.now() - t0) / 1000).toFixed(0)} s → ${r.file}`);
