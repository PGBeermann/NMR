import path from 'node:path';
import { fileURLToPath } from 'node:url';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');

export const config = {
  root,
  port: Number(process.env.PORT) || 3000,
  dataDir: path.resolve(process.env.DATA_DIR || path.join(root, 'data')),
  // NMRShiftDB2 (base de datos abierta de espectros asignados, ~50 000 compuestos)
  nmrshiftdbUrl:
    process.env.NMRSHIFTDB_URL ||
    'https://sourceforge.net/projects/nmrshiftdb2/files/data/nmrshiftdb2withsignals.sd/download',
  // Si es "0", no se descarga ni indexa NMRShiftDB2 automáticamente al arrancar
  nmrshiftdbAuto: process.env.NMRSHIFTDB_AUTO !== '0',
  pubchemBase: 'https://pubchem.ncbi.nlm.nih.gov/rest/pug',
  pubchemViewBase: 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view',
  opsinBase: 'https://www.ebi.ac.uk/opsin/ws',
  httpTimeoutMs: Number(process.env.HTTP_TIMEOUT_MS) || 25000,
  maxFormulaResults: Number(process.env.MAX_FORMULA_RESULTS) || 500,
};
