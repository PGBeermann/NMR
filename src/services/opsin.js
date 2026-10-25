// OPSIN (EBI): convierte nombres IUPAC sistemáticos en estructuras.
// Se usa como respaldo cuando PubChem no reconoce el nombre.
import { config } from '../config.js';
import { fetchJson } from './http.js';

export async function nameToStructure(name) {
  const d = await fetchJson(`${config.opsinBase}/${encodeURIComponent(name)}.json`);
  if (!d || d.status !== 'SUCCESS' || !d.smiles) return null;
  return { smiles: d.smiles, inchikey: d.stdinchikey || null, inchi: d.stdinchi || null };
}
