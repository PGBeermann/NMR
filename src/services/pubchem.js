// Consultas a PubChem (PUG REST y PUG View).
// PubChem agrega espectros experimentales de RMN de varias fuentes (HMDB,
// SpectraBase/Wiley, NMRShiftDB, catálogo Varian…).
import { config } from '../config.js';
import { fetchJson } from './http.js';

const B = config.pubchemBase;
const PROPS = 'Title,IUPACName,MolecularFormula,MolecularWeight,ExactMass,SMILES,ConnectivitySMILES,InChIKey';

export async function cidsByName(name) {
  const d = await fetchJson(`${B}/compound/name/${encodeURIComponent(name)}/cids/JSON`);
  return (d?.IdentifierList?.CID || []).filter(Boolean);
}

export async function cidsBySmiles(smiles) {
  const d = await fetchJson(`${B}/compound/smiles/cids/JSON`, {
    method: 'POST',
    body: new URLSearchParams({ smiles }),
    headers: { 'Content-Type': 'application/x-www-form-urlencoded' },
    cacheKey: `smiles:${smiles}`,
  });
  return (d?.IdentifierList?.CID || []).filter(Boolean);
}

export async function cidsByInchikey(inchikey) {
  const d = await fetchJson(`${B}/compound/inchikey/${encodeURIComponent(inchikey)}/cids/JSON`);
  return (d?.IdentifierList?.CID || []).filter(Boolean);
}

export async function cidsByFormula(formula, max = config.maxFormulaResults) {
  const d = await fetchJson(
    `${B}/compound/fastformula/${encodeURIComponent(formula)}/cids/JSON?MaxRecords=${max}`,
  );
  return (d?.IdentifierList?.CID || []).filter(Boolean);
}

export async function properties(cids) {
  const out = [];
  for (let i = 0; i < cids.length; i += 100) {
    const chunk = cids.slice(i, i + 100);
    const d = await fetchJson(`${B}/compound/cid/${chunk.join(',')}/property/${PROPS}/JSON`);
    for (const p of d?.PropertyTable?.Properties || []) {
      out.push({
        cid: p.CID,
        title: p.Title,
        iupac: p.IUPACName,
        formula: p.MolecularFormula,
        mw: Number(p.MolecularWeight),
        exactMass: Number(p.ExactMass),
        smiles: p.SMILES || p.IsomericSMILES || p.CanonicalSMILES || p.ConnectivitySMILES,
        inchikey: p.InChIKey,
      });
    }
  }
  return out;
}

export async function propertiesByInchikeys(keys) {
  const out = [];
  for (let i = 0; i < keys.length; i += 50) {
    const chunk = keys.slice(i, i + 50);
    const d = await fetchJson(`${B}/compound/inchikey/property/${PROPS}/JSON`, {
      method: 'POST',
      body: new URLSearchParams({ inchikey: chunk.join(',') }),
      headers: { 'Content-Type': 'application/x-www-form-urlencoded' },
      cacheKey: `ik:${chunk.join(',')}`,
    });
    for (const p of d?.PropertyTable?.Properties || []) {
      out.push({ cid: p.CID, title: p.Title, iupac: p.IUPACName, formula: p.MolecularFormula, inchikey: p.InChIKey });
    }
  }
  return out;
}

export async function synonyms(cid, max = 12) {
  const d = await fetchJson(`${B}/compound/cid/${cid}/synonyms/JSON`);
  return (d?.InformationList?.Information?.[0]?.Synonym || []).slice(0, max);
}

function valueToText(v) {
  if (!v) return '';
  if (v.StringWithMarkup) return v.StringWithMarkup.map((s) => s.String).join(' ');
  if (v.Number) return v.Number.join(' ') + (v.Unit ? ` ${v.Unit}` : '');
  if (v.ExternalDataURL) return v.ExternalDataURL[0];
  return '';
}

function parsePeakString(str) {
  // "1.26:560.00, 1.18:214.00, ..." → [{ppm, intensity}]
  const peaks = [];
  for (const part of str.split(/[,;]\s*/)) {
    const m = part.trim().match(/^(-?\d+(?:\.\d+)?)\s*:\s*(-?\d+(?:\.\d+)?)/);
    if (m) peaks.push({ ppm: Number(m[1]), intensity: Number(m[2]) });
  }
  return peaks;
}

/**
 * Devuelve los espectros de RMN 1D publicados en PubChem para un CID.
 * Cada espectro agrupa la información por referencia (fuente).
 */
export async function nmrSpectra(cid) {
  const d = await fetchJson(
    `${config.pubchemViewBase}/data/compound/${cid}/JSON?heading=${encodeURIComponent('1D NMR Spectra')}`,
  );
  if (!d?.Record) return { spectra: [], links: [] };
  const refs = new Map((d.Record.Reference || []).map((r) => [r.ReferenceNumber, r]));
  const spectra = [];
  const links = [];

  function walk(sections) {
    for (const s of sections || []) {
      const heading = s.TOCHeading || '';
      const nucleus = /^1H\b/.test(heading) ? '1H' : /^13C\b/.test(heading) ? '13C' : null;
      if (s.Information) {
        if (!nucleus) {
          for (const info of s.Information) {
            const ref = refs.get(info.ReferenceNumber);
            const url = info.Value?.StringWithMarkup?.[0]?.Markup?.[0]?.URL;
            links.push({
              text: valueToText(info.Value),
              url: url || ref?.URL,
              source: ref?.SourceName,
            });
          }
        } else {
          const byRef = new Map();
          for (const info of s.Information) {
            const key = info.ReferenceNumber;
            if (!byRef.has(key)) byRef.set(key, []);
            byRef.get(key).push(info);
          }
          for (const [refNo, infos] of byRef) {
            const ref = refs.get(refNo) || {};
            const fields = {};
            let peaks = [];
            let thumbnail = null;
            for (const info of infos) {
              const name = info.Name || 'Dato';
              if (/Shifts/i.test(name)) {
                peaks = parsePeakString(valueToText(info.Value));
                fields[name] = valueToText(info.Value);
              } else if (/Thumbnail/i.test(name)) {
                thumbnail = info.Value?.ExternalDataURL?.[0] || null;
              } else {
                fields[name] = valueToText(info.Value);
              }
            }
            const freq = Number(String(fields.Frequency || '').match(/[\d.]+/)?.[0]) || 0;
            spectra.push({
              nucleus,
              source: ref.SourceName || 'PubChem',
              sourceId: ref.SourceID || String(refNo),
              url: ref.URL || null,
              license: ref.LicenseURL || ref.LicenseNote || null,
              frequency: freq,
              solvent: fields.Solvent || null,
              fields,
              peaks,
              thumbnail,
            });
          }
        }
      }
      walk(s.Section);
    }
  }
  walk(d.Record.Section);
  return { spectra, links };
}

export const pubchemUrl = (cid) => `https://pubchem.ncbi.nlm.nih.gov/compound/${cid}`;
