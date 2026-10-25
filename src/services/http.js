import { config } from '../config.js';

// Caché sencilla en memoria para no repetir consultas a servicios externos
const cache = new Map();
const MAX_CACHE = 500;

export class HttpError extends Error {
  constructor(status, message) {
    super(message);
    this.status = status;
  }
}

export async function fetchJson(url, { method = 'GET', body, headers = {}, retries = 2, cacheKey } = {}) {
  const key = cacheKey ?? (method === 'GET' ? url : null);
  if (key && cache.has(key)) return cache.get(key);

  for (let attempt = 0; ; attempt++) {
    const ctrl = new AbortController();
    const timer = setTimeout(() => ctrl.abort(), config.httpTimeoutMs);
    try {
      const res = await fetch(url, {
        method,
        body,
        headers: { Accept: 'application/json', 'User-Agent': 'nmr-analyzer/1.0', ...headers },
        signal: ctrl.signal,
        redirect: 'follow',
      });
      if (res.status === 404) return null;
      if ((res.status === 503 || res.status === 429) && attempt < retries) {
        await new Promise((r) => setTimeout(r, 800 * (attempt + 1)));
        continue;
      }
      if (!res.ok) {
        const text = await res.text().catch(() => '');
        // PubChem devuelve 400 cuando no encuentra un nombre/fórmula
        if (res.status === 400 && /not found|No CID|Invalid/i.test(text)) return null;
        throw new HttpError(res.status, `Error ${res.status} al consultar ${new URL(url).host}`);
      }
      const data = await res.json();
      if (key) {
        if (cache.size >= MAX_CACHE) cache.delete(cache.keys().next().value);
        cache.set(key, data);
      }
      return data;
    } catch (err) {
      if (err.name === 'AbortError') {
        if (attempt < retries) continue;
        throw new HttpError(504, `Tiempo de espera agotado consultando ${new URL(url).host}`);
      }
      if (err instanceof HttpError) throw err;
      if (attempt < retries) {
        await new Promise((r) => setTimeout(r, 500));
        continue;
      }
      throw new HttpError(502, `No se pudo contactar ${new URL(url).host}: ${err.message}`);
    } finally {
      clearTimeout(timer);
    }
  }
}
