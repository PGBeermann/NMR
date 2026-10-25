import { parentPort, workerData } from 'node:worker_threads';
import { buildIndex } from './nmrshiftdb-indexer.js';

buildIndex({ ...workerData, onProgress: (p) => parentPort.postMessage({ type: 'progress', ...p }) })
  .then((r) => parentPort.postMessage({ type: 'done', ...r }))
  .catch((e) => parentPort.postMessage({ type: 'error', message: e.message }));
