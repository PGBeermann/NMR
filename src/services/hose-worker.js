import { parentPort, workerData } from 'node:worker_threads';
import { buildHoseTable } from './hose.js';

try {
  const r = buildHoseTable({ dataDir: workerData.dataDir, onProgress: (p) => parentPort.postMessage({ type: 'progress', ...p }) });
  parentPort.postMessage({ type: 'done', ...r });
} catch (e) {
  parentPort.postMessage({ type: 'error', message: e.message });
}
