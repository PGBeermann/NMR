import express from 'express';
import path from 'node:path';
import { config } from './config.js';
import { api } from './routes/api.js';
import * as nmrshiftdb from './services/nmrshiftdb.js';

const app = express();
app.disable('x-powered-by');
app.use(express.json({ limit: '2mb' }));

app.get('/health', (req, res) => res.json({ ok: true }));
app.use('/api', api);
app.use(express.static(path.join(config.root, 'public'), { maxAge: '1h' }));

app.use('/api', (req, res) => res.status(404).json({ error: 'Ruta no encontrada' }));
// eslint-disable-next-line no-unused-vars
app.use((err, req, res, next) => {
  const status = err.status || 500;
  if (status >= 500) console.error(err);
  res.status(status).json({ error: err.message || 'Error interno' });
});

app.listen(config.port, () => {
  console.log(`Analizador RMN escuchando en http://localhost:${config.port}  (datos: ${config.dataDir})`);
  if (nmrshiftdb.isReady()) {
    console.log(`NMRShiftDB2 lista: ${nmrshiftdb.status.records} compuestos`);
    nmrshiftdb.ensureHose();
  } else if (config.nmrshiftdbAuto) {
    console.log('NMRShiftDB2 no encontrada: se descargará e indexará en segundo plano.');
    nmrshiftdb.startIndexing();
  }
});
