# Analizador RMN ¹H / ¹³C

Aplicación web en **Node.js** para el análisis de espectros de **RMN-¹H** y **RMN-¹³C** de compuestos orgánicos.

## Modos de análisis

1. **Dibujo de la estructura (JSME)**: se dibuja la molécula en el editor JSME y se obtiene:
   - Espectros de ¹H y ¹³C **predichos** y **simulados** como imagen (SVG/PNG), con etiquetas, integrales y curva de integración.
   - Tabla de señales: δ (ppm), multiplicidad, **constantes de acoplamiento J (Hz)** con el H acoplado y el nº de enlaces, integración, átomo asignado, entorno químico, rango típico y explicación (regla n+1).
   - **Asignación con letras**: cada señal de ¹H recibe una letra (a, b, c…) que aparece en la estructura (junto al átomo que porta esos H), sobre el pico del espectro y en la tabla.
   - ¹³C: tipo de carbono (CH₃/CH₂/CH/C), **DEPT-135 y DEPT-90**, intensidades relativas, interpretación por regiones.
   - Lista de líneas del espectro simulado con **intensidades** relativas (ppm y Hz).
   - Resumen: fórmula, masa, índice de insaturación, grupos funcionales, simetría, formato de publicación.
   - **Búsqueda de espectros experimentales** de la misma estructura en PubChem (HMDB, SpectraBase…) y NMRShiftDB2, con **comparación predicho vs. experimental** (error medio por átomo). Los espectros experimentales no se muestran automáticamente: hay botones separados para ¹H y ¹³C.
2. **Fórmula molecular**: lista **todos** los compuestos con esa fórmula en PubChem (nombre común, nombre IUPAC, SMILES, imagen) y marca los que tienen espectros asignados en NMRShiftDB2. Cualquiera se puede analizar con un clic.
3. **Nombre IUPAC** (o común): se resuelve en PubChem; si no existe, el nombre sistemático se convierte en estructura con OPSIN. Luego se buscan y predicen los espectros.

4. **Espectro experimental (CSV) → estructura**: se carga o pega una tabla de picos (`ppm` e intensidad; opcional una columna `Hz`, con la que se calcula la frecuencia del equipo). El programa:
   - agrupa las líneas en multipletes, calcula su centro, la multiplicidad (fusionando la estructura fina < 2.5 Hz), las constantes **J** y la integración relativa **por área** (con la fórmula reparte el nº total de H; el usuario puede corregir los nH y recalcular). El área de cada línea se toma de una columna `Área`/`Integral` del archivo, o se calcula como altura × anchura si hay columna `W½`/`Width`; si sólo hay alturas, se estima con líneas lorentzianas de anchura común y deconvolución de las líneas solapadas;
   - marca las señales del disolvente y propone **fragmentos estructurales** (benceno monosustituido, sustituyente donador/aceptor, OCH₃, acetilo, etilo, isopropilo, OH/NH…);
   - busca **estructuras candidatas** en los ~17 600 espectros de ¹H (o ~58 000 de ¹³C) de NMRShiftDB2 y, si se indica la fórmula, predice los espectros de todos los isómeros de PubChem con esa fórmula (hasta 150) con un **predictor local** (códigos de entorno tipo HOSE construidos a partir de NMRShiftDB2, `data/hose.db`, ~2 min la primera vez) y refina los 3 mejores con el servicio NMRium (que limita el nº de consultas);
   - identifica cada multiplete con una letra y dibuja cada estructura candidata con esas letras sobre los átomos cuyos H se asignan a cada pico;
   - ordena los candidatos por **% de concordancia** = desplazamientos + cobertura de las señales del candidato + integración + nº de señales + **multiplicidad esperada por la regla n+1** (calculada desde la estructura), con reconocimiento de OH/NH intercambiables. Para los desplazamientos se calcula el error δ(candidato) − δ(experimental) de cada señal asignada y **se aceptan sólo los candidatos cuya desviación estándar de esos errores sea ≤ 0.1 ppm** (¹H; 1.5 ppm recomendado para ¹³C). La DE no castiga un desplazamiento sistemático común a todas las señales (p. ej. por el disolvente); el error medio se muestra y se penaliza de forma suave. Además **sólo se muestran los candidatos con concordancia ≥ 90 %** (ambos valores se pueden cambiar en la pantalla o con `sdMax` y `minScore` en la API).
   Ejemplo incluido: `public/examples/fenol-cdcl3-400MHz.tsv`.

En **todos los casos** el análisis se **guarda automáticamente** en una base de datos SQLite (`data/nmr.db`): compuesto, señales, desplazamientos, intensidades, constantes de acoplamiento, integración, multiplicidad, asignación, disolvente, frecuencia, fuente y comparación. La pestaña **Espectros guardados** permite consultarlos, exportarlos y eliminarlos.

### Formatos de exportación
PNG · SVG · **JCAMP-DX 5.01** (curva x/y y tabla de picos, con señales, J y multiplicidades) · CSV (señales + líneas con intensidades) · CSV x/y · JSON (espectro o análisis completo).

## Fuentes de datos
| Fuente | Uso |
|---|---|
| [nmr-processing](https://www.npmjs.com/package/nmr-processing) (NMRium / Zakodium) | Predicción de δ y J; simulación del sistema de espín (efectos de 2.º orden) |
| [PubChem](https://pubchem.ncbi.nlm.nih.gov) PUG REST / PUG View | Nombres, fórmula, búsqueda por fórmula/nombre/estructura, espectros experimentales (HMDB con intensidades, miniaturas de SpectraBase) |
| [NMRShiftDB2](https://nmrshiftdb.nmr.uni-koeln.de) | ~58 000 compuestos con espectros asignados. Se descarga e indexa **localmente** la primera vez que arranca el servidor (≈ 3 min, ~230 MB en el volumen) |
| [OPSIN](https://www.ebi.ac.uk/opsin/) | Nombre IUPAC → estructura cuando PubChem no conoce el nombre |
| [JSME](https://jsme-editor.github.io) | Editor de estructuras en el navegador |
| [OpenChemLib](https://github.com/cheminfo/openchemlib-js) | Estructuras, fórmula, simetría (átomos diastereotópicos), dibujo numerado |

> La predicción usa el servicio web público de NMRium; el servidor necesita salida a Internet. Si NMRium no responde (p. ej. porque limita el nº de consultas), el análisis usa automáticamente el **predictor local** de respaldo: desplazamientos por entornos tipo HOSE de NMRShiftDB2 y constantes J típicas (regla n+1); el espectro queda marcado como «Predicción local (HOSE/NMRShiftDB2, respaldo)». `nmr-processing` tiene licencia CC-BY-NC-SA 4.0 (uso académico/no comercial).

## Ejecución local
Requiere Node.js ≥ 22.13 (recomendado 24).
```bash
npm install
npm start          # http://localhost:3000
npm test           # pruebas rápidas
npm run index-nmrshiftdb   # (opcional) construir el índice de NMRShiftDB2 manualmente
```

## Despliegue en Dokploy (VPS Hostinger)
El repositorio incluye `Dockerfile` (y `docker-compose.yml` como alternativa).

1. En Dokploy: **Create Project** → **Create Service → Application**.
2. **Provider**: GitHub (o *Git* con la URL `https://github.com/PGBeermann/NMR.git`), rama `main`.
3. **Build Type**: `Dockerfile` (ruta `./Dockerfile`).
4. **Environment** (opcional):
   ```
   PORT=3000
   DATA_DIR=/app/data
   NMRSHIFTDB_AUTO=1
   ```
5. **Advanced → Volumes/Mounts**: *Volume mount* con nombre `nmr-data` y ruta de montaje `/app/data` (conserva los espectros guardados y el índice NMRShiftDB2 entre despliegues).
6. **Domains**: añadir el dominio (o el generado por Dokploy), **Container Port `3000`**, HTTPS con Let's Encrypt.
7. **Deploy**. En el primer arranque, la barra superior muestra el progreso de la descarga/indexación de NMRShiftDB2; mientras tanto, la predicción y PubChem ya funcionan.

### Variables de entorno
| Variable | Defecto | Descripción |
|---|---|---|
| `PORT` | 3000 | Puerto HTTP |
| `DATA_DIR` | `./data` | Carpeta de datos persistentes |
| `NMRSHIFTDB_AUTO` | 1 | `0` para no descargar NMRShiftDB2 |
| `NMRSHIFTDB_URL` | SourceForge | URL del archivo `nmrshiftdb2withsignals.sd` |
| `MAX_FORMULA_RESULTS` | 500 | Máximo de compuestos de PubChem por fórmula |
| `HTTP_TIMEOUT_MS` | 25000 | Tiempo máximo por consulta externa |

## API
| Método | Ruta | Descripción |
|---|---|---|
| POST | `/api/analyze` | `{smiles|molfile|cid, name, frequency, solvent}` → análisis completo (se guarda) |
| GET | `/api/search/formula?formula=C4H8O2` | Compuestos con esa fórmula |
| GET | `/api/search/name?q=...` | Resolución de nombre |
| GET | `/api/compounds` · `/api/compounds/:id` | Guardados |
| GET | `/api/spectra/:id/image.svg|png?from=&to=` | Imagen del espectro (con zoom) |
| GET | `/api/spectra/:id/export?format=jdx|jdx-peaks|csv|xy|json` | Exportación |
| POST | `/api/spectrum/analyze` | `{csv, name, nucleus, frequency, solvent, formula, excludeSolvent, rangeOnly, integrals, sdMax=0.1, minScore=90}` → multipletes, interpretación y candidatos (se guarda) |
| GET | `/api/uploads` · `/api/uploads/:id` | Espectros CSV cargados |
| GET | `/api/uploads/:id/image.svg|png` · `/api/uploads/:id/export?format=csv|json` | Imagen y exportación |
| GET | `/api/status` | Estado de las bases de datos |

## Estructura
```
src/server.js              servidor Express
src/routes/api.js          rutas de la API
src/services/analysis.js   flujo de análisis, interpretación y comparación
src/services/chemistry.js  entornos químicos, rangos, multiplicidades, DEPT, disolventes
src/services/pubchem.js    PubChem (PUG REST / PUG View)
src/services/nmrshiftdb*.js índice local de NMRShiftDB2
src/services/plot.js       imágenes SVG/PNG de espectros
src/services/exporters.js  JCAMP-DX, CSV
public/                    interfaz web (JSME)
```
