# Contexto Experto del Repositorio: Tesis Hologenómica del Clado *Datura* (Dual RNA-seq + Metabolómica TA)

Este repositorio implementa y documenta el trabajo de tesis:
> **"Desarrollo de un flujo computacional reproducible para caracterizar la expresión de genes precursores y biosintéticos de la ruta de alcaloides tropánicos en plantas y sus microorganismos asociados en el clado Datura."**

---

## 1. Base de Conocimiento LLMWiki (`llmwiki_tesis_datura/`)

Antes de responder consultas metodológicas, escribir scripts de R/Bash/Python o redactar secciones del manuscrito, consulta siempre el grafo modular **LLMWiki** de 8 nodos ubicado en `llmwiki_tesis_datura/`:

- [`00_INDEX.md`](./llmwiki_tesis_datura/00_INDEX.md): Grafo maestro y arquitectura operativa de 4 módulos.
- [`01_MARCO_Y_OBJETIVOS.md`](./llmwiki_tesis_datura/01_MARCO_Y_OBJETIVOS.md): Planteamiento del problema (interacción holobionte, homología vs. expresión activa), objetivos (OE1–OE3) y bifurcación analítica del OE2 (Rama A dirigida por lecturas vs. Rama B *de novo* por contigs).
- [`02_MUESTRAS_Y_WETLAB.md`](./llmwiki_tesis_datura/02_MUESTRAS_Y_WETLAB.md): Inventario de 12 muestras, clasificación altitudinal (**Subclado A Alta montaña**: `B10`, `B11`, `B15`, `P1`; **Subclado B Tierras bajas**: `B9`, `B12`, `B14`, `B21`; **Clado Datura**: `B17`, `B27`; **Outgroups**: `B19` *Atropa belladonna*, `B20` *Ipomoea purpurea*), panel de genomas de referencia y protocolo Trizol/Illumina NovaSeq 6000 S4.
- [`03_MODULO1_SEPARACION_Y_ENSAMBLAJE.md`](./llmwiki_tesis_datura/03_MODULO1_SEPARACION_Y_ENSAMBLAJE.md): Preprocesamiento (`FastQC`, `MultiQC`, `Trimmomatic`, `RNA-QC-Chain`), filtro de Bloom (`FastqPuri v1.0.1`), co-ensamblaje paralelo en `rnaSPAdes` ($K=49, 73$; $123.099$ transcritos base) y validación estructural (`CD-HIT-EST`, `TransRate`, `Bowtie2`, `SAMtools`).
- [`04_MODULO2_PREDICCION_ORF_Y_ANOTACION.md`](./llmwiki_tesis_datura/04_MODULO2_PREDICCION_ORF_Y_ANOTACION.md): Diferenciación biológica en predicción de ORFs (`TransDecoder` para planta eucariota monocistrónica vs. `Prodigal` para bacteria procariota policistrónica/Shine-Dalgarno), base de referencia de 35 genes nominales de Alcaloides Tropánicos (TA) y anotación con `DIAMOND v2.1.16`.
- [`05_MODULO3_CUANTIFICACION_Y_DESEQ2.md`](./llmwiki_tesis_datura/05_MODULO3_CUANTIFICACION_Y_DESEQ2.md): Cuantificación en `Salmon v1.10.3`, corrección por longitud efectiva en `tximport` ($\log(E(K_{ij})) = \log(q_{ij}) + \log(L_{ij})$), embudo de filtrado en 5 pasos ($123.099 \rightarrow 54.363$ transcritos tras VST y filtro MAD del 75% superior) y expresión diferencial con `DESeq2`.
- [`06_MODULO4_METABOLOMICA_Y_COEXPRESION.md`](./llmwiki_tesis_datura/06_MODULO4_METABOLOMICA_Y_COEXPRESION.md): Curaduría metabolómica en 4 etapas ($82 \rightarrow 79 \rightarrow 73 \rightarrow 33 \rightarrow 15\text{–}25$ moléculas), comparativa de pipelines (Pablo vs. Ronald vs. Tesis actual), parametrización de `WGCNA`, correlación cruzada de *Module Eigengenes* ($\text{Planta} \leftrightarrow \text{Endófito}$) y métodos alternativos para *"Small $n$, Large $p$"* (`sPLS`/`DIABLO`, `MOFA+ v2`, Louvain, Mantel, CIA, `rCCA`, `ARACNE`).
- [`07_AUDITORIA_Y_REFERENCIAS.md`](./llmwiki_tesis_datura/07_AUDITORIA_Y_REFERENCIAS.md): Auditoría de errores críticos en matrices experimentales (`B27` vs. `B20`, 55 columnas vacías, $n=10$), inconsistencias del manuscrito `NewMethods.docx` y control bibliográfico.

---

## 2. Reglas Inquebrantables de Dominio y Bioinformática

1. **Validación de Cruce Muestral (`B27` vs. `B20`)**:
   - Nunca cruces la matriz de expresión génica con la matriz metabolómica sin verificar antes los identificadores de fila (`Sample ID`).
   - La matriz de expresión contiene **`B27`** (*Datura wrightii*), mientras que la tabla de fenotipos inicial contiene **`B20`** (*Ipomoea purpurea*). Exige siempre `stopifnot(identical(rownames(datExpr), rownames(datTraits)))` en scripts de R.
2. **Separación Estricta Eucariota vs. Procariota en Predicción de Proteínas**:
   - Usa exclusivamente **`TransDecoder`** (`LongOrfs -m 100` + `Predict`) para el transcriptoma vegetal (eucariota, monocistrónico).
   - Usa exclusivamente **`Prodigal`** para el transcriptoma bacteriano/endófito (procariota, operones policistrónicos, sitios *Shine-Dalgarno* y codones de inicio alternativos).
3. **Respeto al Embudo Estadístico Exacto ($135.008 \rightarrow 123.099 \rightarrow 54.363 \rightarrow 36.647$, `Normalizacion3.R`)**:
   - El ensamblaje inicial de `rnaSPAdes` produce **135.008 contigs crudos**, que tras `CD-HIT-EST` se reducen en un $8,8\%$ a **123.099 transcritos** ($98,3\%$ de alineamiento global).
   - En `Normalizacion3.R` (`Tabla_Embudo_Filtrado_Plantas_Corregido.csv`), el embudo opera en 7 pasos exactos:
     - Paso 0 (`tximport` post-`CD-HIT-EST`): **123.099** transcritos.
     - Paso 1 (`rowSums > 0`) y Paso 2 (`var > 0`): **122.690** retenidos ($409$ eliminados).
     - Paso 3 (`rowSums >= 25`): **54.435** retenidos ($68.255$ eliminados).
     - Paso 4 (`max <= 50.000`): **54.409** retenidos ($26$ eliminados).
     - Paso 5 (`max <= 30.000`) y Paso 6 (`!duplicated` + `vst(dds_clean, blind = TRUE)`): **54.363** retenidos ($46$ eliminados; este es el conteo exacto **pre-MAD**).
     - Paso 7 (Filtro MAD del $75\%$ superior, `mads > quantile(mads, 0.25)`): **36.647 transcritos finales** ($17.716$ eliminados) exportados en `matriz_vst_filtrada_para_WGCNA.csv`.
4. **Curaduría Metabolómica y Prevención de Corrupción Decimal Previa a Correlación**:
   - **¡Nunca uses `fenotipos_TA_limpios_para_WGCNA.csv`!** Ese archivo sufrió eliminación de puntos decimales (`.X`/`.XX`) en $339$ de las $790$ celdas respecto a `Metaboloma_WGNA_TESIS.xlsx`, multiplicando valores por $10\times$ o $100\times$ y distorsionando $61$ de las $79$ columnas (invirtiendo incluso el signo en `X14T`, `X64T`, `X61T`, `X15T`).
   - Utiliza siempre **`fenotipos_TA_corregidos_para_WGCNA.csv`** ($10 \times 79$) o **`fenotipos_TA_homologados_n8_SinOutgroups.csv`** ($8 \times 79$, muestras `B9, B10, B11, B12, B14, B15, B17, B21`), que contienen los **76 alcaloides tropánicos (`X6T..X81T`) + 3 precursores (`X82P..X84P`)** con sus decimales recuperados desde el XML de `Metaboloma_WGNA_TESIS.xlsx`.
   - Considera los clústeres de fragmentos en fuente / isómeros ($r > 0.85$, $79$ pares confirmados como `X7T`/`X18T`/`X62T`/`X10T` y `X15T`/`X67T`/`X71T`/`X63T`) y las matrices consolidadas de **33 moléculas** (o 15–25 post-MAD/OPLS-DA).
5. **Mitigación de *"Small $n$, Large $p$"* ($n=10$ / $n=9$ sin `B19` vs. $p=36.647$ / $54.363$)**:
   - En `WGCNA`, mantén siempre `networkType = "signed"`, `TOMType = "signed"`, `cor <- WGCNA::cor`, `deepSplit = 2.5`, `reassignThreshold = 0` y `pamRespectsDendro = FALSE`.
   - Excluye el *outgroup* `B19` (*Atropa belladonna*) en `WGCNA` ($n=9$: `B9`, `B10`, `B11`, `B12`, `B14`, `B15`, `B17`, `B21`, `B27`) con `power = 12`, configuración validada empíricamente sobre los **36.647 transcritos** de `matriz_vst_filtrada_para_WGCNA.csv` que produce **22 módulos funcionales** y solo **8,10% de genes en `grey`** (`MEpurple` específico de Alta Montaña, `MEblue`/`MEmagenta` de Tierras Bajas, y `MEbrown`/`MEred`/`MEcyan`/`MEtan` del clado *Datura*), eliminando los $614$ falsos positivos inducidos por `B19`.
   - Complementa siempre los hallazgos con intersección de `DESeq2` y métodos multivariados regularizados (`sPLS`/`DIABLO`, `MOFA+ v2`, Louvain en `igraph`, Test de Mantel o Co-inercia).
6. **Corrección de Inconsistencias Redaccionales**:
   - Al redactar o corregir el manuscrito (`NewMethods.docx`), usa siempre **`BLASTp`/`PSI-BLAST`** (nunca `BLASTn`) al describir la búsqueda con matriz `BLOSUM62` y `word_size = 3`, reporta los **35 identificadores génicos nominales** de la base TA, unifica los *k*-meros $K=49, 73$ en `rnaSPAdes` y explicita la exclusión de los *outgroups* (`B19`/`B20`) al pasar de 12 a 10 librerías.
