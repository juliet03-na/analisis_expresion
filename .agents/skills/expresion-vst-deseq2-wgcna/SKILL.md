---
name: expresion-vst-deseq2-wgcna
description: >-
  Especialista estadístico y bioinformático para los Módulos 3 y 4 de la tesis
  hologenómica en el clado Datura: cuantificación en Salmon, corrección por longitud
  efectiva en tximport, embudo de filtrado estadístico en 5 pasos (de 123.099 a 54.363
  transcritos vía VST y MAD), expresión diferencial en DESeq2 (Subclado A Alta montaña
  vs. Subclado B Tierras bajas), curaduría metabolómica en 4 etapas (de 82 a 33 y 15-25
  moléculas), redes de co-expresión WGCNA, correlación cruzada de Module Eigengenes
  planta-endófito y modelos multivariados para "Small n, Large p" (sPLS/DIABLO, MOFA+ v2,
  Louvain, Mantel, CIA, rCCA, ARACNE). Úsese para escribir o depurar código en R de
  expresión diferencial, normalización, metabolómica o co-expresión.
---

# Cuantificación, Embudo VST, DESeq2, Curaduría Metabolómica y Redes WGCNA / Multi-Ómicas (Módulos 3 y 4)

Esta *skill* implementa los flujos estadísticos de [[05_MODULO3_CUANTIFICACION_Y_DESEQ2]] y [[06_MODULO4_METABOLOMICA_Y_COEXPRESION]].

## 1. Módulo 3: Cuantificación (`Salmon`) y Corrección (`tximport`)

1. **Cuantificación en `Salmon v1.10.3`**:
   - Indexación de `transcripts_clustered.fasta` (`salmon index`).
   - Cuantificación pareada: `salmon quant -i index -l A -1 R1.fq.gz -2 R2.fq.gz --validateMappings -o sample_quant`.
2. **Importación y Offset en R (`tximport`)**:
   - Corrige el sesgo de longitud efectiva $L_{ij}$ según:
     $$\log(E(K_{ij})) = \log(q_{ij}) + \log(L_{ij})$$

## 2. Módulo 3: Embudo Estadístico Exacto (`Normalizacion3.R`: $135.008 \rightarrow 123.099 \rightarrow 54.363 \rightarrow 36.647$)

Al programar o auditar el preprocesamiento en R (`Normalizacion3.R`, resultados en `Tabla_Embudo_Filtrado_Plantas_Corregido.csv`), aplica estrictamente este orden:
1. **Paso 0 — Importación Base (`tximport`)**: **123.099** transcritos del ensamblaje agrupado por `CD-HIT-EST` (provenientes de **135.008** contigs crudos de `rnaSPAdes`, reducción del $8,8\%$ y tasa de alineamiento global del $98,3\%$).
2. **Paso 1 y 2 — Eliminación de ceros y varianza nula**: `rowSums(raw_counts) > 0` y `var > 0` $\rightarrow$ **122.690** retenidos ($409$ eliminados).
3. **Paso 3 — Expresión biológica mínima**: `rowSums(raw_counts) >= 25` $\rightarrow$ **54.435** retenidos ($68.255$ eliminados).
4. **Paso 4 — Limpieza de extremos nivel 1**: `max(raw_counts) <= 50000` $\rightarrow$ **54.409** retenidos ($26$ eliminados).
5. **Paso 5 y 6 — Limpieza de extremos nivel 2, duplicados y VST**: `max(raw_counts) <= 30000`, `!duplicated(rownames)` y `vst(dds_clean, blind = TRUE)` $\rightarrow$ **54.363** transcritos retenidos ($46$ eliminados; conteo exacto **pre-MAD**).
6. **Paso 7 — Filtro de Varianza por MAD (75% superior)**:
   - Calcular `mads <- apply(vst_matrix, 1, mad)` y conservar `vst_matrix[mads > quantile(mads, 0.25), ]`, obteniendo exactamente **36.647 transcritos finales** ($17.716$ eliminados).
   - Exportar `Tabla_Embudo_Filtrado_Plantas_Corregido.csv`, `matriz_vst_filtrada_genes_x_muestras.csv` ($36.647 \times 10$) y `matriz_vst_filtrada_para_WGCNA.csv` ($10 \times 36.647$).

## 3. Módulo 4: Curaduría Metabolómica en 4 Etapas ($82 \rightarrow 79 \rightarrow 33 \rightarrow 15\text{–}25$)

Antes de correlacionar con expresión génica:
1. **Regla de Integridad Decimal (`Metaboloma_WGNA_TESIS.xlsx` vs. `fenotipos_TA_limpios_para_WGCNA.csv`)**:
   - **Prohibido usar `fenotipos_TA_limpios_para_WGCNA.csv`**: se borraron los puntos decimales de $339$ celdas (`.X`/`.XX`), multiplicándolas por $10\times$ o $100\times$ y distorsionando $61/79$ columnas.
   - Usar siempre **`fenotipos_TA_corregidos_para_WGCNA.csv`** ($10 \times 79$) o **`fenotipos_TA_homologados_n8_SinOutgroups.csv`** ($8 \times 79$: `B9, B10, B11, B12, B14, B15, B17, B21`), que contienen **76 alcaloides tropánicos (`X6T..X81T`) + 3 precursores (`X82P..X84P`)** con sus decimales recuperados desde el XML de `Metaboloma_WGNA_TESIS.xlsx`.
2. **Etapa 1 ($79 \rightarrow 76\text{ TA} + 3\text{ P}$)**: Separar los 3 precursores (`X82P`, `X83P`, `X84P`) de los 76 alcaloides tropánicos (`X6T`–`X81T`).
3. **Etapa 2 ($76 \rightarrow 73$)**: Remover fragmentos en fuente (*in-source fragments*) evaluando RT, MS2 redundantes y correlación de Pearson $r > 0.85$ entre picos de igual masa ($79$ pares confirmados en los datos reales, como `X7T`/`X18T`/`X62T`/`X10T` con $r=0.98\text{–}0.9998$ y `X15T`/`X67T`/`X71T`/`X63T` con $r=0.986\text{–}0.9936$).
4. **Etapa 3 ($73 \rightarrow 33$ moléculas)**: Colapsar regioisómeros y estereoisómeros (*Feature clustering* por fórmula exacta: 15 de pico único + 18 grupos de isómeros) en dos matrices evaluables:
   - **Matriz de abundancia acumulada** (suma de áreas de isómeros).
   - **Matriz de isómero principal** (pico representativo de mayor abundancia).
5. **Etapa 4 ($33 \rightarrow 15\text{–}25$ alcaloides clave)**: Filtro MAD (elimina el 40% de baja varianza) + OPLS-DA conservando metabolitos con $\text{VIP} > 1.2 - 1.5$.
6. **Asociaciones de Referencia Validadas en Fase 2 ($n=8$, $101$ pares con $\|r\| \ge 0.70, p \le 0.05$)**:
   - **Eje `X70T`**: Positivo en Tierras Bajas (`MEblue` $r = +0.9213$, `MEmagenta` $r = +0.9146$, `MEturquoise` $r = +0.8893$) vs. negativo en *Datura* (`MEred` $r = -0.9732, p = 0.00005$).
   - **Eje `X52T`**: Negativo en Alta Montaña (`MEpurple` $r = -0.9711, p = 0.00006$) vs. positivo en *Datura* (`MEbrown` $r = +0.8605, p = 0.00609$).
   - **Precursor `X82P`**: Asociado positivamente con `MEmagenta` ($r = +0.9546$ lineal, $r = +0.8758$ en $\log_{10}$, $p = 0.00435$).

## 4. Módulo 4: Parametrización de `WGCNA` y Diagnóstico (`blockwiseModules`)

- **Validación obligatoria antes de correr `WGCNA`**:
  ```r
  # Homologar muestras y evitar el error B27 (expresión) vs B20 (metabolómica)
  common_samples <- intersect(rownames(datExpr), rownames(datTraits))
  datExpr <- datExpr[common_samples, , drop = FALSE]
  datTraits <- datTraits[common_samples, , drop = FALSE]
  stopifnot(identical(rownames(datExpr), rownames(datTraits)))
  ```
- **Parámetros Oficiales Validados Empíricamente (`matriz_vst_filtrada_para_WGCNA.csv`, $36.647$ transcritos)**:
  - **Exclusión de `B19` (*Atropa belladonna*)**: Con `B19` ($n=10$), la divergencia taxonómica ($r = 0.19\text{–}0.33$) secuestra $3.222$ transcritos en un módulo artificial. Al excluir `B19` ($n=9$: `B9`, `B10`, `B11`, `B12`, `B14`, `B15`, `B17`, `B21`, `B27`) con `power = 12` y etapa PAM (`pamRespectsDendro = FALSE`), se obtienen **22 módulos funcionales** y **`grey = 8.10%`** ($2.970$ transcritos), descartando $614$ falsos positivos inducidos por `B19`.
  - **Módulos de Referencia (`Sin B19`, $n=9$)**:
    - `MEpurple` ($736$ genes) y `MEgreenyellow` ($704$ genes): Específicos de **Subclado A (Alta Montaña: `B10`, `B11`, `B15`)**.
    - `MEblue` ($2.345$ genes) y `MEmagenta` ($806$ genes): Específicos de **Subclado B (Tierras Bajas: `B9`, `B12`, `B14`, `B21`)**.
    - `MEturquoise` ($17.694$ genes): Divergencia *Brugmansia* (positivo) vs. *Datura* (negativo).
    - `MEbrown` ($1.736$ genes) y `MEred` ($1.427$ genes): Activación en *Datura* (`B17`, `B27`) en oposición a Alta Montaña (`MEbrown`) o Tierras Bajas (`MEred`).
    - `MEcyan` ($604$ genes, `B17` *D. stramonium*) y `MEtan` ($685$ genes, `B27` *D. wrightii*): Especie-específicos en *Datura*.
  - `cor <- WGCNA::cor` antes de `blockwiseModules()` para prevenir el error `unused argument (weights.x = NULL)`.
  - `power = 12` (o `6` / `16` para comparativa), `networkType = "signed"`, `TOMType = "signed"`, `corType = "pearson"`.
  - `minModuleSize = 30`, `deepSplit = 2.5`, `reassignThreshold = 0`, `mergeCutHeight = 0.25`, `maxBlockSize = 10000` (4 bloques para $36.647$ transcritos).
  - **Criterio módulo `grey`**: $< 20\% - 30\%$ (máximo $40\%$) del total de transcritos (logrado: `8.10%`).
  - **Selección de candidatos**: Correlación Módulo-Alcaloide $\|r\| \ge 0.65 - 0.70$ ($p < 0.05$) + alto MM/GS + cruce con `DESeq2` y anotación `DIAMOND`.

## 5. Integración Hologenómica y Alternativas para *"Small $n$, Large $p$"* ($n=10$)

1. **Correlación Cruzada Planta $\leftrightarrow$ Endófito**: Correlacionar los *Module Eigengenes* (`MEs`) de la Red Vegetal contra los `MEs` de la Red Microbiana.
2. **Métodos Alternativos para $n=10$ vs. $p=54.363$**:
   - `mixOmics` (`sPLS` / `DIABLO`, `rCCA` con transformación CLR) y *Elastic Net*.
   - `MOFA2` (`MOFA+ v2`) para descomposición en factores latentes multi-ómicos.
   - Redes por umbrales duros + detección de comunidades con **Algoritmo de Louvain** en `igraph`.
   - **Test de Mantel** (`vegan::mantel`) y **Análisis de Co-inercia (CIA)** (`ade4::coinertia`).
   - `ARACNE` (Información Mutua) y *Targeted Correlation* (Pearson + FDR sobre los 35 genes candidatos TA).
