---
id_nodo: 05_MODULO3_CUANTIFICACION_Y_DESEQ2
tipo: modulo_operativo
modulo: 3
nodos_relacionados:
  - "[[00_INDEX]]"
  - "[[02_MUESTRAS_Y_WETLAB]]"
  - "[[03_MODULO1_SEPARACION_Y_ENSAMBLAJE]]"
  - "[[06_MODULO4_METABOLOMICA_Y_COEXPRESION]]"
---

# `05_MODULO3_CUANTIFICACION_Y_DESEQ2` — Salmon, tximport, Embudo VST y Expresión Diferencial

- **ID_Nodo**: `05_MODULO3_CUANTIFICACION_Y_DESEQ2`
- **Nodos_Relacionados**: [[00_INDEX]], [[02_MUESTRAS_Y_WETLAB]], [[03_MODULO1_SEPARACION_Y_ENSAMBLAJE]], [[06_MODULO4_METABOLOMICA_Y_COEXPRESION]]

---

## 1. Cuantificación en Salmon y Corrección Matemática en tximport

- **Cuantificación (`Salmon v1.10.3` en contenedor Singularity / Slurm)**: Indexación del ensamblaje agrupado (`salmon index` sobre `transcripts_clustered.fasta`, `123.099` targets válidos) seguida de cuantificación de lecturas pareadas (`salmon quant -l A --validateMappings -p 64`, formato de librería detectado `IU`), generando archivos `quant.sf` por muestra.
- **Métricas Empíricas de Mapeo en `Salmon v1.10.3` (`plantas/plantas/*_quant/`, `Salmon_Resumen_Mapeo_Plantas.csv`)**:

| Muestra (`Sample ID`) | Especie / Clado | Targets en Índice | Lecturas Procesadas | Lecturas Mapeadas | Tasa de Mapeo (`%`) | Formato Librería |
| :---: | :--- | :---: | :---: | :---: | :---: | :---: |
| **`B9`** | *Brugmansia insignis* (Tierras bajas) | $123.099$ | $14.630.017$ | $14.196.968$ | **$97,04\%$** | `IU` |
| **`B10`** | *Brugmansia sanguinea* (Alta montaña) | $123.099$ | $5.140.389$ | $5.022.373$ | **$97,70\%$** | `IU` |
| **`B11`** | *B. sanguinea x vulcanicola* (Alta montaña) | $123.099$ | $4.073.187$ | $3.976.613$ | **$97,63\%$** | `IU` |
| **`B12`** | *Brugmansia suaveolens* (Tierras bajas) | $123.099$ | $22.008.241$ | $21.428.362$ | **$97,37\%$** | `IU` |
| **`B14`** | *Brugmansia versicolor* (Tierras bajas) | $123.099$ | $15.531.420$ | $15.152.257$ | **$97,56\%$** | `IU` |
| **`B15`** | *Brugmansia vulcanicola* (Alta montaña) | $123.099$ | $6.886.496$ | $6.719.218$ | **$97,57\%$** | `IU` |
| **`B17`** | *Datura stramonium* (Clado *Datura*) | $123.099$ | $1.286.033$ | $1.246.537$ | **$96,93\%$** | `IU` |
| **`B19`** | *Atropa belladonna* (*Outgroup*) | $123.099$ | $127.327$ | $122.351$ | **$96,09\%$** | `IU` |
| **`B21`** | *Brugmansia aurea* (Tierras bajas) | $123.099$ | $9.347.912$ | $9.092.888$ | **$97,27\%$** | `IU` |
| **`B27`** | *Datura wrightii* (Clado *Datura*) | $123.099$ | $1.514.510$ | $1.478.894$ | **$97,65\%$** | `IU` |

- **Corrección de sesgo por longitud de transcrito (`tximport` en R)**: Los conteos observados dependen de la profundidad de secuenciación y del tamaño del gen; `tximport` aplica un ajuste (*offset*) basado en la longitud efectiva ($L_{ij}$) calculada por `Salmon` para que el tamaño del transcrito no distorsione el análisis estadístico:

$$\log(E(K_{ij})) = \log(q_{ij}) + \log(L_{ij})$$

Donde $K_{ij}$ representa los conteos del gen $i$ en la muestra $j$, $q_{ij}$ es la abundancia proporcional real y $L_{ij}$ es la longitud efectiva del transcrito.

---

## 2. Embudo Exacto de Filtrado y Normalización de Transcritos (`Normalizacion3.R`)

El script oficial `Normalizacion3.R` ejecuta el embudo sobre las 10 librerías vegetales (`Tabla_Embudo_Filtrado_Plantas_Corregido.csv`), resolviendo con exactitud matemática la relación entre los **135.008 contigs crudos**, los **123.099 transcritos post-clustering**, los **54.363 transcritos pre-MAD** y los **36.647 transcritos finales post-MAD**:

| Paso | Operación / Función en R (`Normalizacion3.R`) | Criterio de Filtrado | Propósito Estadístico y Biológico | Transcritos Retenidos | Transcritos Eliminados |
| :---: | :--- | :--- | :--- | :---: | :---: |
| **0** | Importación Base (`tximport` post-`CD-HIT-EST`) | Reducción desde $135.008$ contigs crudos ($-8,8\%$). | Integrar matriz de conteos corrigiendo sesgo por longitud efectiva $L_{ij}$. | **123.099** | $0$ |
| **1** | Eliminación de suma cero | `raw_counts[rowSums(raw_counts) > 0, ]` | Descartar transcritos sin expresión detectada en las 10 librerías. | **122.690** | $409$ |
| **2** | Eliminación de varianza cero | `apply(raw_counts, 1, var) > 0` | Descartar transcritos completamente constantes. | **122.690** | $0$ |
| **3** | Expresión biológica mínima | `raw_counts[rowSums(raw_counts) >= 25, ]` | Eliminar ruido técnico y transcritos marginales con suma $< 25$ lecturas. | **54.435** | $68.255$ |
| **4** | Filtro de extremos nivel 1 | `apply(raw_counts, 1, max) <= 50000` | Remover artefactos de hiper-expresión ribosomal/cloroplástica ($>50.000$). | **54.409** | $26$ |
| **5** | Filtro de extremos nivel 2 | `apply(raw_counts, 1, max) <= 30000` | Limpieza estricta de valores extremos ($>30.000$). | **54.363** | $46$ |
| **6** | Remoción de duplicados y VST | `!duplicated(rownames)` + `vst(dds_clean, blind = TRUE)` | Comprimir la sobredispersión de la Binomial Negativa para lograr homocedasticidad. | **54.363** | $0$ |
| **7** | Filtro de Varianza (MAD) | `vst_matrix[mads > quantile(mads, 0.25), ]` | Retener los transcritos sobre el cuantil $25\%$ de Desviación Absoluta de la Mediana. | **36.647** (Final) | $17.716$ |

### Matrices Generadas tras el Embudo (`Normalizacion3.R`)
- **`Tabla_Embudo_Filtrado_Plantas_Corregido.csv`**: Registro exacto de retención y pérdida de genes por etapa ($123.099 \rightarrow 54.363 \rightarrow 36.647$).
- **`matriz_vst_filtrada_genes_x_muestras.csv`** / **`Matriz_VST_Visualizacion_Limpia.xlsx`**: $36.647\text{ Transcritos en filas} \times 10\text{ Muestras en columnas}$.
- **`matriz_vst_filtrada_para_WGCNA.csv`**: Matriz transpuesta ($10\text{ Muestras en filas} \times 36.647\text{ Transcritos en columnas}$).

---

## 3. Motor de Expresión Diferencial (`DESeq2`) y Resultados de Agrupamiento

### Secuencia Algorítmica en `DESeq2`
1. Importación con `DESeqDataSetFromTximport`.
2. Modelado mediante **Distribución Binomial Negativa** para manejar la sobredispersión de conteos de RNA-Seq.
3. **Contracción de la dispersión (*Shrinkage*)** hacia la curva media de expresión para evitar falsos positivos por ruido en genes de conteo bajo.
4. Prueba de razón de verosimilitud (LRT) / contrastes de $\log_2(\text{Fold Change})$ vs. $-\log_{10}(p\text{-valor})$ visualizados mediante *Volcano Plots* entre **Subclado A (Alta montaña)** y **Subclado B (Tierras bajas)**.

### Patrones Biológicos Validados en Control de Calidad (PCA, Dendrograma Euclidiano y Heatmap de Pearson)
- **Set Vegetal**: Separación neta en 3 clústeres:
  1. *Brugmansia* de Tierras bajas (`B9`, `B12`, `B14`, `B21`): Alta cohesión interna en correlación de Pearson sobre VST ($r = 0.7146\text{–}0.8824$; máximo entre `B9` y `B12` con $r = 0.8824$).
  2. *Brugmansia* de Tierras altas (`B10`, `B11`, `B15`, `P1`): Alta cohesión interna ($r = 0.7761\text{–}0.8580$; máximo entre `B10` y `B11` con $r = 0.8580$).
  3. Clado *Datura* (`B17`, `B27`, correlación mutua $r = 0.7608$) junto con el grupo externo *Atropa belladonna* (`B19`) situado por fuera (correlación baja frente a todas las muestras del clado: $r = 0.1929\text{–}0.3311$).
- **Métricas Empíricas de `matriz_vst_filtrada_para_WGCNA.csv` ($10 \times 36.647$)**:
  - **Integridad (`goodSamplesGenes`)**: $0$ valores `NA`/`Inf`, $0$ transcritos con varianza nula, rango VST $[0.8729, 16.7079]$, mediana de $\text{SD} = 1.8845$, mediana de $\text{MAD} = 1.3978$.
  - **Varianza explicada en PCA**:
    - Con `B19` ($n=10$): $\text{PC1} = 54.99\%$, $\text{PC2} = 13.18\%$, $\text{PC3} = 12.27\%$.
    - Sin `B19` ($n=9$, solo *Datura*/*Brugmansia*): $\text{PC1} = 54.00\%$, $\text{PC2} = 18.39\%$ (incrementa $+5.21\%$ el poder discriminante entre Alta Montaña y Tierras Bajas), $\text{PC3} = 7.08\%$.
- **Set Bacteriano**: Replica exactamente la misma arquitectura filogenética y ecológica del hospedero, dividiendo en el PCA la microbiota asociada en nichos de tierras bajas (`B9`, `B12`, `B14`, `B21`), tierras altas (`B10`, `B11`, `B15`, `P1`) y grupo externo.
