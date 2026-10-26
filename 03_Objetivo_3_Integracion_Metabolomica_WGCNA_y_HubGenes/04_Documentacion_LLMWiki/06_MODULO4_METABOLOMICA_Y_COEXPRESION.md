---
id_nodo: 06_MODULO4_METABOLOMICA_Y_COEXPRESION
tipo: modulo_operativo
modulo: 4
nodos_relacionados:
  - "[[00_INDEX]]"
  - "[[04_MODULO2_PREDICCION_ORF_Y_ANOTACION]]"
  - "[[05_MODULO3_CUANTIFICACION_Y_DESEQ2]]"
  - "[[07_AUDITORIA_Y_REFERENCIAS]]"
---

# `06_MODULO4_METABOLOMICA_Y_COEXPRESION` — Curaduría MS, WGCNA, Redes Duales y sPLS/MOFA+

- **ID_Nodo**: `06_MODULO4_METABOLOMICA_Y_COEXPRESION`
- **Nodos_Relacionados**: [[00_INDEX]], [[04_MODULO2_PREDICCION_ORF_Y_ANOTACION]], [[05_MODULO3_CUANTIFICACION_Y_DESEQ2]], [[07_AUDITORIA_Y_REFERENCIAS]]

---

## 1. Protocolo de Curaduría y Reducción de Datos Metabolómicos

- **Problema inicial de la matriz cruda y estado de `Metaboloma_WGNA_TESIS.xlsx`**:
  - El archivo crudo original constaba de **82 compuestos** y presentaba **55 columnas vacías/con ceros**.
  - En el archivo depurado **`Metaboloma_WGNA_TESIS.xlsx`** (hoja `MetabolomaTA`, $10\text{ muestras} \times 79\text{ compuestos}$, sin ceros ni celdas vacías), las variables se encuentran codificadas como **76 alcaloides tropánicos (`6T` a `81T`)** y **3 precursores (`82P`, `83P`, `84P`)**.
  - **Corrección obligatoria del punto decimal (`fenotipos_TA_corregidos_para_WGCNA.csv`)**: En `Metaboloma_WGNA_TESIS.xlsx`, $339$ de las $790$ celdas ($42,9\%$) estaban almacenadas como texto con 1 o 2 decimales (`.X` y `.XX`, ej. `33728376.94`). En `fenotipos_TA_limpios_para_WGCNA.csv` se habían eliminado los puntos `.`, multiplicando esas celdas por $10\times$ o $100\times$ y distorsionando $61$ de las $79$ columnas (ver [[07_AUDITORIA_Y_REFERENCIAS]]). Toda correlación debe ejecutarse sobre **`fenotipos_TA_corregidos_para_WGCNA.csv`** o **`fenotipos_TA_homologados_n8_SinOutgroups.csv`**.
- **Flujo de curaduría en 4 etapas** (adaptado del modelo de Ronald a los datos del proyecto):
  1. **Separación de precursores no-alcaloides (de 79 a 76 TA + 3 P)**: Distinción explícita entre los $76$ alcaloides tropánicos (`X6T`–`X81T`) y los $3$ precursores (`X82P`, `X83P`, `X84P`).
  2. **Remoción de fragmentos en fuente (*in-source fragments*) (reduce a 73)**: Evaluación de tiempos de retención (RT), espectros MS2 redundantes y matriz de correlación de Pearson entre picos de igual masa ($r > 0.85$).
     - *Validación empírica en `Metaboloma_WGNA_TESIS.xlsx`*: Se confirmaron **79 pares de picos con $r > 0.85$**, destacando los clústeres de fragmentación/isomería:
       - **Clúster `X18T`/`X7T`/`X62T`/`X10T`**: `X7T` $\leftrightarrow$ `X18T` ($r = 0.9998$), `X18T` $\leftrightarrow$ `X62T` ($r = 0.9996$), `X7T` $\leftrightarrow$ `X62T` ($r = 0.9989$), `X7T` $\leftrightarrow$ `X10T` ($r = 0.9820$).
       - **Clúster `X71T`/`X15T`/`X67T`/`X63T` (Escopolamina)**: `X15T` $\leftrightarrow$ `X67T` ($r = 0.9936$), `X15T` $\leftrightarrow$ `X71T` ($r = 0.9913$), `X63T` $\leftrightarrow$ `X71T` ($r = 0.9890$), `X67T` $\leftrightarrow$ `X71T` ($r = 0.9866$).
       - **Clúster `X19T`/`X46T`/`X58T`**: `X19T` $\leftrightarrow$ `X46T` ($r = 0.9869$), `X46T` $\leftrightarrow$ `X58T` ($r = 0.9814$).
       - **Clúster `X59T`/`X60T`**: `X59T` $\leftrightarrow$ `X60T` ($r = 0.9782$).
  3. **Colapso de regioisómeros y estereoisómeros (*Feature clustering*) (consolida en 33 moléculas)**: Agrupamiento jerárquico por fórmula molecular exacta e identidad química, consolidando los $73$ picos en **33 moléculas** ($15$ de pico único y $18$ grupos de isómeros). Se generan dos matrices evaluables:
     - **Matriz de abundancia acumulada**: Suma de las áreas de todos los isómeros de la molécula.
     - **Matriz de isómero principal**: Retiene únicamente el pico representativo de mayor abundancia descartando ruido menor.
  4. **Filtrado de varianza y discriminación biológica (retiene 15–25 alcaloides clave)**: Aplicación de filtro MAD (eliminando el $40\%$ de baja varianza en MetaboAnalyst) y selección supervisada mediante OPLS-DA conservando metabolitos con puntaje $\text{VIP} > 1.2 - 1.5$ capaces de separar grupos biológicos.

---

## 2. Comparativa de Referencia: Pipeline Pablo vs. Pipeline Ronald vs. Tesis Actual

| Característica | Pipeline Pablo (Referencia Expresión) | Pipeline Ronald (Referencia Metabolómica) | Configuración Tesis Actual |
| :--- | :--- | :--- | :--- |
| **Enfoque Principal** | Integración de `WGCNA` con expresión diferencial clásica (`DESeq2`/`edgeR`) y búsqueda de ortogrupos. | Análisis multi-ómico correlacionando `WGCNA` directamente con metaboloma curado en Solanaceae. | Hologenómica (Dual RNA-seq Planta + Endófitos) integrando `WGCNA`, `DESeq2` y curaduría de 33 moléculas TA. |
| **Datos Fenotípicos** | Valores morfológicos y resumen PCA de alcaloides. | 102 metabolitos en tesis; 47 en heatmap `WGCNA` (25–27 tropanos, 10 glicoalcaloides SGAs, 10 flavonoides). | Matriz curada de alcaloides tropánicos (33 moléculas consolidadas / 15–25 post-filtro MAD) y PCA de alcaloides. |
| **Parámetros WGCNA** | `power = 14`, `deepSplit = 2.0`, `mergeCutHeight = 0.15` (redes más grandes y fusionadas). | `power = 6`, `deepSplit = 2.5`, `mergeCutHeight = 0.25` (redes segmentadas y específicas). | `power = 6` (o `12–16` según validación), `TOMType = "signed"`, `corType = "pearson"`, `deepSplit = 2.5`, `minModuleSize = 20–30`, `mergeCutHeight = 0.25`, `reassignThreshold = 0`. |
| **Selección de Candidatos** | Intersección estricta: Módulo clave + Alto MM/GS + $\log_2\text{FC}$ significativo + FDR bajo. | Correlación directa de *Module Eigengenes* con la acumulación de alcaloides específicos. | Correlación Módulo-Alcaloide ($\|r\| \ge 0.65 - 0.70$, $p < 0.05$) + Alto MM/GS + cruce con expresión diferencial (`DESeq2`) y anotación `DIAMOND`. |

---

## 3. Parametrización, Scripts y Matriz de Validación de WGCNA

### Flujo de Scripts en R
- **`WGCNA_Phase1_Prueba1.R`**: Recibe `matriz2_vst_filtrada` (prueba reducida de $2.000$ a $5.000$ transcritos en $8\text{ GB}$ RAM) $\rightarrow$ evalúa topología libre de escala (`Topologia.pdf`) y construye árbol jerárquico (`Archivo_Rdata_WGCNA`).
- **`WGCNA_Phase1_Prueba2.R` y `WGCNA_4_Mejorar_Modulos.R`**: Limpia columnas vacías de fenotipos y cruza módulos contra alcaloides $\rightarrow$ genera `WGCNA4_correlacion_mejorado.pdf`.
- **`WGCNA_All_Data.R`**: Procesa los $54.363$ transcritos en bloques (`WGCNA.RData`).

### Matriz de Parametrización y Diagnóstico en `blockwiseModules()`

| Parámetro en `blockwiseModules()` | Configuración | Función Matemática y Biológica | Criterio de Éxito / Diagnóstico | Acción Correctiva si Falla |
| :--- | :--- | :--- | :--- | :--- |
| **`power` ($\beta$)** | `6` (evaluar `12–16` o `18` por $n<20$) | Penaliza correlaciones débiles para alcanzar topología libre de escala. | Curva `pickSoftThreshold()` alcanza $R^2 \ge 0.80 - 0.85$ sin caer conectividad a cero. | Si no alcanza $0.85$ por heterogeneidad muestral, usar potencia recomendada para redes con signo (`12–16`). |
| **`networkType` / `TOMType`** | `"signed"` | Agrupa genes solo si co-varían en la misma dirección (se encienden y apagan juntos). | Transcritos con correlación opuesta no se mezclan en el mismo módulo. | Mantener fijo en `"signed"`. |
| **`minModuleSize`** | `20` (pruebas) / `30–50` (54k genes) | Número mínimo de transcritos para constituir un módulo independiente. | Evitar micro-módulos sin sentido biológico. | Si resultan $>40$ módulos, elevar el umbral a `50` o `100`. |
| **`deepSplit`** | `2.5` (o `2`) | Sensibilidad de corte de ramas en el dendrograma (rango $0$ a $4$). | Bloques definidos y compactos en el dendrograma. | Si el árbol queda hiper-fragmentado, bajar a `2.0`. |
| **`mergeCutHeight`** | `0.25` | Umbral de disimilitud ($1 - r$) que fusiona módulos con *Eigengenes* correlacionados en $r \ge 0.75$. | Ausencia de módulos contiguos con patrones idénticos en el heatmap. | Si persisten módulos redundantes, incrementar el corte a `0.30`. |
| **`maxBlockSize`** | `30000` (HPC) / `10000` (16 GB RAM) | Límite de genes por bloque para el cálculo de la matriz TOM en memoria. | Ejecución sin error `cannot allocate vector` en R. | Ajustar según RAM disponible. |
| **Módulo `grey`** | Post-ejecución (`table(moduleColors)`) | Agrupa transcritos sin patrón de co-expresión coordinado. | Debe contener $< 20\% - 30\%$ (máximo $40\%$) de los $54.363$ transcritos. | Si domina el gráfico, reducir `power` o relajar `minModuleSize`. |
| **Heatmap Módulo-Rasgo** | Correlación Pearson $r$ y $p$-valor ($t$-Student) | Relaciona cada *Module Eigengene* con la concentración de alcaloides (Rojo $\approx +1$, Azul $\approx -1$, Blanco $\approx 0$). | Magnitud $\|r\| \ge 0.65 - 0.70$ y significancia estadística $p < 0.05$. | Aplicar curaduría de metabolitos (reducción a 33 / 15–25 moléculas) o métodos regularizados para $n=10$. |

---

## 4. Integración Hologenómica y Métodos Alternativos para *"Small $n$, Large $p$"*

- **Integración Cruzada de *Eigengenes* ($\text{Planta} \leftrightarrow \text{Endófito}$)**: Se construyen de forma independiente la Red Vegetal y la Red Microbiana, y posteriormente se correlacionan los *Module Eigengenes* de los módulos clave de la planta contra los *Module Eigengenes* de los endófitos para identificar activación o represión metabólica coordinada entre el hospedero y su microbiota.
- **Batería de Métodos Estadísticos Alternativos a WGCNA (Mitigación de $n=10$ muestras frente a $p=54.363$ transcritos)**:
  - **`sPLS` / `DIABLO` (*sparse Partial Least Squares*) y *Elastic Net***: Métodos multivariados supervisados con penalización que seleccionan el subconjunto exacto de genes que maximiza la covarianza entre transcriptoma, metaboloma y microbioma.
  - **`MOFA+ v2` (*Multi-Omics Factor Analysis*)**: Equivalente directo a los *Module Eigengenes* diseñado para $n$ pequeño; descompone simultáneamente las matrices de planta, microbioma y metaboloma en Factores latentes ("súper módulos").
  - **Redes por Umbrales Duros + Algoritmo de Louvain (en `igraph`)**: Construye grafos de correlación y detecta comunidades/módulos sin forzar artificialmente una topología libre de escala cuando `pickSoftThreshold` es inestable.
  - **Test de Mantel y Análisis de Co-inercia (CIA)**: Evalúan la correlación estadística y covarianza entre matrices de distancia completas, siendo altamente robustos con tamaños muestrales reducidos ($n=10$).
  - **`rCCA` (*Regularized Canonical Correlation Analysis*) y Ordenación Canónica (con transformación CLR)**: Maximizan la correlación cruzada entre subconjuntos de transcritos y perfiles cuantitativos de metabolitos.
  - **`ARACNE`**: Reconstruye redes de co-expresión regulatorias utilizando Información Mutua en lugar de correlación lineal de Pearson.
  - ***Targeted Correlation***: Evalúa exclusivamente la correlación de Pearson (con corrección FDR) sobre los genes candidatos ya conocidos de la ruta TA ([[04_MODULO2_PREDICCION_ORF_Y_ANOTACION]]).

---

## 5. Resultados Empíricos de Fase 1 (`WGCNA` sobre $36.647$ Transcritos: Con `B19` vs. Sin `B19`)

### 5.1 Diagnóstico de Topología Libre de Escala (`pickSoftThreshold` en Red `signed`)
- Debido al fuerte contraste filogenético entre *Brugmansia* (Tierras Bajas vs. Alta Montaña), *Datura* y *Atropa* con $n \le 10$, miles de transcritos co-varían en bloque por clado, generando una pendiente positiva en $\log_{10}(p(k))$ vs. $\log_{10}(k)$ y alta conectividad media en $\beta = 6$ ($k_{\text{medio}} = 8.374,08$ con `B19`; $6.727,46$ sin `B19`).
- Siguiendo el protocolo de contingencia para redes `signed` con $n < 15$, elevar la potencia a **`power = 12`** sin `B19` reduce la conectividad media a **`2.880,79`** ($1.900,95$ en $\beta = 16$), limpiando correlaciones espurias de fondo.

### 5.2 Impacto de Excluir el *Outgroup* `B19` (*Atropa belladonna*) y Activar Etapa PAM (`pamRespectsDendro = FALSE`)

| Indicador | Con `B19` ($n=10$, $\beta=6$, sin PAM) | **Sin `B19` ($n=9$, $\beta=12$, con PAM)** | Diagnóstico Metodológico |
| :--- | :---: | :---: | :--- |
| **Transcritos evaluados** | $36.647$ ($4$ bloques $\le 10.000$) | **$36.647$ ($4$ bloques $\le 10.000$)** | $0$ transcritos pierden varianza al excluir `B19`. |
| **Módulos funcionales (excl. `grey`)** | $4$ | **$22$** | Resolución óptima por nicho ecológico y especie. |
| **Porcentaje en módulo `grey`** | $67,39\%$ ($24.697$ transcritos) | **$8,10\%$ ($2.970$ transcritos)** | Cumple holgadamente el umbral $< 20\% - 30\%$. |
| **Sesgo por *outgroup* `B19`** | El módulo `brown` ($3.222$ genes) captura exclusivamente la divergencia de `B19` ($\text{ME} = -2.84$). | De esos $3.222$ genes, **$614$ ($19,06\%$)** caen a `grey` (falsos positivos inducidos por `B19`) y **$2.608$** se redistribuyen en módulos reales del clado *Datura*. | Excluir `B19` es obligatorio antes de correlacionar con alcaloides tropánicos. |

### 5.3 Catálogo de Módulos Clave y Significado Biológico (Red Depurada Sin `B19`, $n=9$)

| Módulo (`Sin B19`) | N° Transcritos (`%`) | Perfil del *Module Eigengene* (`ME`) por Muestra | Interpretación Biológica para Cruce con Alcaloides (Fase 2) |
| :--- | :---: | :--- | :--- |
| **`MEturquoise`** | $17.694$ ($48,28\%$) | Positivo en *Brugmansia* (`B9`: $+1.05$, `B12`: $+1.03$, `B14`: $+0.57$, `B21`: $+0.47$); muy negativo en *Datura* (`B17`: $-1.62$, `B27`: $-1.62$). | **Macro-módulo de divergencia morfológica/metabólica *Brugmansia* (arbórea) vs. *Datura* (herbácea)**. |
| **`MEblue`** | $2.345$ ($6,40\%$) | Alto en **Tierras Bajas** (`B9`: $+1.00$, `B12`: $+0.82$, `B14`: $+1.02$, `B21`: $+0.98$); negativo en **Alta Montaña** (`B10`: $-0.90$, `B11`: $-0.93$) y *Datura* ($-1.13$). | **Módulo maestro del Subclado B (Tierras Bajas)**. |
| **`MEbrown`** | $1.736$ ($4,74\%$) | Alto en *Datura* (`B17`: $+1.30$, `B27`: $+1.52$); muy negativo en **Alta Montaña** (`B10`: $-1.15$, `B11`: $-1.10$, `B15`: $-1.13$). | **Oposición transcripcional Clado *Datura* vs. Subclado A (Alta Montaña)**. |
| **`MEyellow`** | $1.573$ ($4,29\%$) | Alto en **Alta Montaña** (`B10`: $+1.06$, `B11`: $+1.07$, `B15`: $+0.75$) y `B9`/`B12`; negativo en `B14`/`B21` y *Datura* ($-1.45$ a $-1.47$). | **Co-expresión compartida en *Brugmansia* andinas vs. *Datura***. |
| **`MEred`** | $1.427$ ($3,89\%$) | Alto en *Datura* (`B17`: $+1.49$, `B27`: $+1.43$) y positivo en Alta Montaña; negativo en **Tierras Bajas** (`B9`: $-1.14$, `B12`: $-1.04$, `B14`: $-0.78$, `B21`: $-0.74$). | **Oposición transcripcional Clado *Datura* vs. Subclado B (Tierras Bajas)**. |
| **`MEmagenta`** | $806$ ($2,20\%$) | Muy alto en **Tierras Bajas** (`B14`: $+1.67$, `B9`: $+1.03$, `B12`: $+0.93$); negativo en Alta Montaña y *Datura*. | **Módulo específico de *Brugmansia* de Tierras Bajas (con pico en *B. versicolor* `B14`)**. |
| **`MEpurple`** | $736$ ($2,01\%$) | **Exclusivamente alto en Alta Montaña** (`B10`: $+1.38$, `B11`: $+1.41$, `B15`: $+1.11$); negativo en todas las de Tierras Bajas ($-0.40$ a $-0.57$) y *Datura* ($-0.90$ a $-1.03$). | **Módulo maestro de adaptación a Alta Montaña (Subclado A)**. |
| **`MEgreenyellow`** | $704$ ($1,92\%$) | Positivo en Alta Montaña (`B10`: $+0.89$, `B11`: $+0.94$, `B15`: $+0.83$) y `B14`/`B21`; negativo en `B9`/`B12` y *Datura*. | **Módulo secundario de tierras altas y transición**. |
| **`MEtan` vs. `MEcyan`** | $685$ y $604$ | `MEtan` es exclusivo de **`B27` (*D. wrightii*: $+2.33$)**; `MEcyan` es exclusivo de **`B17` (*D. stramonium*: $+2.47$)**. | **Módulos especie-específicos dentro del género *Datura***. |

---

## 6. Resultados Empíricos de Fase 2: Correlación Módulos (`WGCNA Sin B19`) vs. Alcaloides Tropánicos ($n=8$)

Al homologar las muestras comunes del clado *Datura*/*Brugmansia* entre expresión y metabolómica (**$n=8$: `B9`, `B10`, `B11`, `B12`, `B14`, `B15`, `B17`, `B21`**, archivo `fenotipos_TA_homologados_n8_SinOutgroups.csv`) y cruzar los **22 módulos funcionales** contra los **79 compuestos** ($76$ alcaloides `X6T..X81T` + $3$ precursores `X82P..X84P`), se identificaron **101 asociaciones estadísticamente significativas ($\|r\| \ge 0.70, p \le 0.05$)** (`WGCNA_SinB19_03_Top_Asociaciones_Significativas.csv`):

| Módulo (`WGCNA Sin B19`) | Nicho / Significado Biológico | Compuesto Asociado | Tipo | $r$ Pearson ($\log_{10}$) | $r$ Pearson (Lineal) | $p$-valor |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| **`MEred`** ($1.427$ genes) | Clado *Datura* (`+`) vs. Tierras Bajas (`-`) | **`X70T`** | Alcaloide TA | **`-0.9732`** | `-0.8205` | **`0.00005`** |
| **`MEpurple`** ($736$ genes) | **Módulo maestro Alta Montaña (`B10, B11, B15`)** | **`X52T`** | Alcaloide TA | **`-0.9711`** | `-0.8728` | **`0.00006`** |
| **`MEdarkred`** ($82$ genes) | Contraste intra-Tierras Bajas (`B14`) | **`X54T`** | Alcaloide TA | **`+0.9539`** | `+0.5691` | **`0.00024`** |
| **`MEroyalblue`** ($86$ genes) | Específico *Datura* (`B17`) / Tierras Bajas | **`X51T`** / **`X46T`** | Alcaloide TA | **`+0.9398`** / **`-0.8928`** | `+0.8134` / **`-0.9692`** | **`0.00052`** |
| **`MElightcyan`** ($187$ genes) | Contraste `B15` vs. `B11`/`B17` | **`X35T`** | Alcaloide TA | **`-0.9356`** | `-0.7377` | **`0.00064`** |
| **`MEblue`** ($2.345$ genes) | **Módulo maestro Tierras Bajas (`B9, B12, B14, B21`)** | **`X70T`** | Alcaloide TA | **`+0.9213`** | `+0.7807` | **`0.00115`** |
| **`MEgrey60`** ($183$ genes) | Específico `B15` (*B. vulcanicola*) | **`X58T`** / **`X44T`** | Alcaloide TA | **`+0.9193`** / **`+0.8640`** | `+0.6835` / `+0.5939` | **`0.00123`** |
| **`MEmagenta`** ($806$ genes) | **Específico *Brugmansia* Tierras Bajas** | **`X70T`** / **`X60T`** / **`X14T`** | Alcaloide TA | **`+0.9146`** / **`+0.8763`** / **`-0.8842`** | **`+0.9271`** / `+0.7371` / `-0.7818` | **`0.00146`** |
| **`MEmagenta`** ($806$ genes) | **Específico *Brugmansia* Tierras Bajas** | **`X82P`** | **Precursor** | **`+0.8758`** | **`+0.9546`** | **`0.00435`** |
| **`MEmidnightblue`** ($469$ genes) | Alta Montaña (`B15, B10, B11`) + `B14` | **`X17T`** | Alcaloide TA | **`-0.9072`** | `-0.6879` | **`0.00186`** |
| **`MElightgreen`** ($164$ genes) | Específico `B17` / `B15` / `B14` | **`X44T`** | Alcaloide TA | **`+0.9068`** | `+0.8665` | **`0.00189`** |
| **`MEdarkgreen`** ($81$ genes) | Específico `B21` (*Brugmansia aurea*) | **`X22T`** / **`X15T`** (Clúster `X71T`) | Alcaloide TA | **`+0.9041`** / **`-0.8884`** | `+0.3162` / **`-0.9744`** | **`0.00205`** |
| **`MEturquoise`** ($17.694$ genes) | Divergencia *Brugmansia* (`+`) vs. *Datura* (`-`) | **`X70T`** | Alcaloide TA | **`+0.8893`** | `+0.7251` | **`0.00312`** |
| **`MEbrown`** ($1.736$ genes) | Clado *Datura* (`+`) vs. Alta Montaña (`-`) | **`X52T`** | Alcaloide TA | **`+0.8605`** | `+0.8014` | **`0.00609`** |

### Síntesis de Patrones Metabólico-Transcripcionales Validados
1. **Eje `X70T` (Diferenciación Tierras Bajas vs. *Datura*)**: Presenta correlación fuertemente positiva con los tres módulos de *Brugmansia* de Tierras Bajas (`MEblue` $r = +0.9213$, `MEmagenta` $r = +0.9146$, `MEturquoise` $r = +0.8893$) y correlación fuertemente negativa con el módulo opuesto `MEred` ($r = -0.9732, p = 0.00005$).
2. **Eje `X52T` (Diferenciación Alta Montaña vs. *Datura*)**: Presenta correlación fuertemente negativa con el módulo maestro de Alta Montaña `MEpurple` ($r = -0.9711, p = 0.00006$) y positiva con el módulo opuesto `MEbrown` ($r = +0.8605, p = 0.00609$).
3. **Acoplamiento con el Precursor `X82P`**: El módulo `MEmagenta` (*Brugmansia* de Tierras Bajas) se asocia estrechamente con la acumulación del precursor `X82P` ($r = +0.9546$ lineal, $r = +0.8758$ en $\log_{10}$, $p = 0.00435$) junto con los alcaloides `X70T` y `X60T`.

---

## 7. Intersección Expresión Diferencial (Alta Montaña vs. Tierras Bajas) + `WGCNA` (*Hub Genes* Candidatos)

Siguiendo el criterio híbrido de selección de candidatos (Módulo clave + alto `kME` + alto `GS` + expresión diferencial significativa entre Subclado A Alta Montaña `B10, B11, B15` y Subclado B Tierras Bajas `B9, B12, B14, B21`), el script `DESeq_y_HubGenes_Candidatos.py` identificó:
- **3.096 transcritos diferenciales** ($\|\log_2\text{FC}\| \ge 1.5, p \le 0.01$): **644 sobre-expresados en Alta Montaña** y **2.452 sobre-expresados en Tierras Bajas**.
- **Polaridad direccional del 100% en los módulos ecológicos**:
  - Módulos con **100% de sus DEGs sobre-expresados en Alta Montaña (Subclado A)**: `purple` ($269/269$, $\overline{\log_2\text{FC}} = +3.78$), `red` ($117/117$, $\overline{\log_2\text{FC}} = +4.12$), `yellow` ($77/77$, $\overline{\log_2\text{FC}} = +2.53$), `greenyellow` ($29/29$, $\overline{\log_2\text{FC}} = +2.94$).
  - Módulos con **100% de sus DEGs sobre-expresados en Tierras Bajas (Subclado B)**: `blue` ($304/304$, $\overline{\log_2\text{FC}} = -3.46$), `brown` ($207/207$, $\overline{\log_2\text{FC}} = -3.66$), `magenta` ($92/92$, $\overline{\log_2\text{FC}} = -3.44$), `salmon` ($47/47$, $\overline{\log_2\text{FC}} = -3.16$).
- **1.440 *Hub Genes* Candidatos de Alta Prioridad** exportados en `WGCNA_SinB19_04_Top_HubGenes_Candidatos_DEGs_TA.csv` que cumplen simultáneamente `kME >= 0.85`, $\|\log_2\text{FC}\| \ge 1.5$ ($p \le 0.01, \text{FDR} < 0.05$) y $\|\text{GS}\| \ge 0.75$ con los alcaloides clave (`X70T`, `X52T`, `X82P`, `X71T`).
