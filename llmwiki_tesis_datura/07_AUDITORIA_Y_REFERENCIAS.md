---
id_nodo: 07_AUDITORIA_Y_REFERENCIAS
tipo: nodo_transversal
nodos_relacionados:
  - "[[00_INDEX]]"
  - "[[02_MUESTRAS_Y_WETLAB]]"
  - "[[03_MODULO1_SEPARACION_Y_ENSAMBLAJE]]"
  - "[[06_MODULO4_METABOLOMICA_Y_COEXPRESION]]"
---

# `07_AUDITORIA_Y_REFERENCIAS` — Inconsistencias de Datos (`B27`/`B20`, $n=10$) y Citas

- **ID_Nodo**: `07_AUDITORIA_Y_REFERENCIAS`
- **Nodos_Relacionados**: [[00_INDEX]], [[02_MUESTRAS_Y_WETLAB]], [[03_MODULO1_SEPARACION_Y_ENSAMBLAJE]], [[06_MODULO4_METABOLOMICA_Y_COEXPRESION]]

---

## 1. Problemas Críticos Detectados en las Matrices de Datos Experimentales

> [!CAUTION]
> Estos tres puntos deben verificarse de forma obligatoria en cualquier script de R/Python antes de ejecutar correlaciones multi-ómicas o redes de co-expresión.

1. **Desajuste de IDs y orden de filas entre Expresión y Metabolómica (`B27` vs. `B20`, $n=8/9$)**:
   - La matriz de expresión génica (`matriz_vst_filtrada_para_WGCNA.csv`) tiene las muestras en el orden `B9, B10, B11, B12, B14, B15, B17, B19, B21, B27` (incluye **`B27`** *Datura wrightii*), mientras que la matriz metabolómica (`Metaboloma_WGNA_TESIS.xlsx`) viene en el orden `B21, B9, B12, B14, B10, B11, B15, B19, B17, B20` (incluye **`B20`** *Ipomoea purpurea* y no contiene `B27`).
   - Al homologar las filas del clado *Datura*/*Brugmansia* excluyendo los *outgroups* (`B19` y `B20`) y `B27`, quedan **8 muestras perfectamente alineadas** (`B9, B10, B11, B12, B14, B15, B17, B21`, exportadas en `fenotipos_TA_homologados_n8_SinOutgroups.csv`).
2. **Corrupción de puntos decimales entre `Metaboloma_WGNA_TESIS.xlsx` y `fenotipos_TA_limpios_para_WGCNA.csv` (Hallazgo Forense Crítico)**:
   - El archivo `Metaboloma_WGNA_TESIS.xlsx` contiene **79 compuestos** sin ceros ni celdas vacías (**76 alcaloides tropánicos** `6T` a `81T` y **3 precursores** `82P`, `83P`, `84P`).
   - En el XML interno de `Metaboloma_WGNA_TESIS.xlsx`, **$339$ de las $790$ celdas ($42,9\%$)** están almacenadas como **texto (`t="s"`)** con punto decimal (`151` celdas con 1 decimal `.X` y `188` celdas con 2 decimales `.XX`, ej. `33728376.94`), mientras que **$451$ celdas ($57,1\%$)** son enteros de 9–11 dígitos y hay **$0$ celdas con 3 decimales (`.XXX`)**. Esto ocurrió porque Excel en configuración regional en español interpretó los puntos seguidos de 3 dígitos (`.XXX`) como separadores de miles, y dejó como texto las celdas con 1 o 2 decimales (`.X` y `.XX`).
   - Al generar `fenotipos_TA_limpios_para_WGCNA.csv`, se eliminaron los puntos `.` de las $339$ celdas de texto (`33728376.94` $\rightarrow$ `3372837694`), multiplicando artificialmente las celdas de 2 decimales por **$100\times$** y las de 1 decimal por **$10\times$** dentro de la misma columna.
   - **Impacto cuantificado**: Esta pérdida del punto decimal distorsionó **$61$ de las $79$ columnas ($77,2\%$)** con $r < 0.95$ respecto al Excel ($42$ columnas con $r < 0.80$; $25$ columnas con $r < 0.50$), llegando a **invertir el signo de la correlación** en `X14T` ($r = -0.5819$), `X64T` ($r = -0.3823$), `X61T` ($r = -0.2786$), `X15T` ($r = -0.2753$), `X79T` ($r = -0.1938$), `X72T` ($r = -0.1455$) y `X17T` ($r = -0.0724$).
   - **Corrección**: Se reconstruyó la matriz numérica exacta desde el XML del Excel en **`fenotipos_TA_corregidos_para_WGCNA.csv`** y **`fenotipos_TA_homologados_n8_SinOutgroups.csv`**. Adicionalmente, si se dispone del archivo crudo de espectrometría de masas previo a su apertura en Excel en español, pueden restaurarse también los decimales `.XXX` de las $451$ celdas enteras.
3. **Limitación de potencia estadística (*Small $n$, Large $p$*)**:
   - El ensayo cuenta con **$n=10$ muestras efectivas** en el cruce, por debajo del mínimo de $15$ muestras recomendado por el manual de `WGCNA`.
   - Esto genera falsos positivos e inestabilidad en `pickSoftThreshold()`, haciendo indispensable respaldar los hallazgos de `WGCNA` con el cruce de expresión diferencial (`DESeq2`) y métodos regularizados (`sPLS`/`DIABLO`, `MOFA+`, Louvain o Mantel).
4. **Resolución del conteo exacto del embudo (`54.363` pre-MAD vs. `36.647` post-MAD) y sesgo del *outgroup* `B19`**:
   - A partir de `Normalizacion3.R` y los archivos `quant.sf` de `plantas/plantas/*_quant/`, quedó demostrado que **54.363** es el número exacto de transcritos **antes del filtro MAD** (tras `rowSums >= 25`, `max <= 50.000` y `max <= 30.000` sobre los $123.099$ contigs post-`CD-HIT-EST` provenientes de $135.008$ contigs crudos), mientras que **36.647** es el número exacto de transcritos **después de aplicar `mads > quantile(mads, 0.25)`** en `matriz_vst_filtrada_para_WGCNA.csv` ($10\text{ muestras} \times 36.647\text{ transcritos}$, dividida en 4 bloques con `maxBlockSize = 10000`).
   - En `Salmon v1.10.3`, la muestra `B19` (*Atropa belladonna*) presenta solo **127.327 lecturas procesadas** ($122.351$ mapeadas, $96,09\%$) frente a $4,07\text{–}22,01$ millones en *Brugmansia*.
   - **Demostración cuantitativa del sesgo por `B19`**: Al mantener `B19`, su baja cobertura y divergencia taxonómica ($r = 0.19\text{–}0.33$) fuerzan la aparición de un módulo artificial (`brown`, $3.222$ transcritos con $\text{ME} = -2.84$ en `B19`). Al excluir `B19` ($n=9$, `power = 12`, `pamRespectsDendro = FALSE`), **$614$ de esos transcritos ($19,06\%$)** se descartan como falsos positivos (`grey`), el porcentaje global de `grey` baja a **`8,10%`** y emergen **$22$ módulos funcionales** con separación limpia entre Alta Montaña (`MEpurple`), Tierras Bajas (`MEblue`, `MEmagenta`) y *Datura* (`MEbrown`, `MEred`, `MEcyan`, `MEtan`).
5. **Blindaje técnico en `WGCNA1_Data.R` (`stats::cor` vs. `WGCNA::cor`)**:
   - En R moderno, `blockwiseModules()` falla con `unused argument (weights.x = NULL)` si otro paquete o R base antepone `stats::cor`. Se incorporó `cor <- WGCNA::cor` antes de `blockwiseModules()` y se parametrizaron explícitamente `deepSplit = 2.5`, `reassignThreshold = 0`, `REMOVER_OUTGROUP_B19 = TRUE` y `POTENCIA_FIJA = 12`.

---

## 2. Inconsistencias Redaccionales en el Manuscrito (`NewMethods.docx`)

| Ítem Auditado | Inconsistencia Detectada en `NewMethods.docx` | Acción Correctiva Requerida |
| :--- | :--- | :--- |
| **Versiones y *k*-meros de `rnaSPAdes`** | La sección 2.4 reporta `v3.15.5` con *k*-meros automáticos, mientras que la sección 2.7.1 reporta `v3.15.4` con *k*-meros empíricos fijos $K=49$ y $73$. | Unificar la versión exacta ejecutada en el HPC y dejar explícitos los *k*-meros $K=49, 73$ en ambas secciones. |
| **Número de librerías secuenciadas** | La Tabla 3 lista $12$ muestras, pero la sección 2.7.1 y los apuntes de filtrado (umbral $7/10$ muestras y $n=10$) indican el procesamiento de $10$ librerías. | Aclarar en el texto en qué paso se separan los *outgroups* `B19` (*Atropa belladonna*) y `B20` (*Ipomoea purpurea*). |
| **Conteo de genes en base de datos TA** | El texto declara $32$ secuencias proteicas, pero enumera $35$ identificadores génicos nominales. | Sincronizar el número declarado con el catálogo real de 35 genes nominales ([[04_MODULO2_PREDICCION_ORF_Y_ANOTACION]]). |
| **Error de herramienta en búsqueda de homólogos** | La sección 2.5.1 cita `BLASTn` (nucleótidos) pero describe parámetros exclusivos de proteínas (`PSI-BLAST`, `word size = 3`, matriz `BLOSUM62`, `gap open = 11`, `gap extend = 1`). | Corregir en la sección 2.5.1 el término `BLASTn` por `BLASTp`/`PSI-BLAST`. |
| **Omisión de `Prodigal`, `TransDecoder`, `DESeq2` y `WGCNA`** | El archivo `NewMethods.docx` se detiene en la cuantificación con `Salmon`. | Redactar e incorporar formalmente las secciones de predicción de ORFs (`TransDecoder`/`Prodigal`), filtrado y normalización (`tximport`/`DESeq2`/VST), curaduría metabolómica (de 82 a 33 moléculas) y redes de co-expresión (`WGCNA`) documentadas en el cuaderno de laboratorio. |

---

## 3. Estado de Referencias Bibliográficas

### Completas en `NewMethods.docx`
- Altschul et al. (1990, 1997)
- Andrews (2010)
- Benson et al. (2018)
- Bolger et al. (2014)
- Buchfink et al. (2021)
- Bushmanova et al. (2019)
- Capella-Gutiérrez et al. (2009)
- Conesa et al. (2016)
- Danecek et al. (2021)
- Eddy (2011)
- Ewels et al. (2016)
- Fu et al. (2012)
- Katoh & Standley (2013)
- Langmead & Salzberg (2012)
- Pérez-Rubio et al. (2019)
- Smith-Unna et al. (2016)

### Incompletas con *placeholders* en `NewMethods.docx`
- Bedewitz et al. (2018), Chen et al. (2023), Cokelaer et al. (2023), Gao et al. (2025), Gruening et al. (2023), Guo et al. (2019), Hacquard et al. (2024), Hernández-Duarte et al. (2025), Jiang et al. (2024), Kavanagh et al. (2008), Li et al. (2016), Malik et al. (2023), Naik et al. (2018), Schlesinger et al. (2021), Singh et al. (2024), Trapnell et al. (2012), Wang et al. (2025), Zhou et al. (2018).

### Citadas en planteamiento/apuntes pendientes de agregar a la bibliografía
- Fitzpatrick et al. (2020), Xu et al. (2021), Jain et al. (2024), Chen et al. (2022), Taj et al. (2023), Shakya et al. (2019), Di Leo et al. (2025), Engelmann et al. (2019) y las referencias de software para `Salmon` (Patro et al., 2017), `tximport` (Soneson et al., 2015), `DESeq2` (Love et al., 2014), `WGCNA` (Langfelder & Horvath, 2008), `TransDecoder` (Haas et al., 2013) y `Prodigal` (Hyatt et al., 2010).
