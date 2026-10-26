---
name: auditoria-matrices-y-manuscrito
description: >-
  Auditor de calidad de datos experimentales, matrices de R y redacción del manuscrito
  (NewMethods.docx) de la tesis en el clado Datura. Úsese cuando el usuario pida revisar
  inconsistencias entre tablas de expresión y fenotipos (caso B27 vs. B20, columnas con
  ceros, n=10), auditar o redactar secciones de metodología del documento Word
  (NewMethods.docx), o completar y verificar referencias bibliográficas.
---

# Auditoría de Matrices Experimentales, Manuscrito (`NewMethods.docx`) y Referencias

Esta *skill* aplica las reglas de control de calidad documentadas en [[07_AUDITORIA_Y_REFERENCIAS]].

## 1. Checklist de Auditoría sobre Matrices de Datos

1. **Desajuste `B27` vs. `B20` y Orden de Filas ($n=8$ sin *outgroups*)**:
   - En expresión (`matriz_vst_filtrada_para_WGCNA.csv`) las filas están en orden `B9, B10, B11, B12, B14, B15, B17, B19, B21, B27`, mientras que en metabolómica (`Metaboloma_WGNA_TESIS.xlsx`) están en orden `B21, B9, B12, B14, B10, B11, B15, B19, B17, B20`.
   - Al homologar las muestras del clado *Datura*/*Brugmansia* sin *outgroups* (`B19`/`B20`) ni `B27`, quedan **8 muestras alineadas**: `B9, B10, B11, B12, B14, B15, B17, B21` (`fenotipos_TA_homologados_n8_SinOutgroups.csv`).
2. **Corrupción de Puntos Decimales en `fenotipos_TA_limpios_para_WGCNA.csv`**:
   - En `Metaboloma_WGNA_TESIS.xlsx` ($10\text{ muestras} \times 79\text{ compuestos}$: $76$ alcaloides `6T..81T` y $3$ precursores `82P..84P`), $339$ de las $790$ celdas tenían 1 o 2 decimales (`.X`/`.XX`) y quedaron como texto en Excel en español (mientras que las celdas con `.XXX` fueron tomadas como separador de miles).
   - Al crear `fenotipos_TA_limpios_para_WGCNA.csv` se borraron los puntos `.`, multiplicando esas $339$ celdas por $10\times$ o $100\times$ y distorsionando $61$ de las $79$ columnas (invirtiendo incluso el signo en `X14T`, `X64T`, `X61T`, `X15T`).
   - Exigir siempre el uso de **`fenotipos_TA_corregidos_para_WGCNA.csv`** o **`fenotipos_TA_homologados_n8_SinOutgroups.csv`**.
3. **Advertencia de Potencia Estadística ($n=10$)**:
   - Documentar en qué paso se separan los *outgroups* (`B19` *Atropa belladonna* y `B20` *Ipomoea purpurea*) pasando de 12 muestras en inventario a $n=10$ muestras efectivas.
   - Respaldar siempre `WGCNA` ($n=10 < 15$ recomendados) con `DESeq2` y métodos penalizados/multivariados (`sPLS`/`DIABLO`, `MOFA+ v2`, Louvain, Mantel/CIA).

## 2. Correcciones Obligatorias al Redactar o Editar `NewMethods.docx`

Al generar texto académico para la sección de Materiales y Métodos, corrige sistemáticamente estas 5 inconsistencias del borrador:
1. **Unificar versión y *k*-meros de `rnaSPAdes`**: Armonizar la sección 2.4 (`v3.15.5`, *k*-meros automáticos) con la sección 2.7.1 (`v3.15.4`, *k*-meros empíricos fijos $K=49, 73$).
2. **Trazabilidad de 12 vs. 10 librerías**: Explicar cómo las 12 muestras de la Tabla 3 se reducen a 10 librerías en el filtro de expresión ($7/10$ muestras) al separar los *outgroups* (`B19`/`B20`).
3. **Sincronizar conteo de la base TA**: Corregir "32 secuencias proteicas" para que coincida con los **35 identificadores génicos nominales** listados en [[04_MODULO2_PREDICCION_ORF_Y_ANOTACION]].
4. **Corregir `BLASTn` por `BLASTp`/`PSI-BLAST`**: En la sección 2.5.1, sustituir `BLASTn` por `PSI-BLAST`/`BLASTp`, dado que los parámetros empleados (`word_size = 3`, matriz `BLOSUM62`, `gap_open = 11`, `gap_extend = 1`) corresponden exclusivamente a secuencias proteicas.
5. **Incorporar las secciones omitidas post-Salmon**: Redactar formalmente:
   - Predicción de ORFs diferenciada por reino (`TransDecoder` para planta eucariota vs. `Prodigal` para microorganismos procariotas).
   - Importación con `tximport` (corrección por longitud efectiva $L_{ij}$), embudo estadístico ($123.099 \rightarrow 54.363$ transcritos tras VST y filtro MAD del 75% superior) y expresión diferencial en `DESeq2`.
   - Curaduría metabolómica en 4 etapas ($82 \rightarrow 79 \rightarrow 73 \rightarrow 33 \rightarrow 15\text{–}25$ moléculas).
   - Construcción de redes `WGCNA`, correlación cruzada de *Module Eigengenes* planta-endófito y métodos multivariados para *"Small $n$, Large $p$"*.
