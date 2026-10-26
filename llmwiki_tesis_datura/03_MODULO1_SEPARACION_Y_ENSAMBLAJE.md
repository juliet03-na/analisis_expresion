---
id_nodo: 03_MODULO1_SEPARACION_Y_ENSAMBLAJE
tipo: modulo_operativo
modulo: 1
nodos_relacionados:
  - "[[00_INDEX]]"
  - "[[01_MARCO_Y_OBJETIVOS]]"
  - "[[02_MUESTRAS_Y_WETLAB]]"
  - "[[04_MODULO2_PREDICCION_ORF_Y_ANOTACION]]"
---

# `03_MODULO1_SEPARACION_Y_ENSAMBLAJE` — FastqPuri, rnaSPAdes y Validación Estructural

- **ID_Nodo**: `03_MODULO1_SEPARACION_Y_ENSAMBLAJE`
- **Nodos_Relacionados**: [[00_INDEX]], [[01_MARCO_Y_OBJETIVOS]], [[02_MUESTRAS_Y_WETLAB]], [[04_MODULO2_PREDICCION_ORF_Y_ANOTACION]]

---

## 1. Preprocesamiento y Segregación Taxonómica *In Silico*

- **Diagnóstico Inicial**: `FastQC v0.11.9` (Andrews, 2010) y `MultiQC v1.25` (Ewels et al., 2016).
- **Limpieza y Recorte**: `Trimmomatic v0.39` (Bolger et al., 2014) con parámetros:
  - `ILLUMINACLIP: adapters.fa:2:30:10`
  - `LEADING: 3`
  - `TRAILING: 3`
  - `SLIDINGWINDOW: 4:15`
  - `MINLEN: 140`
- **Remoción de ARNr Residual**: `RNA-QC-Chain` (Zhou et al., 2018).
- **Separación Planta vs. Microbiota (Filtro de Bloom)**: `FastqPuri v1.0.1` (Pérez-Rubio et al., 2019; Engelmann et al., 2019).
  - **Creación de índice probabilístico**: Función `makeBloom` sobre los genomas de referencia de Solanaceae y Virus RefSeq ([[02_MUESTRAS_Y_WETLAB]]).
  - **Clasificación de lecturas**: Función `trimFilter(PE)` con opción `BLOOM`.
  - **Partición resultante**:
    - **Set Hospedero (*Datura*/*Brugmansia*)**: Lecturas que alinean contra el panel Solanaceae.
    - **Set Microorganismos Asociados (Bacterias, Hongos y otros endófitos)**: Lecturas sin alineamiento contra plantas ni virus.

---

## 2. Co-ensamblaje *de novo* y Control de Calidad Estructural

- **Ensamblador**: `rnaSPAdes` (integrado en `SPAdes v3.15.4` / `v3.15.5`, Bushmanova et al., 2019) ejecutado de forma paralela e independiente sobre el conjunto de lecturas del hospedero y el conjunto microbiano.
- **Parámetros de ejecución en HPC**: Modo `--rna`, `-t 64` hilos, `-m 230` GB de RAM, y *k*-meros iterativos $K=49, 73$ ajustados a lecturas de $150\text{ pb}$ (generando **135.008 contigs crudos** en el ensamblaje inicial, que tras `CD-HIT-EST` se consolidan en un archivo de $165\text{ MB}$ en `FASTA` con **123.099 transcritos**, reducción de redundancia del $8,8\%$ y tasa de alineamiento global de **98,3%** en `Bowtie2`).

### Etapas de Validación Estructural del Ensamblaje

| Etapa de Validación | Software (Versión) | Parámetros y Criterios Técnicos | Propósito Metodológico |
| :--- | :--- | :--- | :--- |
| **1. Reducción de Redundancia** | `CD-HIT-EST` (`v4.8.1`) | `-c 0.95 -aS 0.80 -M 16000 -T 8`. Ordena secuencias de mayor a menor longitud; la más larga es la representativa. | Colapsar errores de secuenciación, isoformas con variación mínima y polimorfismos ($\ge 95\%$ identidad) asegurando que el fragmento corto quede cubierto en $\ge 80\%$ (`transcripts_clustered.fasta`). |
| **2. Métricas de Contigüidad** | `TransRate` (`v1.0.3`) | `--assembly` (sin umbral estricto). | Calcular número total de transcritos, tamaño promedio de contigs, longitud máxima y estadístico N50. |
| **3. Representatividad** | `Bowtie2` (`v2.4.5`) | Mapeo de lecturas crudas contra el co-ensamblaje depurado; umbral $> 95\%$ de *overall alignment rate*. | Verificar que la diversidad transcriptómica de cada muestra no se perdió durante el clustering con `CD-HIT`. |
| **4. Mapeo Espacial** | `SAMtools` (`v1.15`) | `samtools stats` y `samtools coverage`. | Evaluar porcentaje de identidad, cobertura de bases y profundidad media para soportar expresión diferencial. |
