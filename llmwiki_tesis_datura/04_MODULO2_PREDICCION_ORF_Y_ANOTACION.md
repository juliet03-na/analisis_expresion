---
id_nodo: 04_MODULO2_PREDICCION_ORF_Y_ANOTACION
tipo: modulo_operativo
modulo: 2
nodos_relacionados:
  - "[[00_INDEX]]"
  - "[[03_MODULO1_SEPARACION_Y_ENSAMBLAJE]]"
  - "[[06_MODULO4_METABOLOMICA_Y_COEXPRESION]]"
---

# `04_MODULO2_PREDICCION_ORF_Y_ANOTACION` — TransDecoder vs. Prodigal, Base TA y DIAMOND

- **ID_Nodo**: `04_MODULO2_PREDICCION_ORF_Y_ANOTACION`
- **Nodos_Relacionados**: [[00_INDEX]], [[03_MODULO1_SEPARACION_Y_ENSAMBLAJE]], [[06_MODULO4_METABOLOMICA_Y_COEXPRESION]]

---

## 1. Predicción de Regiones Codificantes (Diferenciación Eucariota vs. Procariota)

- **Justificación en el flujo**: Dado que `WGCNA` y `DESeq2` entregan listados de transcritos anónimos altamente correlacionados con los alcaloides, es obligatorio predecir sus péptidos codificados para anotarlos funcionalmente con `DIAMOND` y confirmar a qué enzima de la ruta pertenecen.

| Atributo | Transcriptoma Vegetal (Eucariota) | Transcriptoma Bacteriano (Procariota) |
| :--- | :--- | :--- |
| **Software Asignado** | `TransDecoder` | `Prodigal` |
| **Justificación Biológica** | Diseñado para transcritos eucariotas monocistrónicos. | `TransDecoder` falla en procariotas porque ignora estructuras de operones (ARNm policistrónicos), asume inicio exclusivo en `ATG` y no reconoce sitios de unión ribosomal *Shine-Dalgarno*. |
| **Ejecución y Parámetros** | 1) `TransDecoder.LongOrfs -t transcripts_clustered.fasta -m 100` (evalúa los 6 marcos de lectura posibles con longitud mínima de 100 aminoácidos).<br/>2) `TransDecoder.Predict` (aplica modelo de cadenas de Markov de 5° orden). | Ejecución estándar de predicción de genes procariotas sobre contigs bacterianos. |
| **Archivos de Salida** | `.pep` (proteínas predichas completas y parciales para `DIAMOND`), `.cds` (secuencias codificantes nucleotídicas), `.gff3` (coordenadas estructurales) y `.bed` (pistas de visualización). | Archivos de proteínas predichas (`.faa`/`.pep`) y secuencias codificantes para anotación funcional. |

---

## 2. Base de Datos de Referencia de Alcaloides Tropánicos (TA) y Anotación

- **Rutas metabólicas cubiertas**:
  - Ciclo del citrato
  - Metabolismo de arginina y prolina
  - Biosíntesis de poliaminas
  - Formación del anillo tropano
  - Síntesis de fenilalanina, tirosina, triptófano
  - Ácido trópico, hiosciamina y escopolamina
- **Catálogo de genes incluidos (multifasta indexado — 35 identificadores nominales)**:
  `ARG2P`, `PMT`, `OTC/ArgF`, `ASS`, `ASL`, `ADC2`, `AIH`, `CPA`, `Speb`, `MPO`, `PYKS`, `CYP82M3`, `TRI`, `LS`, `CYP80F1`, `HDH`, `H6Hhi`, `H6Hdio`, `H6H`, `AT4`, `fldH`, `PPAR`, `UGT1`, `GDH`, `aspC`, `HPA3`, `argB`, `argC`, `argD`, `argJ`, `argE`, `DELTA-OAT`, `proB`, `proA`, `PYCR1`.
- **Construcción del compendio**:
  Secuencias de GenBank (Benson et al., 2018) basadas en *B. arborea* y Naik et al. (2018) $\rightarrow$ búsqueda `PSI-BLAST` (`max_target_seqs = 500`, `word_size = 3`, `BLOSUM62`, `gap_open = 11`, `gap_extend = 1`) $\rightarrow$ alineamiento múltiple con `MAFFT v7.526` (Katoh & Standley, 2013).
- **Anotación Funcional y Mapeo Dirigido (`DIAMOND`)**:
  `DIAMOND v2.1.16` (Buchfink et al., 2021) ejecutado tanto sobre lecturas crudas (`blastx`, Rama A) como sobre proteínas predichas por `TransDecoder`/`Prodigal` (`blastp`, Rama B), con:
  - `E-value <= 0.001`
  - `--threads 32`
  - `--outfmt 6` (porcentaje de identidad, `bitscore`, `scovhsp`)
  - `--al`
