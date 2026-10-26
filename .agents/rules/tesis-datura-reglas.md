# Reglas Operativas del Proyecto: Tesis Hologenómica *Datura*

Estas reglas rigen toda generación de código, análisis estadístico, interpretación biológica y redacción académica dentro del repositorio.

## 1. Integridad de Matrices y Muestras
- **Subclado A (Alta montaña / Tierras altas)**: `B10` (*B. sanguinea*), `B11` (*B. sanguinea x vulcanicola*), `B15` (*B. vulcanicola*), `P1` (*B. arborea*).
- **Subclado B (Tierras bajas)**: `B9` (*B. insignis*), `B12` (*B. suaveolens*), `B14` (*B. versicolor*), `B21` (*B. aurea*).
- **Clado *Datura***: `B17` (*D. stramonium*), `B27` (*D. wrightii* — presente en expresión, ausente en tabla fenotípica cruda).
- **Outgroups**: `B19` (*Atropa belladonna*, outgroup clado *Datura*), `B20` (*Ipomoea purpurea*, outgroup Solanaceae — presente en fenotipos, ausente en expresión).
- **Regla de Oro en R y Python**:
  - Antes de ejecutar `WGCNA`, `sPLS`, `MOFA+` o cualquier correlación expresión-metaboloma, verifica la intersección exacta y el orden de filas (`B9, B10, B11, B12, B14, B15, B17, B21` para $n=8$ sin *outgroups*).
  - **Prohibido usar `fenotipos_TA_limpios_para_WGCNA.csv`** (contiene $339$ celdas multiplicadas por $10\times$ o $100\times$ al borrarse el punto decimal respecto a `Metaboloma_WGNA_TESIS.xlsx`). Usa siempre **`fenotipos_TA_corregidos_para_WGCNA.csv`** o **`fenotipos_TA_homologados_n8_SinOutgroups.csv`**.

## 2. Trazabilidad de Software y Parámetros Oficiales
- **Limpieza**: `FastQC v0.11.9`, `MultiQC v1.25`, `Trimmomatic v0.39` (`ILLUMINACLIP: adapters.fa:2:30:10`, `LEADING: 3`, `TRAILING: 3`, `SLIDINGWINDOW: 4:15`, `MINLEN: 140`), `RNA-QC-Chain`.
- **Separación Taxonómica**: `FastqPuri v1.0.1` (`makeBloom` + `trimFilter(PE)` con `BLOOM` contra panel Solanaceae + Virus RefSeq).
- **Co-ensamblaje y Validación**: `rnaSPAdes` (`--rna -t 64 -m 230`, $K=49, 73$), `CD-HIT-EST v4.8.1` (`-c 0.95 -aS 0.80 -M 16000 -T 8`), `TransRate v1.0.3`, `Bowtie2 v2.4.5` ($>95\%$ *overall alignment rate*), `SAMtools v1.15`.
- **Predicción de ORFs y Anotación TA**:
  - Eucariota (Vegetal): `TransDecoder` (`LongOrfs -m 100` + `Predict`).
  - Procariota (Bacteriano): `Prodigal`.
  - Alineamiento y Homología: `PSI-BLAST` (`max_target_seqs = 500`, `word_size = 3`, `BLOSUM62`, `gap_open = 11`, `gap_extend = 1`), `MAFFT v7.526`, `DIAMOND v2.1.16` (`blastx` para Rama A en lecturas; `blastp` para Rama B en proteínas predichas; `E-value <= 0.001 --threads 32 --outfmt 6 --al`).
- **Cuantificación, Embudo y Redes**: `Salmon v1.10.3` (`-l A --validateMappings`), `tximport` (offset $\log(L_{ij})$), `DESeq2` (`vst` modo `blind`, filtro MAD 75% superior $\rightarrow 54.363$ transcritos), `WGCNA` (`TOMType = "signed"`, `corType = "pearson"`, `deepSplit = 2.5`, `mergeCutHeight = 0.25`).

## 3. Mantenimiento del Grafo LLMWiki
- Cualquier actualización metodológica o nuevo resultado experimental debe reflejarse manteniendo enlaces bidireccionales `[[ID_NODO]]` entre los 8 archivos de `llmwiki_tesis_datura/`.
