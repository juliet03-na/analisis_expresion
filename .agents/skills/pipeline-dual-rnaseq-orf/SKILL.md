---
name: pipeline-dual-rnaseq-orf
description: >-
  Guía técnica y generador de comandos/scripts para los Módulos 1 y 2 de la tesis
  hologenómica en el clado Datura: separación taxonómica con FastqPuri (filtro de Bloom),
  co-ensamblaje de novo en rnaSPAdes, validación estructural (CD-HIT-EST, TransRate,
  Bowtie2, SAMtools), predicción diferenciada de ORFs (TransDecoder para planta eucariota
  vs. Prodigal para bacteria procariota) y anotación funcional contra los 35 genes de la
  ruta de alcaloides tropánicos (TA) con DIAMOND. Úsese cuando el usuario solicite scripts
  HPC/Slurm, comandos o explicaciones de los Módulos 1 o 2, o de las Ramas A (lecturas)
  y B (contigs).
---

# Pipeline Dual RNA-seq: Separación Taxonómica, Ensamblaje y Predicción de ORFs (Módulos 1 y 2)

Esta *skill* implementa los protocolos exactos de [[03_MODULO1_SEPARACION_Y_ENSAMBLAJE]] y [[04_MODULO2_PREDICCION_ORF_Y_ANOTACION]].

## 1. Módulo 1: Preprocesamiento y Separación Taxonómica (`FastqPuri`)

1. **Diagnóstico y Limpieza**:
   - `FastQC v0.11.9` y `MultiQC v1.25`.
   - `Trimmomatic v0.39`:
     ```bash
     trimmomatic PE -threads 16 raw_R1.fq.gz raw_R2.fq.gz \
       clean_R1.fq.gz unpaired_R1.fq.gz clean_R2.fq.gz unpaired_R2.fq.gz \
       ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:140
     ```
   - Remoción de ARNr residual con `RNA-QC-Chain`.
2. **Filtro de Bloom (`FastqPuri v1.0.1`)**:
   - Construir el índice probabilístico (`makeBloom`) sobre el **Panel Solanaceae** (*Brugmansia arborea* `PRJNA765960`, *Datura stramonium* `GCA_018107945.1` + 2° ensamblaje, *Atropa belladonna* `PRJNA766188`, *Solanum lycopersicum* `GCF_036512215.1`) y **Virus RefSeq**.
   - Ejecutar `trimFilterPE` con opción `BLOOM` para particionar en:
     - **Set Hospedero** (alinea con Solanaceae).
     - **Set Microorganismos Asociados** (no alinea ni con plantas ni con virus).

## 2. Módulo 1: Co-ensamblaje *de novo* y Validación Estructural

1. **Co-ensamblaje en `rnaSPAdes` (`v3.15.4` / `v3.15.5`)**:
   - Ejecutar de manera independiente para Set Hospedero y Set Microbiano:
     ```bash
     rnaspades.py --rna -t 64 -m 230 -k 49,73 \
       -1 host_R1.fq.gz -2 host_R2.fq.gz -o rnaspades_host_out
     ```
   - El ensamblaje base vegetal produce **123.099 transcritos** ($165\text{ MB}$).
2. **Validación Estructural en 4 Etapas**:
   - **Reducción de redundancia (`CD-HIT-EST v4.8.1`)**:
     ```bash
     cd-hit-est -i transcripts.fasta -o transcripts_clustered.fasta -c 0.95 -aS 0.80 -M 16000 -T 8
     ```
   - **Contigüidad (`TransRate v1.0.3`)**: `transrate --assembly transcripts_clustered.fasta`
   - **Representatividad (`Bowtie2 v2.4.5`)**: Verificar *overall alignment rate* $> 95\%$.
   - **Cobertura (`SAMtools v1.15`)**: `samtools stats` y `samtools coverage`.

## 3. Módulo 2: Predicción Diferenciada de ORFs (Eucariota vs. Procariota)

> [!IMPORTANT]
> Nunca uses `TransDecoder` sobre el ensamblaje bacteriano. `TransDecoder` asume ARNm monocistrónicos e inicio exclusivo en `ATG`, ignorando operones policistrónicos y sitios *Shine-Dalgarno*.

- **Transcriptoma Vegetal (Eucariota) $\rightarrow$ `TransDecoder`**:
  ```bash
  TransDecoder.LongOrfs -t transcripts_clustered.fasta -m 100
  TransDecoder.Predict -t transcripts_clustered.fasta
  ```
  *(Genera archivos `.pep`, `.cds`, `.gff3` y `.bed`).*
- **Transcriptoma Bacteriano (Procariota) $\rightarrow$ `Prodigal`**:
  ```bash
  prodigal -i microbial_transcripts_clustered.fasta -a microbial_proteins.faa -d microbial_cds.fna -p meta
  ```

## 4. Módulo 2: Catálogo de 35 Genes TA y Anotación con `DIAMOND`

- **35 identificadores génicos nominales de la base TA**:
  `ARG2P`, `PMT`, `OTC/ArgF`, `ASS`, `ASL`, `ADC2`, `AIH`, `CPA`, `Speb`, `MPO`, `PYKS`, `CYP82M3`, `TRI`, `LS`, `CYP80F1`, `HDH`, `H6Hhi`, `H6Hdio`, `H6H`, `AT4`, `fldH`, `PPAR`, `UGT1`, `GDH`, `aspC`, `HPA3`, `argB`, `argC`, `argD`, `argJ`, `argE`, `DELTA-OAT`, `proB`, `proA`, `PYCR1`.
- **Construcción de la base**: Secuencias GenBank (*B. arborea* + Naik et al., 2018) $\rightarrow$ `PSI-BLAST` (`max_target_seqs = 500`, `word_size = 3`, `BLOSUM62`, `gap_open = 11`, `gap_extend = 1`) $\rightarrow$ `MAFFT v7.526`.
- **Mapeo en `DIAMOND v2.1.16`**:
  - **Rama A (Lecturas crudas filtradas)**: `diamond blastx -d TA_db.dmnd -q reads.fq -e 0.001 -p 32 --outfmt 6 qseqid sseqid pident length bitscore evalue scovhsp --al aligned_reads.fq`
  - **Rama B (Proteínas predichas `.pep`/`.faa`)**: `diamond blastp -d TA_db.dmnd -q predicted.pep -e 0.001 -p 32 --outfmt 6 qseqid sseqid pident length bitscore evalue scovhsp --al aligned_proteins.faa`
