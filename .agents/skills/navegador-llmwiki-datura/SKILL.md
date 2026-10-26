---
name: navegador-llmwiki-datura
description: >-
  Gestor arquitectónico del grafo de conocimiento LLMWiki (llmwiki_tesis_datura/) de 8
  nodos modulares. Úsese cuando el usuario pida actualizar, consultar, enlazar, expandir
  o auditar los nodos de la wiki de la tesis en modo WIKILLM / LLMWiki.
---

# Navegador y Mantenedor del Grafo LLMWiki (`llmwiki_tesis_datura/`)

Esta *skill* mantiene la coherencia estructural y semántica de los **8 nodos modulares** en `llmwiki_tesis_datura/`:

## 1. Topología Oficial de 8 Nodos

1. `00_INDEX.md` (`00_INDEX`): Grafo maestro y arquitectura de 4 módulos operativos + 3 nodos transversales.
2. `01_MARCO_Y_OBJETIVOS.md` (`01_MARCO_Y_OBJETIVOS`): Problema (interacción holobionte, homología vs. expresión activa, reproducibilidad), pregunta, OE1–OE3 y bifurcación del OE2 (Rama A lecturas vs. Rama B contigs).
3. `02_MUESTRAS_Y_WETLAB.md` (`02_MUESTRAS_Y_WETLAB`): 12 muestras, nichos altitudinales (Subclado A Alta montaña vs. Subclado B Tierras bajas, Clado *Datura*, *Outgroups*), genomas de referencia y protocolo Trizol/NovaSeq 6000 S4.
4. `03_MODULO1_SEPARACION_Y_ENSAMBLAJE.md` (`03_MODULO1_SEPARACION_Y_ENSAMBLAJE`): `FastQC`, `MultiQC`, `Trimmomatic`, `RNA-QC-Chain`, `FastqPuri` (filtro de Bloom), `rnaSPAdes` ($K=49, 73$) y validación estructural (`CD-HIT-EST`, `TransRate`, `Bowtie2`, `SAMtools`).
5. `04_MODULO2_PREDICCION_ORF_Y_ANOTACION.md` (`04_MODULO2_PREDICCION_ORF_Y_ANOTACION`): `TransDecoder` (eucariota) vs. `Prodigal` (procariota), catálogo de 35 genes TA (`PSI-BLAST` + `MAFFT`) y anotación con `DIAMOND`.
6. `05_MODULO3_CUANTIFICACION_Y_DESEQ2.md` (`05_MODULO3_CUANTIFICACION_Y_DESEQ2`): `Salmon`, offset matemático en `tximport`, embudo de 5 pasos ($123.099 \rightarrow 54.363$ transcritos vía VST y MAD) y expresión diferencial en `DESeq2`.
7. `06_MODULO4_METABOLOMICA_Y_COEXPRESION.md` (`06_MODULO4_METABOLOMICA_Y_COEXPRESION`): Curaduría MS en 4 etapas ($82 \rightarrow 33 \rightarrow 15\text{–}25$), comparativa Pablo vs. Ronald vs. Tesis actual, `WGCNA`, redes duales planta-endófito y métodos para *"Small $n$, Large $p$"*.
8. `07_AUDITORIA_Y_REFERENCIAS.md` (`07_AUDITORIA_Y_REFERENCIAS`): Problemas en matrices (`B27`/`B20`, 55 columnas vacías, $n=10$), inconsistencias en `NewMethods.docx` y estado de referencias bibliográficas.

## 2. Convenciones de Formato WIKILLM

- **Frontmatter YAML**: Todo archivo debe iniciar con `id_nodo`, `tipo` y `nodos_relacionados`.
- **Wikilinks**: Usar siempre enlaces estilo `[[ID_NODO]]` (por ejemplo, `[[05_MODULO3_CUANTIFICACION_Y_DESEQ2]]`) para conectar conceptos entre módulos.
- **Fórmulas y Umbrales**: Conservar la notación matemática en LaTeX y las tablas comparativas parametrizadas sin resumir ni omitir criterios numéricos.
