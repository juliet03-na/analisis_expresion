---
id_nodo: 02_MUESTRAS_Y_WETLAB
tipo: nodo_transversal
nodos_relacionados:
  - "[[00_INDEX]]"
  - "[[03_MODULO1_SEPARACION_Y_ENSAMBLAJE]]"
  - "[[05_MODULO3_CUANTIFICACION_Y_DESEQ2]]"
  - "[[07_AUDITORIA_Y_REFERENCIAS]]"
---

# `02_MUESTRAS_Y_WETLAB` — Muestras, Nichos Altitudinales, Genomas y ARN

- **ID_Nodo**: `02_MUESTRAS_Y_WETLAB`
- **Nodos_Relacionados**: [[00_INDEX]], [[03_MODULO1_SEPARACION_Y_ENSAMBLAJE]], [[05_MODULO3_CUANTIFICACION_Y_DESEQ2]]

---

## 1. Inventario de Muestras y Clasificación por Nicho Ecológico (Subclados)

- **Origen**: Colección viva en invernadero del grupo Max Planck Tandem GEME (Universidad Nacional de Colombia), recolectada bajo permiso reportado en Hernández-Duarte et al. (2025) y preservada a $-80\text{ °C}$.

| No. | Sample ID | Especie | Tejido | Clasificación Taxonómica / Nicho Altitudinal (PCA y VST) | Presencia en Matrices (`Expresión` vs. `Metabolómica`) |
| :---: | :---: | :--- | :---: | :--- | :--- |
| 1 | `B10` | *Brugmansia sanguinea* | Hoja | **Subclado A**: Alta montaña (Tierras altas) | Ambas (cruce homologado $n=8$) |
| 2 | `B11` | *Brugmansia sanguinea x vulcanicola* | Hoja | **Subclado A**: Alta montaña (Tierras altas) | Ambas (cruce homologado $n=8$) |
| 3 | `B15` | *Brugmansia vulcanicola* | Hoja | **Subclado A**: Alta montaña (Tierras altas) | Ambas (cruce homologado $n=8$) |
| 4 | `P1` | *Brugmansia arborea* | Hoja | **Subclado A**: Alta montaña (Tierras altas) | Ausente en matrices finales exportadas (`n=10`) |
| 5 | `B9` | *Brugmansia insignis* | Hoja | **Subclado B**: Tierras bajas | Ambas (cruce homologado $n=8$) |
| 6 | `B12` | *Brugmansia suaveolens* | Hoja | **Subclado B**: Tierras bajas | Ambas (cruce homologado $n=8$) |
| 7 | `B14` | *Brugmansia versicolor* | Hoja | **Subclado B**: Tierras bajas | Ambas (cruce homologado $n=8$) |
| 8 | `B21` | *Brugmansia aurea* | Hoja | **Subclado B**: Tierras bajas | Ambas (cruce homologado $n=8$) |
| 9 | `B17` | *Datura stramonium* | Hoja | **Clado Datura** (agrupa con `B27` en VST) | Ambas (cruce homologado $n=8$) |
| 10 | `B27` | *Datura wrightii* | Hoja | **Clado Datura** | **Solo en Expresión** (`matriz_vst_filtrada_para_WGCNA.csv`) |
| 11 | `B19` | *Atropa belladonna* | Hoja | **Outgroup clado Datura** (control externo en PCA/VST) | Ambas (se excluye en `WGCNA` para evitar sesgo taxonómico) |
| 12 | `B20` | *Ipomoea purpurea* | Hoja | **Outgroup familia Solanaceae** | **Solo en Metabolómica** (`Metaboloma_WGNA_TESIS.xlsx`) |

> [!WARNING]
> **Trazabilidad exacta de tamaños muestrales ($12 \rightarrow 10 \rightarrow 9 \rightarrow 8$)**:
> - **12 muestras en inventario biológico**: `B9, B10, B11, B12, B14, B15, P1, B21, B17, B27, B19, B20`.
> - **10 muestras en `matriz_vst_filtrada_para_WGCNA.csv` (Expresión)**: `B9, B10, B11, B12, B14, B15, B17, B19, B21, B27` (sin `P1` ni `B20`).
> - **10 muestras en `Metaboloma_WGNA_TESIS.xlsx` (Metabolómica)**: `B21, B9, B12, B14, B10, B11, B15, B19, B17, B20` (en distinto orden de filas, con `B20` y sin `B27`).
> - **9 muestras en Red `WGCNA` depurada (Expresión sin *outgroup* `B19`)**: `B9, B10, B11, B12, B14, B15, B17, B21, B27`.
> - **8 muestras en Cruce Fase 2 Expresión $\leftrightarrow$ Metabolómica (sin *outgroups* `B19`/`B20` ni `B27`)**: `B9, B10, B11, B12, B14, B15, B17, B21` (`fenotipos_TA_homologados_n8_SinOutgroups.csv`). Ver [[07_AUDITORIA_Y_REFERENCIAS]].

---

## 2. Genomas de Referencia para Descontaminación y Filtrado

- **Panel Solanaceae**:
  - *Brugmansia arborea* (`PRJNA765960`)
  - *Datura stramonium* (`GCA_018107945.1` y un segundo ensamblaje de *Datura*)
  - *Atropa belladonna* (`PRJNA766188`)
  - *Solanum lycopersicum* (`GCF_036512215.1`)
- **Panel Viral**: Base de datos Virus RefSeq de NCBI.

---

## 3. Protocolo Húmedo (Extracción de ARN y Secuenciación)

### Extracción con Trizol/Cloroformo
Basada en Bedewitz et al. (2018), Schlesinger et al. (2021) y Trapnell et al. (2012):
1. Maceración en $\text{N}_2$ líquido + $1\text{ mL}$ Trizol por $50\text{–}100\text{ mg}$ de hoja $\rightarrow$ vórtex e incubación $5\text{ min}$ a temperatura ambiente (TA).
2. Cloroformo $1:5$ respecto a Trizol $\rightarrow$ inversión manual $15\text{ s}$ $\rightarrow$ incubación $2\text{–}3\text{ min}$ a TA $\rightarrow$ centrifugación a $12.000\text{ rpm}$ por $15\text{ min}$ a $4\text{ °C}$.
3. Precipitación de fase acuosa con $0,5\text{ mL}$ de isopropanol ($10\text{ min}$ a TA) $\rightarrow$ centrifugación a $12.500\text{ rpm}$ por $15\text{ min}$ a $4\text{ °C}$ $\rightarrow$ doble lavado del pellet con etanol frío al $75\%$.

### Control de Calidad y Secuenciación
- **Control de calidad**: Pureza $260/280$ y $260/230$ en NanoDrop, cuantificación en Qubit e integridad en gel de agarosa.
- **Secuenciación**: Preparación de bibliotecas (depleción de ARNr, síntesis de ADNc, enriquecimiento PCR) y secuenciación ejecutadas por NovoGene en plataforma **Illumina NovaSeq 6000 S4** (*paired-end* $150\text{ pb}$, $6\text{–}8\text{ GB}$ por muestra en formato `FASTQ`).
