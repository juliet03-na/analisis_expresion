---
id_nodo: 01_MARCO_Y_OBJETIVOS
tipo: nodo_transversal
nodos_relacionados:
  - "[[00_INDEX]]"
  - "[[03_MODULO1_SEPARACION_Y_ENSAMBLAJE]]"
  - "[[05_MODULO3_CUANTIFICACION_Y_DESEQ2]]"
  - "[[06_MODULO4_METABOLOMICA_Y_COEXPRESION]]"
---

# `01_MARCO_Y_OBJETIVOS` — Problema, Objetivos, Ramas Analíticas y Subclados

- **ID_Nodo**: `01_MARCO_Y_OBJETIVOS`
- **Nodos_Relacionados**: [[00_INDEX]], [[03_MODULO1_SEPARACION_Y_ENSAMBLAJE]], [[05_MODULO3_CUANTIFICACION_Y_DESEQ2]], [[06_MODULO4_METABOLOMICA_Y_COEXPRESION]]

---

## 1. Síntesis del Problema y Justificación

- **Interacción holobionte**: El intercambio de señales y metabolitos exige analizar de forma integrada los componentes de la planta hospedera y su microbiota asociada (Fitzpatrick et al., 2020).
- **Homología vs. Expresión activa**: Aunque se han reportado endófitos en Solanaceae con genes homólogos capaces de sintetizar hiosciamina y escopolamina (Naik et al., 2018), la presencia del gen no prueba su actividad funcional, requiriendo evidencia transcripcional simultánea en planta y microorganismos.
- **Desafío bioinformático y reproducibilidad**: La sobrerrepresentación de lecturas vegetales (Taj et al., 2023) y la diversidad de herramientas analíticas (Shakya et al., 2019) obligan a implementar un flujo computacional reproducible con trazabilidad de versiones y parámetros (Di Leo et al., 2025; Gruening et al., 2023).

---

## 2. Pregunta y Objetivos

### Pregunta de Investigación
> ¿Cómo desarrollar un flujo computacional reproducible que integre la discriminación del origen de los transcritos entre la planta hospedera y sus microorganismos asociados, para caracterizar la expresión de genes precursores y biosintéticos de la ruta de alcaloides tropánicos y su correlación con los perfiles metabolómicos en especies del clado *Datura*?

### Objetivo General
Desarrollar un flujo computacional reproducible para caracterizar la expresión de genes del metabolismo de precursores y de la biosíntesis de alcaloides tropánicos en plantas y microorganismos asociados del clado *Datura*.

### Objetivos Específicos (OE)
- **OE1**: Implementar una estrategia computacional para discriminar los transcritos de la planta hospedera y de sus microorganismos asociados en datos metatranscriptómicos de especies del clado *Datura* ([[03_MODULO1_SEPARACION_Y_ENSAMBLAJE]]).
- **OE2**: Determinar los patrones de expresión de genes precursores y biosintéticos de la ruta de alcaloides tropánicos entre la planta hospedera y sus microorganismos asociados en especies del clado *Datura* ([[04_MODULO2_PREDICCION_ORF_Y_ANOTACION]], [[05_MODULO3_CUANTIFICACION_Y_DESEQ2]]).
- **OE3**: Correlacionar los perfiles de expresión de genes precursores y biosintéticos de la ruta con la asociación de perfiles metabolómicos previamente reportados en especies del clado *Datura* ([[06_MODULO4_METABOLOMICA_Y_COEXPRESION]]).

---

## 3. Bifurcación Analítica del OE2 (Lecturas vs. Ensamblaje)

| Estrategia | Dependencia de Ensamblaje | Flujo de Procesamiento | Propósito Biológico |
| :--- | :--- | :--- | :--- |
| **Rama A: A nivel de lecturas (Dirigida)** | No dependiente del ensamblaje. | Lecturas filtradas $\rightarrow$ Alineamiento directo (`DIAMOND blastx`) contra genes de referencia y precursores $\rightarrow$ Conteos $\rightarrow$ Perfiles de expresión de biosíntesis. | Reproducir el flujo de referencia sobre genes conocidos (clados de alta vs. baja montaña reportados por Hannah; validación de Naik et al.). |
| **Rama B: A nivel de contigs (*De novo*)** | Dependiente del co-ensamblaje. | Co-ensamblaje *de novo* (`rnaSPAdes`) vegetal y microbiano $\rightarrow$ Contigs depurados (`CD-HIT-EST`) $\rightarrow$ Cuantificación (`Salmon`) y predicción proteica (`TransDecoder`/`Prodigal`). | Determinar nuevas secuencias de referencia de la ruta y alimentar redes no supervisadas (`WGCNA`) para proponer nuevos candidatos. |
