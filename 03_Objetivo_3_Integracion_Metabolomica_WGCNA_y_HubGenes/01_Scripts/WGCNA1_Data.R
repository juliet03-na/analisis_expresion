# ==============================================================================
# FASE 1 WGCNA: RED COMPLETA (36.647 / 54.363 TRANSCRITOS) EN 16 GB RAM
# Tesis Hologenómica Clado Datura — Módulo 4 (Co-expresión)
# ==============================================================================
library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

t_inicio <- Sys.time()

# ------------------------------------------------------------------------------
# 0. Parámetros de configuración general
# ------------------------------------------------------------------------------
ARCHIVO_ENTRADA      <- "matriz_vst_filtrada_para_WGCNA.csv"
REMOVER_OUTGROUP_B19 <- TRUE    # TRUE recomendado: evita que B19 (Atropa belladonna) secuestre >3.200 genes
POTENCIA_FIJA        <- 12      # Potencia beta recomendada para red signed con n=9-10 (cambiar a 6 para comparar)
EJECUTAR_PICK_POWER  <- TRUE    # Generar gráfico diagnóstico de topología libre de escala

# ------------------------------------------------------------------------------
# 1. Cargar la matriz curada completa (Sin recortes arbitrarios)
# ------------------------------------------------------------------------------
message("1. Cargando matriz VST completa desde: ", ARCHIVO_ENTRADA)
datExpr <- read.csv(ARCHIVO_ENTRADA, row.names = 1, check.names = FALSE)

cat(sprintf("Matriz cargada: %d muestras x %d transcritos\n", nrow(datExpr), ncol(datExpr)))
cat("Muestras detectadas:", paste(rownames(datExpr), collapse = ", "), "\n")

if (REMOVER_OUTGROUP_B19 && "B19" %in% rownames(datExpr)) {
  message("Excluyendo outgroup B19 (Atropa belladonna) del análisis de co-expresión...")
  datExpr <- datExpr[rownames(datExpr) != "B19", , drop = FALSE]
  cat(sprintf("Dimensiones tras excluir B19: %d muestras x %d transcritos\n", nrow(datExpr), ncol(datExpr)))
}

# ------------------------------------------------------------------------------
# 2. Control de calidad de genes y muestras (goodSamplesGenes)
# ------------------------------------------------------------------------------
message("2. Verificando calidad de muestras y transcritos (goodSamplesGenes)...")
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0) {
    message(paste("Removiendo genes con varianza cero o demasiados NA:", sum(!gsg$goodGenes)))
  }
  if (sum(!gsg$goodSamples) > 0) {
    message(paste("Removiendo muestras con baja calidad:", sum(!gsg$goodSamples)))
  }
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
}

# ------------------------------------------------------------------------------
# 2.1 Dendrograma de agrupamiento de muestras (Detección de outliers / B19)
# ------------------------------------------------------------------------------
sampleTree <- hclust(dist(datExpr), method = "average")
pdf("WGCNA_00_Clustering_Muestras.pdf", width = 10, height = 6)
par(cex = 0.9, mar = c(0, 4, 2, 0))
plot(
  sampleTree,
  main = "Agrupamiento Jerárquico de Muestras (Distancia Euclidiana sobre VST)",
  sub = "",
  xlab = ""
)
dev.off()
message("Dendrograma de muestras guardado en 'WGCNA_00_Clustering_Muestras.pdf'")

# ------------------------------------------------------------------------------
# 2.2 Diagnóstico de Topología Libre de Escala (pickSoftThreshold)
# ------------------------------------------------------------------------------
if (EJECUTAR_PICK_POWER) {
  message("Evaluando potencia suave (pickSoftThreshold) para red 'signed'...")
  powers <- c(1:10, seq(from = 12, to = 20, by = 2))
  sft <- pickSoftThreshold(
    datExpr,
    powerVector = powers,
    networkType = "signed",
    corFnc = "cor",
    verbose = 2
  )
  write.csv(sft$fitIndices, "WGCNA_01_Tabla_SoftThreshold.csv", row.names = FALSE)

  pdf("WGCNA_01_Topologia_Libre_Escala.pdf", width = 10, height = 5)
  par(mfrow = c(1, 2))
  cex1 <- 0.9
  # Ajuste R^2 a topología libre de escala
  plot(
    sft$fitIndices[, 1],
    -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    xlab = "Potencia Suave (Soft Threshold beta)",
    ylab = "Ajuste Topología Libre de Escala, R^2 con signo",
    type = "n",
    main = "Independencia de Escala (Red Signed)"
  )
  text(
    sft$fitIndices[, 1],
    -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    labels = powers,
    cex = cex1,
    col = "red"
  )
  abline(h = 0.80, col = "blue", lty = 2)
  abline(h = 0.85, col = "darkgreen", lty = 2)

  # Conectividad media
  plot(
    sft$fitIndices[, 1],
    sft$fitIndices[, 5],
    xlab = "Potencia Suave (Soft Threshold beta)",
    ylab = "Conectividad Media",
    type = "n",
    main = "Conectividad Media"
  )
  text(
    sft$fitIndices[, 1],
    sft$fitIndices[, 5],
    labels = powers,
    cex = cex1,
    col = "red"
  )
  dev.off()
  message("Diagnóstico de potencia guardado en 'WGCNA_01_Topologia_Libre_Escala.pdf'")
}

# ------------------------------------------------------------------------------
# 3. Limpieza de memoria y blindaje del namespace de correlación
# ------------------------------------------------------------------------------
gc()
# Evita el error 'unused argument (weights.x = NULL)' causado por stats::cor
cor_original <- stats::cor
cor <- WGCNA::cor

# ------------------------------------------------------------------------------
# 4. Construcción de módulos en bloques optimizados para 16 GB RAM
# ------------------------------------------------------------------------------
message(sprintf("4. Iniciando blockwiseModules para %d transcritos con power = %d...", ncol(datExpr), POTENCIA_FIJA))
net <- blockwiseModules(
  datExpr,
  power = POTENCIA_FIJA,          # Potencia beta seleccionada (6 por defecto; evaluar 12-16)
  networkType = "signed",         # Agrupa genes que co-varían en la misma dirección
  TOMType = "signed",
  corType = "pearson",
  deepSplit = 2.5,                # Sensibilidad de corte de ramas (documentado en Módulo 4)
  minModuleSize = 30,             # Tamaño mínimo adecuado para matrices de >35k-54k
  reassignThreshold = 0,          # Umbral de reasignación entre módulos
  mergeCutHeight = 0.25,          # Fusión de módulos con correlación de Eigengenes >= 0.75
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,

  # --- PARÁMETROS CRÍTICOS PARA 16 GB RAM ---
  maxBlockSize = 10000,           # Divide los 36.647 transcritos en 4 bloques manejables
  saveTOMs = TRUE,                # Guarda matrices TOM en disco para liberar RAM
  saveTOMFileBase = "TOM_bloques_16GB",
  # ------------------------------------------

  verbose = 3
)

# Restaurar stats::cor
cor <- cor_original

# ------------------------------------------------------------------------------
# 5. Extraer módulos, perfiles representativos (Eigengenes) y diagnóstico
# ------------------------------------------------------------------------------
moduleColors <- labels2colors(net$colors)
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)

tabla_modulos <- as.data.frame(table(Color_Modulo = moduleColors))
tabla_modulos <- tabla_modulos[order(-tabla_modulos$Freq), ]
tabla_modulos$Porcentaje <- round(100 * tabla_modulos$Freq / ncol(datExpr), 2)

pct_grey <- ifelse(
  "grey" %in% tabla_modulos$Color_Modulo,
  tabla_modulos$Porcentaje[tabla_modulos$Color_Modulo == "grey"],
  0
)

cat("\n==============================================================================\n")
cat("RESUMEN DE MÓDULOS DETECTADOS POR WGCNA\n")
cat("==============================================================================\n")
cat(sprintf("Total de transcritos analizados : %d\n", ncol(datExpr)))
cat(sprintf("Número de bloques procesados    : %d\n", length(net$dendrograms)))
cat(sprintf("Total de módulos (excl. grey)   : %d\n", sum(tabla_modulos$Color_Modulo != "grey")))
cat(sprintf("Porcentaje en módulo 'grey'     : %.2f%% (criterio óptimo: < 30%% - 40%%)\n", pct_grey))
print(head(tabla_modulos, 15), row.names = FALSE)

# Exportar tablas CSV de resultados para trazabilidad
write.csv(tabla_modulos, "WGCNA_02_Resumen_Modulos.csv", row.names = FALSE)
write.csv(MEs, "WGCNA_02_Module_Eigengenes.csv", row.names = TRUE)

asignacion_genes <- data.frame(
  Transcript_ID = colnames(datExpr),
  Module_Number = net$colors,
  Module_Color = moduleColors
)
write.csv(asignacion_genes, "WGCNA_02_Asignacion_Transcritos_Modulos.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 6. Guardar el entorno para Fase 2 y Visualización
# ------------------------------------------------------------------------------
save(datExpr, net, moduleColors, MEs, tabla_modulos, file = "WGCNA_Phase1_RedCompleta_16GB.RData")

t_fin <- Sys.time()
duracion <- round(difftime(t_fin, t_inicio, units = "mins"), 2)
message(sprintf("¡Fase 1 completada exitosamente en %s minutos! Guardado en 'WGCNA_Phase1_RedCompleta_16GB.RData'", duracion))
