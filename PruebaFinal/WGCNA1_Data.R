# ==============================================================================
# FASE 1 WGCNA: RED COMPLETA (54.000 TRANSCRITOS) EN 16 GB RAM
# ==============================================================================
library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# 1. Cargar la matriz curada completa (Sin recortes arbitrarios)
message("Cargando matriz VST completa...")
datExpr <- read.csv("matriz_vst_filtrada_para_WGCNA.csv", row.names = 1, check.names = FALSE)

cat(sprintf("Matriz cargada: %d muestras x %d transcritos\n", nrow(datExpr), ncol(datExpr)))

# 2. Control de calidad de genes y muestras
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0) message(paste("Removiendo genes con baja calidad:", sum(!gsg$goodGenes)))
  if (sum(!gsg$goodSamples) > 0) message(paste("Removiendo muestras con baja calidad:", sum(!gsg$goodSamples)))
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# 3. Limpieza de memoria previa al cálculo pesado
gc()

# 4. Construcción de módulos en bloques optimizados para 16 GB
message("Iniciando blockwiseModules para 54.000 genes...")
net <- blockwiseModules(
  datExpr,
  power = 6,                      # Potencia beta seleccionada
  networkType = "signed",
  TOMType = "signed",
  corType = "pearson",
  minModuleSize = 30,             # Tamaño mínimo adecuado para matrices de >50k
  mergeCutHeight = 0.25,          # Fusión de módulos con correlación > 0.75
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  
  # --- PARÁMETROS CRÍTICOS PARA 16 GB RAM ---
  maxBlockSize = 10000,           # Divide los ~54k genes en ~6 bloques manejables
  saveTOMs = TRUE,                # Guarda matrices TOM en disco para liberar RAM
  saveTOMFileBase = "TOM_bloques_16GB",
  # ------------------------------------------
  
  verbose = 3
)

# 5. Extraer módulos y perfiles representativos (Eigengenes)
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs

# 6. Guardar el entorno para Fase 2
save(datExpr, net, moduleColors, MEs, file = "WGCNA_Phase1_RedCompleta_16GB.RData")
message("¡Fase 1 completada exitosamente! Guardado en 'WGCNA_Phase1_RedCompleta_16GB.RData'")
