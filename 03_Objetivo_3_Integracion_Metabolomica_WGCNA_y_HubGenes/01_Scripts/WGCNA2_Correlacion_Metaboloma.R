# ==============================================================================
# FASE 2 WGCNA: CORRELACIÓN MÓDULOS DE CO-EXPRESIÓN vs. METABOLOMA TA CORREGIDO
# Tesis Hologenómica Clado Datura — Módulo 4
# ==============================================================================
library(WGCNA)
options(stringsAsFactors = FALSE)

# 1. Cargar entorno de Fase 1 (WGCNA_Phase1_RedCompleta_16GB.RData)
message("1. Cargando resultados de Fase 1...")
load("WGCNA_Phase1_RedCompleta_16GB.RData")

# 2. Cargar matriz metabolómica corregida desde Metaboloma_WGNA_TESIS.xlsx
# IMPORTANTE: No usar 'fenotipos_TA_limpios_para_WGCNA.csv' (contiene distorsión decimal 10x/100x)
ARCHIVO_METABOLOMA <- "fenotipos_TA_corregidos_para_WGCNA.csv"
message("2. Cargando matriz metabolómica corregida: ", ARCHIVO_METABOLOMA)
datTraitsRaw <- read.csv(ARCHIVO_METABOLOMA, row.names = 1, check.names = FALSE)

# 3. Homologación estricta de muestras (evita el error B27 vs B20 y desorden de filas)
common_samples <- intersect(rownames(MEs), rownames(datTraitsRaw))
message(sprintf("3. Muestras comunes homologadas (n = %d): %s",
                length(common_samples), paste(common_samples, collapse = ", ")))

MEs_aligned   <- MEs[common_samples, , drop = FALSE]
datExpr_align <- datExpr[common_samples, , drop = FALSE]
datTraits     <- datTraitsRaw[common_samples, , drop = FALSE]

stopifnot(identical(rownames(MEs_aligned), rownames(datTraits)))
stopifnot(identical(rownames(datExpr_align), rownames(datTraits)))

# Transformación log10 de áreas de pico MS para estabilizar asimetría
datTraits_log10 <- log10(datTraits)

# 4. Cálculo de Correlación de Pearson (Módulo - Alcaloide) y p-valores (t-Student)
nSamples <- nrow(datExpr_align)
moduleTraitCor <- WGCNA::cor(MEs_aligned, datTraits_log10, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

write.csv(round(moduleTraitCor, 4), "WGCNA_Fase2_Matriz_Correlacion_R.csv")
write.csv(round(moduleTraitPvalue, 5), "WGCNA_Fase2_Matriz_Pvalores_R.csv")

# 5. Generar Heatmap Módulo-Alcaloide (Top 30 compuestos con mayor asociación)
max_cor_per_trait <- apply(abs(moduleTraitCor), 2, max)
top_traits <- names(sort(max_cor_per_trait, decreasing = TRUE))[1:min(30, ncol(moduleTraitCor))]

subCor <- moduleTraitCor[, top_traits, drop = FALSE]
subP   <- moduleTraitPvalue[, top_traits, drop = FALSE]

textMatrix <- paste0(signif(subCor, 2), "\n(", signif(subP, 1), ")")
dim(textMatrix) <- dim(subCor)

pdf("WGCNA_04_Heatmap_Modulo_Alcaloide_Corregido.pdf", width = 16, height = 10)
par(mar = c(6, 10, 3, 2))
labeledHeatmap(
  Matrix = subCor,
  xLabels = colnames(subCor),
  yLabels = rownames(subCor),
  ySymbols = rownames(subCor),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.55,
  zlim = c(-1, 1),
  main = sprintf("Relación Módulo-Alcaloide Tropánico (n = %d muestras homologadas, Top 30 compuestos)", nSamples)
)
dev.off()

message("¡Fase 2 completada exitosamente! Heatmap exportado en 'WGCNA_04_Heatmap_Modulo_Alcaloide_Corregido.pdf'")
