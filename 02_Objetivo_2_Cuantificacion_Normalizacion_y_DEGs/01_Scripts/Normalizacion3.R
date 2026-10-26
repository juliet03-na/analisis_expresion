# ==============================================================================
# PIPELINE DE PREPROCESAMIENTO: tximport + Filtrado Exacto + VST (+ MAD)
# SET: PLANTAS (Datureae / Brugmansia) - OPTIMIZADO PARA WGCNA
# ==============================================================================
# NOTA DE ENSAMBLAJE (Basado en evaluación CD-HIT-EST / TransRate):
# - Number of Contigs (Crudos): 135,008
# - Number of Contigs post-cluster: 123,099 (Reducción del 8.8%)
# - Tasa de alineamiento global: 98.3%
# ==============================================================================

library(tximport)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(writexl)

# -------------------------------------------------------------
# 1. Definir directorios e Importar con tximport
# -------------------------------------------------------------
if (dir.exists("plantas/plantas")) {
  base_dir <- normalizePath("plantas/plantas", winslash = "/")
} else {
  base_dir <- "/home/camila/Tesis/Conteos_Salmon/plantas"
  setwd(base_dir)
}

sample_folders <- c("B9_quant", "B10_quant", "B11_quant", "B12_quant",
                    "B14_quant", "B15_quant", "B17_quant", "B19_quant",
                    "B21_quant", "B27_quant")
sample_names <- gsub("_quant", "", sample_folders)
quant_files <- file.path(base_dir, sample_folders, "quant.sf")
names(quant_files) <- sample_names

stopifnot(all(file.exists(quant_files)))

colData <- data.frame(sample = sample_names, row.names = sample_names)

message("Importando quant.sf con tximport...")
txi <- tximport(quant_files, type = "salmon", txOut = TRUE, dropInfReps = TRUE)
dds <- DESeqDataSetFromTximport(txi, colData = colData, design = ~ 1)

# -------------------------------------------------------------
# 2. EMBUDO DE FILTRADO (Basado en el diagrama de flujo)
# -------------------------------------------------------------
message("Iniciando filtrado según diagrama de flujo...")
raw_counts <- counts(dds, normalized = FALSE)

# NOTA: El Paso 0 ("Transcriptoma ensamblado") capturará automáticamente
# los 123,099 transcritos post-cluster asumiendo que tu índice de Salmon 
# se construyó con ese archivo FASTA filtrado.
pasos <- c("0. Transcriptoma ensamblado (Post-cluster CD-HIT-EST)",
           "1. Eliminación suma de conteos = 0",
           "2. Eliminación varianza = 0",
           "3. Eliminación suma de conteos < 25",
           "4. Eliminación conteos > 50.000",
           "5. Eliminación conteos > 30.000",
           "6. Remover duplicados")

retenidos <- numeric(7)
eliminados <- numeric(7)

# Paso 0: Transcritos iniciales (Debería reflejar ~123,099 según métricas)
retenidos[1] <- nrow(raw_counts)
eliminados[1] <- 0

# Paso 1: Suma = 0
raw_counts <- raw_counts[rowSums(raw_counts) > 0, ]
retenidos[2] <- nrow(raw_counts); eliminados[2] <- retenidos[1] - retenidos[2]

# Paso 2: Varianza = 0
varianzas <- apply(raw_counts, 1, var)
raw_counts <- raw_counts[varianzas > 0, ]
retenidos[3] <- nrow(raw_counts); eliminados[3] <- retenidos[2] - retenidos[3]

# Paso 3: Suma < 25
raw_counts <- raw_counts[rowSums(raw_counts) >= 25, ]
retenidos[4] <- nrow(raw_counts); eliminados[4] <- retenidos[3] - retenidos[4]

# Paso 4: Conteos > 50.000
max_counts <- apply(raw_counts, 1, max)
raw_counts <- raw_counts[max_counts <= 50000, ]
retenidos[5] <- nrow(raw_counts); eliminados[5] <- retenidos[4] - retenidos[5]

# Paso 5: Conteos > 30.000
max_counts2 <- apply(raw_counts, 1, max)
raw_counts <- raw_counts[max_counts2 <= 30000, ]
retenidos[6] <- nrow(raw_counts); eliminados[6] <- retenidos[5] - retenidos[6]

# Paso 6: Remover duplicados
raw_counts <- raw_counts[!duplicated(rownames(raw_counts)), ]
retenidos[7] <- nrow(raw_counts); eliminados[7] <- retenidos[6] - retenidos[7]

# Exportar tabla del embudo
tabla_embudo <- data.frame(Paso_Aplicado = pasos, Genes_Retenidos = retenidos, Genes_Eliminados = eliminados)
write.csv(tabla_embudo, "Tabla_Embudo_Filtrado_Plantas_Corregido.csv", row.names = FALSE)
message("Generado: Tabla_Embudo_Filtrado_Plantas_Corregido.csv")

# -------------------------------------------------------------
# 3. TRANSFORMACIÓN VST Y FILTRO MAD 
# -------------------------------------------------------------
message("Calculando VST y aplicando filtro MAD...")
dds_clean <- dds[rownames(raw_counts), ]
vsd <- vst(dds_clean, blind = TRUE)
vst_matrix <- assay(vsd)

# Filtro MAD (Retener el 75% superior con mayor variación; como >25% tienen MAD=0,
# mads > quantile(mads, 0.25) retiene exactamente los 36.647 transcritos con MAD > 0)
mads <- apply(vst_matrix, 1, mad)
umbral_mad <- quantile(mads, 0.25)
vst_matrix_mad <- vst_matrix[mads > umbral_mad, ]

# Registrar también el Paso 7 (Post-MAD = 36.647) en la tabla del embudo
tabla_embudo <- rbind(
  tabla_embudo,
  data.frame(
    Paso_Aplicado = "7. Filtro de Varianza MAD sobre VST (mads > quantile 0.25)",
    Genes_Retenidos = nrow(vst_matrix_mad),
    Genes_Eliminados = nrow(vst_matrix) - nrow(vst_matrix_mad)
  )
)
write.csv(tabla_embudo, "Tabla_Embudo_Filtrado_Plantas_Corregido.csv", row.names = FALSE)

# -------------------------------------------------------------
# 4. EXPORTACIÓN DE MATRICES
# -------------------------------------------------------------
write.csv(vst_matrix_mad, "matriz_vst_filtrada_genes_x_muestras.csv")
message("Generado: matriz_vst_filtrada_genes_x_muestras.csv")

write.csv(t(vst_matrix_mad), "matriz_vst_filtrada_para_WGCNA.csv")
message("Generado: matriz_vst_filtrada_para_WGCNA.csv")

df_excel <- data.frame(Transcrito = rownames(vst_matrix_mad), vst_matrix_mad, check.names = FALSE)
write_xlsx(df_excel, "Matriz_VST_Visualizacion_Limpia.xlsx")
message("Generado: Matriz_VST_Visualizacion_Limpia.xlsx")

# -------------------------------------------------------------
# 5. GRÁFICOS DE CONTROL DE CALIDAD (QC)
# -------------------------------------------------------------
message("Generando PDFs de Control de Calidad...")

pca_data <- plotPCA(vsd, intgroup = "sample", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))
pdf("QC_PCA_muestras_filtradas.pdf", width = 7, height = 6)
print(
  ggplot(pca_data, aes(PC1, PC2, label = name)) +
    geom_point(size = 3.5, color = "#2c7bb6") +
    geom_text(vjust = -0.8, size = 3.5) +
    xlab(paste0("PC1: ", percentVar[1], "% de la varianza")) +
    ylab(paste0("PC2: ", percentVar[2], "% de la varianza")) +
    ggtitle("PCA de Muestras (VST)") +
    theme_minimal()
)
dev.off()

pdf("QC_Dendrograma_muestras_filtradas.pdf", width = 7, height = 5)
sample_dists <- dist(t(vst_matrix_mad))
plot(hclust(sample_dists, method = "average"),
     main = "Dendrograma de Muestras (VST-MAD)", xlab = "Muestras", sub = "")
dev.off()

pdf("QC_Heatmap_correlacion_filtrado.pdf", width = 8, height = 7)
pheatmap(cor(vst_matrix_mad, method = "pearson"),
         main = "Correlación de Pearson entre Muestras (VST-MAD)")
dev.off()

message("¡Proceso completado exitosamente! Todos los archivos han sido generados.")
