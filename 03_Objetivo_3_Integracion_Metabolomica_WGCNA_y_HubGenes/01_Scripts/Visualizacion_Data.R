# ==============================================================================
# VISUALIZACIÓN DE DENDROGRAMAS MULTIBLOQUE Y TAMAÑO DE MÓDULOS (16 GB)
# ==============================================================================
library(WGCNA)
load("WGCNA_Phase1_RedCompleta_16GB.RData")

pdf("WGCNA_03_Dendrogramas_TodosLosBloques.pdf", width = 12, height = 9)

# 1. Graficar el árbol de cada bloque generado
num_bloques <- length(net$dendrograms)

for (b in seq_len(num_bloques)) {
  plotDendroAndColors(
    dendro = net$dendrograms[[b]],
    colors = moduleColors[net$blockGenes[[b]]],
    groupLabels = paste("Módulos Bloque", b),
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05,
    main = sprintf(
      "Dendrograma de Co-expresión - Bloque %d de %d (%d transcritos)",
      b, num_bloques, length(net$blockGenes[[b]])
    )
  )
}

# 2. Gráfico de barras con el tamaño de los principales módulos
par(mar = c(7, 5, 4, 2))
conteo_modulos <- sort(table(moduleColors), decreasing = TRUE)
top_n <- min(25, length(conteo_modulos))
barplot(
  conteo_modulos[1:top_n],
  col = names(conteo_modulos)[1:top_n],
  las = 2,
  ylab = "Número de Transcritos",
  main = sprintf("Top %d Módulos por Número de Transcritos", top_n),
  cex.names = 0.85
)

dev.off()
message("Dendrogramas exportados exitosamente a 'WGCNA_03_Dendrogramas_TodosLosBloques.pdf'")
