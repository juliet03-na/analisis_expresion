# ==============================================================================
# VISUALIZACIÓN DE DENDROGRAMAS MULTIBLOQUE (16 GB)
# ==============================================================================
library(WGCNA)
load("WGCNA_Phase1_RedCompleta_16GB.RData")

pdf("WGCNA_03_Dendrogramas_TodosLosBloques.pdf", width = 12, height = 9)

# Graficar el árbol de cada bloque generado
num_bloques <- length(net$dendrograms)

for (b in 1:num_bloques) {
  plotDendroAndColors(
    dendro = net$dendrograms[[b]], 
    colors = moduleColors[net$blockGenes[[b]]],
    groupLabels = paste("Módulos Bloque", b),
    dendroLabels = FALSE, 
    hang = 0.03,
    addGuide = TRUE, 
    guideHang = 0.05,
    main = paste("Dendrograma de Coexpresión - Bloque", b, "de", num_bloques)
  )
}

dev.off()
message("Dendrogramas exportados a 'WGCNA_03_Dendrogramas_TodosLosBloques.pdf'")
