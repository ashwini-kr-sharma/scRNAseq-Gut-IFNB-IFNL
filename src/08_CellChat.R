# Some notes on installing CellChat

# devtools::install_github("sqjin/CellChat")

# Doing this will show the list of all dependencies which also needs to be installed
# Install all of those manually first

# Some of them depend on systemfonts, whose installations throws error. Downgrade the two below
# devtools::install_version("systemfonts", version = "0.3.2", repos = "http://cran.us.r-project.org")
# devtools::install_version("svglite", version = "1.2.2", repos = "http://cran.us.r-project.org")

# devtools::install_github("sqjin/CellChat")
# No it should install properly

library(CellChat)
library(patchwork)
library(Seurat)
path = "/icgc/dkfzlsdf/analysis/B080/sharma/boulantLab/sc_interferon/"

# Read data
seu = readRDS(paste0(path, "analysis/HJHG3BGXF_Seurat_Combined.RDS"))

# Ligand-receptor interaction database
CellChatDB = CellChatDB.human
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)

# Create CellChat object
cc = lapply(seu, function(x){ 
 x = createCellChat(object = x$SeuratObject, group.by = "seurat_clusters") 
 x@DB = CellChatDB
 x = subsetData(x)
 x = identifyOverExpressedGenes(x)
 x = identifyOverExpressedInteractions(x)
 #x = projectData(cellchat, PPI.human)
 x = computeCommunProb(x, raw.use = TRUE) # use raw.data=FALSE, if uncommenting the above line
 x = filterCommunication(x, min.cells = 10)
 x = computeCommunProbPathway(x)
 x = aggregateNet(x)
 return(x)
})

# Plotting
tmp = cc$`Ileum|Beta_Lambda`

# All in one
groupSize <- as.numeric(table(tmp@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(tmp@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(tmp@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# Per cluster
mat <- tmp@net$weight
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

# Specific pathways
pathways.show <- c("GRN") 
netVisual_aggregate(tmp, signaling = pathways.show,  layout = "circle")
netVisual_aggregate(tmp, signaling = pathways.show,  layout = "chord")
netVisual_heatmap(tmp, signaling = pathways.show, color.heatmap = "Reds")

netAnalysis_contribution(tmp, signaling = pathways.show)

netVisual_bubble(tmp, sources.use = c(1:8), targets.use = c(1:8), remove.isolate = FALSE)

