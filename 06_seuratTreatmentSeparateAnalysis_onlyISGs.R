library(Seurat)
library(patchwork)
library(tidyverse)
library(rio)
library(RColorBrewer)

path = "/icgc/dkfzlsdf/analysis/B080/sharma/boulantLab/sc_interferon/"

#-------------------------------------------------------------------------------
# Signature gene lists
#-------------------------------------------------------------------------------

# ISG gene signature | PMC5963606
isg1 = c("IFNAR1", "IFNAR2", "IFNLR1", "IL10RB", 
         "CXCL10", "DDX60", "EPSTI1", "GBP1", "HERC5", "HERC6", "IFI27", "IFI44", 
         "IFI44L", "IFI6", "IFIT1", "IFIT2", "USP18", "CGAS", "IFIT3", "IFIT5", "ISG15",
         "ISG20","LAMP3", "LY6E", "MX1", "OAS1", "OAS2", "CH25H", "OAS3", "OASL", "PLSCR1", 
         "RSAD2", "RTP4", "SIGLEC1", "SOCS1", "SPATS2L", "ISGF3")

# ISGs from https://doi.org/10.1038/s41590-019-0323-3 (Supplementary Table 1)
isg2 = import("https://static-content.springer.com/esm/art%3A10.1038%2Fs41590-019-0323-3/MediaObjects/41590_2019_323_MOESM3_ESM.xlsx", sheet = 3)
isg2 = isg2[,1:4]
colnames(isg2) = isg2[1,]
isg2 = isg2[-1,]
isg2 = isg2$ISG

isg = unique(c(isg1, isg2))
rm(isg1, isg2)

#-------------------------------------------------------------------------------
# Identify data to read
#-------------------------------------------------------------------------------

anno = data.frame(rbind( 
  c("Int1", "Ileum|Mock"),
  c("Int2", "Ileum|Beta"),
  c("Int3", "Ileum|Lambda"),
  c("Int4", "Ileum|Beta_Lambda"),
  c("Int5", "Colon|Mock"),
  c("Int6", "Colon|Beta"),
  c("Int7", "Colon|Lambda"),
  c("Int8", "Colon|Beta_Lambda")
),
stringsAsFactors = F)
colnames(anno) = c("type", "treatment")

f = list.files(paste0(path, "analysis/postFilter_merge_counts/Seurat"), recursive = T, full.names = T)
id = sapply(strsplit(basename(f), "_", fixed=T), function(x)x[1])

# Read all the Seurat objects
seu = lapply(f, readRDS)

if(identical(id, anno$type)){
  names(seu) = anno$treatment
}
rm(anno, f, id)

#-------------------------------------------------------------------------------
# Preprocess Seurat function
#-------------------------------------------------------------------------------

preProcessSeurat = function(dat)
{
  # Subsetting to only ISGs
  dat = subset(dat, features = isg)
  
  # Seurat preprocessing
  dat <- NormalizeData(object = dat)
  dat <- FindVariableFeatures(object = dat, nfeatures = 5000)
  dat <- ScaleData(object = dat)
  dat <- RunPCA(object = dat, npcs = 30)
  dat <- FindNeighbors(object = dat)
  dat <- FindClusters(object = dat, resolution = 0.4)
  dat <- RunUMAP(object = dat, dim = 1:15)
  
  return(list(SeuratObject = dat, Markers = NA)) 
}

# Read in all objects
seu = lapply(seu, preProcessSeurat)

#-------------------------------------------------------------------------------
# Plot cell types
#-------------------------------------------------------------------------------

# Ileum

p1 = DimPlot(seu$`Ileum|Mock`$SeuratObject, group.by = "seurat_clusters", pt.size = 0.1, label = T, cols = c(brewer.pal(8, "Set2"),"black")) + NoLegend() + labs(title = "", subtitle = "Ileum|Mock") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()
p2 = DimPlot(seu$`Ileum|Beta`$SeuratObject, group.by = "seurat_clusters", pt.size = 0.1, label = T, cols = c(brewer.pal(8, "Set2"),"black")) + NoLegend() + labs(title = "", subtitle = "Ileum|Beta") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()
p3 = DimPlot(seu$`Ileum|Lambda`$SeuratObject, group.by = "seurat_clusters", pt.size = 0.1,label = T, cols = c(brewer.pal(8, "Set2"),"black")) + NoLegend() + labs(title = "", subtitle = "Ileum|Lambda") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()
p4 = DimPlot(seu$`Ileum|Beta_Lambda`$SeuratObject, group.by = "seurat_clusters", pt.size = 0.1, label = T, cols = c(brewer.pal(8, "Set2"),"black")) + NoLegend() + labs(title = "", subtitle = "Ileum|Beta+Lambda") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()

p5 = FeaturePlot(seu$`Ileum|Mock`$SeuratObject, features = "nFeature_RNA", pt.size = 0.1) + labs(title = "", subtitle = "Ileum|Mock - ISG score") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()
p6 = FeaturePlot(seu$`Ileum|Beta`$SeuratObject, features = "nFeature_RNA", pt.size = 0.1) + labs(title = "", subtitle = "Ileum|Beta - ISG score") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()
p7 = FeaturePlot(seu$`Ileum|Lambda`$SeuratObject, features = "nFeature_RNA", pt.size = 0.1) + labs(title = "", subtitle = "Ileum|Lambda - ISG score") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()
p8 = FeaturePlot(seu$`Ileum|Beta_Lambda`$SeuratObject, features = "nFeature_RNA", pt.size = 0.1) + labs(title = "", subtitle = "Ileum|Beta+Lambda - ISG score") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()

p_ileum_denovo = p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + plot_layout(ncol = 4, nrow = 2)

# Colon

p1 = DimPlot(seu$`Colon|Mock`$SeuratObject, group.by = "seurat_clusters", pt.size = 0.1, label = T, cols = c(brewer.pal(8, "Set2"),"black")) + NoLegend() + labs(title = "", subtitle = "Colon|Mock") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()
p2 = DimPlot(seu$`Colon|Beta`$SeuratObject, group.by = "seurat_clusters", pt.size = 0.1, label = T, cols = c(brewer.pal(8, "Set2"),"black")) + NoLegend() + labs(title = "", subtitle = "Colon|Beta") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()
p3 = DimPlot(seu$`Colon|Lambda`$SeuratObject, group.by = "seurat_clusters", pt.size = 0.1,label = T, cols = c(brewer.pal(8, "Set2"),"black")) + NoLegend() + labs(title = "", subtitle = "Colon|Lambda") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()
p4 = DimPlot(seu$`Colon|Beta_Lambda`$SeuratObject, group.by = "seurat_clusters", pt.size = 0.1, label = T, cols = c(brewer.pal(8, "Set2"),"black")) + NoLegend() + labs(title = "", subtitle = "Colon|Beta+Lambda") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()

p5 = FeaturePlot(seu$`Colon|Mock`$SeuratObject, features = "nFeature_RNA", pt.size = 0.1) + labs(title = "", subtitle = "Ileum|Mock - ISG score") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()
p6 = FeaturePlot(seu$`Colon|Beta`$SeuratObject, features = "nFeature_RNA", pt.size = 0.1) + labs(title = "", subtitle = "Ileum|Beta - ISG score") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()
p7 = FeaturePlot(seu$`Colon|Lambda`$SeuratObject, features = "nFeature_RNA", pt.size = 0.1) + labs(title = "", subtitle = "Ileum|Lambda - ISG score") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()
p8 = FeaturePlot(seu$`Colon|Beta_Lambda`$SeuratObject, features = "nFeature_RNA", pt.size = 0.1) + labs(title = "", subtitle = "Ileum|Beta+Lambda - ISG score") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()

p_colon_denovo = p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + plot_layout(ncol = 4, nrow = 2)

p = p_ileum_denovo / p_colon_denovo

ggsave(filename = paste0(path, "results/merged_after_kb/TreatmentSeparate/cell_types_only_using_ISGs.pdf"), plot = p, width = 9, height = 6)

rm(p1, p2, p3, p4, p5, p6, p7, p8, p_ileum_denovo, p_colon_denovo)
