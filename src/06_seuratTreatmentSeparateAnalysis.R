library(Seurat)
library(CytoTRACE)
library(progeny)
library(patchwork)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(WriteXLS)
library(reshape2)
library(rio)
library(ggpubr)
library(ggplotify)

path = "/icgc/dkfzlsdf/analysis/B080/sharma/boulantLab/sc_interferon/"
#source(paste0(path,"src/merged_after_kb/modifiedDotPlot.R"))

#-------------------------------------------------------------------------------
# Create directories
#-------------------------------------------------------------------------------

# Create output directory 
if( ! dir.exists(paste0(path, "results/merged_after_kb/TreatmentSeparate"))){
    dir.create(paste0(path, "results/merged_after_kb/TreatmentSeparate"),recursive = T)
}

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
  # dat <- CreateSeuratObject(counts = dat)
  #if(ncol(dat) > 8000){dat = subset(dat, cells = sample(Cells(dat), 8000))}
  
  # Seurat preprocessing
  dat <- NormalizeData(object = dat)
  dat <- FindVariableFeatures(object = dat, nfeatures = 5000)
  dat <- ScaleData(object = dat)
  dat <- RunPCA(object = dat, npcs = 30)
  dat <- FindNeighbors(object = dat)
  dat <- FindClusters(object = dat, resolution = 0.4)
  dat <- RunUMAP(object = dat, dim = 1:15)
  
  # Find top 50 markers per cell type
  markers = FindAllMarkers(dat, logfc.threshold = 0.7, test.use = "roc")
  markers = data.frame(markers %>% group_by(cluster) %>% top_n(n = 50, wt = myAUC) %>% filter(myAUC > 0.7))
  
  return(list(SeuratObject = dat, Markers = markers)) 
}

# Read in all objects
seu = lapply(seu, preProcessSeurat)

#-------------------------------------------------------------------------------
# Signature gene lists
#-------------------------------------------------------------------------------

# ISG gene signature | PMC5963606
# isg1 = c("IFNAR1", "IFNAR2", "IFNLR1", "IL10RB", 
#         "CXCL10", "DDX60", "EPSTI1", "GBP1", "HERC5", "HERC6", "IFI27", "IFI44", 
#         "IFI44L", "IFI6", "IFIT1", "IFIT2", "USP18", "CGAS", "IFIT3", "IFIT5", "ISG15",
#         "ISG20","LAMP3", "LY6E", "MX1", "OAS1", "OAS2", "CH25H", "OAS3", "OASL", "PLSCR1", 
#         "RSAD2", "RTP4", "SIGLEC1", "SOCS1", "SPATS2L", "ISGF3")

# ISGs from https://doi.org/10.1038/s41590-019-0323-3 (Supplementary Table 1)
isg2 = import("https://static-content.springer.com/esm/art%3A10.1038%2Fs41590-019-0323-3/MediaObjects/41590_2019_323_MOESM3_ESM.xlsx", sheet = 3)
isg2 = isg2[,1:4]
colnames(isg2) = isg2[1,]
isg2 = isg2[-1,]
isg2 = isg2$ISG
#isg2 = isg2$ISG[which(isg2$`Used as Bait` == "+")]

# JAK STAT pathway genes
jak.stat = c("JAK1", "JAK2", "TYK2", 
             "STAT1",  "STAT2", "STAT3", "STAT5", 
             "IRF1", "IRF2", "IRF3", "IRF4", "IRF5", "IRF6", "IRF7", "IRF8", "IRF9")

# Ileum marker genes from Sergio et.al
markers.il = c("OLFM4", "AGR2", "LDHB", "C1QBP", "HIST1H4C", "H2AFZ",
               "MKI67", "UBE2S", "CENPF", "CCNB1", "TOP2A", "ZG16",
               "CLCA1", "SPINK4", "FCGBP", "GCG", "PYY", "NTS", "CHGA",
               "SAA2", "PDZK1IP1", "CYP3A4", "CA4", "FABP1", "PHGR1",
               "TM4SF20", "MUC13","APOA1", "FABP6", "APOB", "RBP2",
               "GUCA2A", "CCL2", "CXCL1", "SAA1", "MMP7")

# Colon marker genes from Sergio et.al
markers.co = c("GNB2L1", "EEF1B2", "EEF1A1", "AKAP12","HIST1H4C", "NPM1",
               "CENPF", "PYY", "CHGA", "TPH1", "CHGB", "NUPR1", "G0S2",
               "CCND2", "KRT18", "INSIG1", "CKB", "NEAT1", "SCD", "SLC26A3",
               "GUCA2A", "SLC40A1", "CEACAM7", "CA4", "PLCG2", "MMP7",
               "CXCL1", "IL32", "MUC1", "CLCA4","MUC13", "LAMA3", "ABCB1")

#-------------------------------------------------------------------------------
# Add ISG module expression score
#-------------------------------------------------------------------------------

for(i in 1:length(seu)){
  seu[[i]]$SeuratObject = AddModuleScore(object = seu[[i]]$SeuratObject, 
                                         features = list(isg2), name = "ISGscore")
}
rm(i)

#-------------------------------------------------------------------------------
# Add JAK-STAT module expression score
#-------------------------------------------------------------------------------

for(i in 1:length(seu)){
  seu[[i]]$SeuratObject = AddModuleScore(object = seu[[i]]$SeuratObject, 
                                         features = list(jak.stat), name = "JAK.STAT")
}
rm(i)

#-------------------------------------------------------------------------------
# Label transfer with Sergio et.al https://doi.org/10.15252/msb.202110232
#-------------------------------------------------------------------------------

# Download Sergio et.al data
if( ! file.exists(paste0(path, "data/COVID19_July.rda"))){
  download.file(url = "https://ndownloader.figshare.com/files/26328535", 
                destfile = paste0(path, "data/COVID19_July.rda"))
}

load("/icgc/dkfzlsdf/analysis/B080/sharma/boulantLab/sc_interferon/data/COVID19_July.rda")

il2d = subset(Illeum_H_T, idents = "Illeum_Mock")
tmp = CreateSeuratObject(il2d@assays$RNA@counts)
if(identical(rownames(tmp@meta.data), rownames(il2d@meta.data))){
  tmp@meta.data = cbind(tmp@meta.data, il2d@meta.data[,c("Smillie", "Ileum_Tissue", "Ileum_org", "CellTypes")])
  tmp = PercentageFeatureSet(tmp, pattern = "^MT-", col.name = "percent.mt")
  tmp@meta.data$organoid_model = "2D"
  il2d = tmp
  rm(tmp)
}

co2d = subset(Colon_H_T, idents = "Colon_Mock")
tmp = CreateSeuratObject(co2d@assays$RNA@counts)
if(identical(rownames(tmp@meta.data), rownames(co2d@meta.data))){
  tmp@meta.data = cbind(tmp@meta.data, co2d@meta.data[,c("Smillie", "Ileum_Tissue", "Ileum_org", "CellTypes")])
  tmp = PercentageFeatureSet(tmp, pattern = "^MT-", col.name = "percent.mt")
  tmp@meta.data$organoid_model = "2D"
  co2d = tmp
  rm(tmp)
}

rm(Colon_H_T, Illeum_H_T)

##########################
# Label transfer function
##########################

labelPredict = function(reference, query)
{
  reference <- NormalizeData(object = reference)
  reference <- FindVariableFeatures(object = reference, nfeatures = 5000)
  reference <- ScaleData(object = reference)
  
  # Label transfer predictions
  pred = FindTransferAnchors(reference = reference, query = query, dims = 1:30)
  pred = TransferData(anchorset = pred, refdata = reference$CellTypes, dims = 1:30)
  
  # Adding predictions to metadata
  query = AddMetaData(query, metadata = pred)
  
  return(query)
}

# Perform predictions
for(i in names(seu))
{
  if(grepl("Ileum", i)){
    seu[[i]]$SeuratObject = labelPredict(reference = il2d, query = seu[[i]]$SeuratObject)
  }else if(grepl("Colon", i)){
    seu[[i]]$SeuratObject = labelPredict(reference = co2d, query = seu[[i]]$SeuratObject)
  }
}
rm(i)

#-------------------------------------------------------------------------------
# CytoTRACE analysis
#-------------------------------------------------------------------------------

cytodat = lapply(seu, function(x){
  CytoTRACE(as.matrix(GetAssayData(x$SeuratObject, slot = "counts")))
})

cytoVal = lapply(names(seu), function(x){
  a = cytodat[[x]]$CytoTRACE
  b = seu[[x]]$SeuratObject$predicted.id
  df = NULL
  if(identical(names(a), names(b))){
    df = data.frame(cytotrace = a, clusters = b, type = x, stringsAsFactors = F)
  }
  return(df)
})

cytoVal = do.call("rbind", cytoVal)

# I have no idea why I can't do this inside the lapply above !! 
# Very weird that meta data slot does not get updated inside lapply above
seu$`Ileum|Mock`$SeuratObject$cytotrace = cytodat$`Ileum|Mock`$CytoTRACE
seu$`Ileum|Beta`$SeuratObject$cytotrace = cytodat$`Ileum|Beta`$CytoTRACE
seu$`Ileum|Lambda`$SeuratObject$cytotrace = cytodat$`Ileum|Lambda`$CytoTRACE
seu$`Ileum|Beta_Lambda`$SeuratObject$cytotrace = cytodat$`Ileum|Beta_Lambda`$CytoTRACE

seu$`Colon|Mock`$SeuratObject$cytotrace = cytodat$`Colon|Mock`$CytoTRACE
seu$`Colon|Beta`$SeuratObject$cytotrace = cytodat$`Colon|Beta`$CytoTRACE
seu$`Colon|Lambda`$SeuratObject$cytotrace = cytodat$`Colon|Lambda`$CytoTRACE 
seu$`Colon|Beta_Lambda`$SeuratObject$cytotrace = cytodat$`Colon|Beta_Lambda`$CytoTRACE 

p = ggplot(cytoVal, aes(x = clusters, y = cytotrace)) + 
  labs(x="", y="Differentiation state (1 = lowest or stem like, 0 = highest or differentiated)") + theme_bw(base_size = 9) +
  #geom_violin(trim = T, color = "black", fill = "#1b9e77", draw_quantiles = 0.5, lwd=0.3) + 
  geom_boxplot(outlier.size = 0.1, lwd=0.3) +
  facet_wrap(~type, scales = "free", ncol = 4) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour="black"), 
        axis.line = element_blank(),
        #strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        #axis.text = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "right")

ggsave(filename = paste0(path, "results/merged_after_kb/TreatmentSeparate/CytoTRACE_analysis.pdf"), plot = p, width = 6, height = 5)
rm(p)

saveRDS(object = seu, file = paste0(path, "analysis/postFilter_merge_counts/all_SeuratObjects_separate.RDS"))

#-------------------------------------------------------------------------------
# Marker gene analysis
#-------------------------------------------------------------------------------

markerlist = setNames(vector("list", length(seu)), names(seu))
for(i in names(seu)){
  
  if(length(grep("Ileum", i)) > 0){
    type = markers.il
  }else{
    type = markers.co
  }
  
  vals = seu[[i]]$Markers
  vals$isMarker = "No"
  vals$isMarker[vals$gene %in% type] = "Yes"
  vals = vals[c(8,1,4:7,9)]
  
  markerlist[[i]] = vals
  rm(type, vals)
}
rm(i)

WriteXLS(markerlist, ExcelFileName = paste0(path, "results/merged_after_kb/TreatmentSeparate/GeneMarkers.xls"),
         AdjWidth = TRUE, BoldHeaderRow = TRUE, FreezeRow = 1)

rm(markerlist)

#-------------------------------------------------------------------------------
# Plot predicted cell types
#-------------------------------------------------------------------------------

#######
# Ileum
#######
p1 = DimPlot(seu$`Ileum|Mock`$SeuratObject, group.by = "seurat_clusters", pt.size = 0.1, label = T, cols = c(brewer.pal(8, "Set2"),"black")) + NoLegend() + labs(title = "", subtitle = "Ileum|Mock") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()
p2 = DimPlot(seu$`Ileum|Beta`$SeuratObject, group.by = "seurat_clusters", pt.size = 0.1, label = T, cols = c(brewer.pal(8, "Set2"),"black")) + NoLegend() + labs(title = "", subtitle = "Ileum|Beta") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()
p3 = DimPlot(seu$`Ileum|Lambda`$SeuratObject, group.by = "seurat_clusters", pt.size = 0.1,label = T, cols = c(brewer.pal(8, "Set2"),"black")) + NoLegend() + labs(title = "", subtitle = "Ileum|Lambda") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()
p4 = DimPlot(seu$`Ileum|Beta_Lambda`$SeuratObject, group.by = "seurat_clusters", pt.size = 0.1, label = T, cols = c(brewer.pal(8, "Set2"),"black")) + NoLegend() + labs(title = "", subtitle = "Ileum|Beta+Lambda") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()

p_ileum_denovo = p1 + p2 + p3 + p4 + plot_layout(ncol = 4)

p5 = DimPlot(seu$`Ileum|Mock`$SeuratObject, group.by = "predicted.id", pt.size = 0.1, label = F, cols = c(brewer.pal(9, "Set1"),"black"))  +  labs(title = "", subtitle = "Ileum|Mock") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()
p6 = DimPlot(seu$`Ileum|Beta`$SeuratObject, group.by = "predicted.id", pt.size = 0.1, label = F, cols = c(brewer.pal(9, "Set1"),"black")) + NoLegend() +labs(title = "", subtitle = "Ileum|Beta") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()
p7 = DimPlot(seu$`Ileum|Lambda`$SeuratObject, group.by = "predicted.id", pt.size = 0.1,label = F, cols = c(brewer.pal(9, "Set1"),"black")) + NoLegend() +labs(title = "", subtitle = "Ileum|Lambda") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()
p8 = DimPlot(seu$`Ileum|Beta_Lambda`$SeuratObject, group.by = "predicted.id", pt.size = 0.1, label = F, cols = c(brewer.pal(9, "Set1"),"black")) + NoLegend() + labs(title = "", subtitle = "Ileum|Beta+Lambda") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()

p_ileum_predict = p5 + p6 + p7 + p8 + plot_layout(ncol = 4)

p9  = FeaturePlot(seu$`Ileum|Mock`$SeuratObject, features = "ISGscore1", pt.size = 0.1) + labs(title = "", subtitle = "Ileum|Mock - ISG score") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()
p10 = FeaturePlot(seu$`Ileum|Beta`$SeuratObject, features = "ISGscore1", pt.size = 0.1) + labs(title = "", subtitle = "Ileum|Beta - ISG score") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()
p11 = FeaturePlot(seu$`Ileum|Lambda`$SeuratObject, features = "ISGscore1", pt.size = 0.1) + labs(title = "", subtitle = "Ileum|Lambda - ISG score") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()
p12 = FeaturePlot(seu$`Ileum|Beta_Lambda`$SeuratObject, features = "ISGscore1", pt.size = 0.1) + labs(title = "", subtitle = "Ileum|Beta+Lambda - ISG score") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()

p_ileum_ISG = p9 + p10 + p11 + p12 + plot_layout(ncol = 4)

p13 = FeaturePlot(seu$`Ileum|Mock`$SeuratObject, features = "JAK.STAT1", pt.size = 0.1) + labs(title = "", subtitle = "Ileum|Mock - JAK STAT score") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()
p14 = FeaturePlot(seu$`Ileum|Beta`$SeuratObject, features = "JAK.STAT1", pt.size = 0.1) + labs(title = "", subtitle = "Ileum|Beta - JAK STAT score") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()
p15 = FeaturePlot(seu$`Ileum|Lambda`$SeuratObject, features = "JAK.STAT1", pt.size = 0.1) + labs(title = "", subtitle = "Ileum|Lambda - JAK STAT score") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()
p16 = FeaturePlot(seu$`Ileum|Beta_Lambda`$SeuratObject, features = "JAK.STAT1", pt.size = 0.1) + labs(title = "", subtitle = "Ileum|Beta+Lambda - JAK STAT score") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()

p_ileum_JAK.STAT = p13 + p14 + p15 + p16 + plot_layout(ncol = 4)
  
p_ileum = p_ileum_denovo / p_ileum_predict / p_ileum_ISG / p_ileum_JAK.STAT

ggsave(filename = paste0(path, "results/merged_after_kb/TreatmentSeparate/Ileum_celltypes.pdf"), plot = p_ileum, width = 9, height = 9)

rm(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p_ileum_denovo, p_ileum_predict, p_ileum_ISG, p_ileum_JAK.STAT)

########
# Colon
########
p1 = DimPlot(seu$`Colon|Mock`$SeuratObject, group.by = "seurat_clusters", pt.size = 0.1, label = T, cols = c(brewer.pal(8, "Set2"),"black")) + NoLegend() + labs(title = "", subtitle = "Colon|Mock") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()
p2 = DimPlot(seu$`Colon|Beta`$SeuratObject, group.by = "seurat_clusters", pt.size = 0.1, label = T, cols = c(brewer.pal(8, "Set2"),"black")) + NoLegend() + labs(title = "", subtitle = "Colon|Beta") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()
p3 = DimPlot(seu$`Colon|Lambda`$SeuratObject, group.by = "seurat_clusters", pt.size = 0.1,label = T, cols = c(brewer.pal(8, "Set2"),"black")) + NoLegend() + labs(title = "", subtitle = "Colon|Lambda") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()
p4 = DimPlot(seu$`Colon|Beta_Lambda`$SeuratObject, group.by = "seurat_clusters", pt.size = 0.1, label = T, cols = c(brewer.pal(8, "Set2"),"black")) + NoLegend() + labs(title = "", subtitle = "Colon|Beta+Lambda") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()

p_colon_denovo = p1 + p2 + p3 + p4 + plot_layout(ncol = 4)

p5 = DimPlot(seu$`Colon|Mock`$SeuratObject, group.by = "predicted.id", pt.size = 0.1, label = F, cols = c(brewer.pal(9, "Set1"),"black"))  +  labs(title = "", subtitle = "Colon|Mock") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()
p6 = DimPlot(seu$`Colon|Beta`$SeuratObject, group.by = "predicted.id", pt.size = 0.1, label = F, cols = c(brewer.pal(9, "Set1"),"black")) + NoLegend() +labs(title = "", subtitle = "Colon|Beta") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()
p7 = DimPlot(seu$`Colon|Lambda`$SeuratObject, group.by = "predicted.id", pt.size = 0.1,label = F, cols = c(brewer.pal(9, "Set1"),"black")) + NoLegend() +labs(title = "", subtitle = "Colon|Lambda") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()
p8 = DimPlot(seu$`Colon|Beta_Lambda`$SeuratObject, group.by = "predicted.id", pt.size = 0.1, label = F, cols = c(brewer.pal(9, "Set1"),"black")) + NoLegend() + labs(title = "", subtitle = "Colon|Beta+Lambda") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()

p_colon_predict = p5 + p6 + p7 + p8 + plot_layout(ncol = 4)

p9  = FeaturePlot(seu$`Colon|Mock`$SeuratObject, features = "ISGscore1", pt.size = 0.1) + labs(title = "", subtitle = "Colon|Mock - ISG score") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()
p10 = FeaturePlot(seu$`Colon|Beta`$SeuratObject, features = "ISGscore1",  pt.size = 0.1) + labs(title = "", subtitle = "Colon|Beta - ISG score") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()
p11 = FeaturePlot(seu$`Colon|Lambda`$SeuratObject, features = "ISGscore1",  pt.size = 0.1) + labs(title = "", subtitle = "Colon|Lambda - ISG score") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()
p12 = FeaturePlot(seu$`Colon|Beta_Lambda`$SeuratObject, features = "ISGscore1",  pt.size = 0.1) + labs(title = "", subtitle = "Colon|Beta+Lambda - ISG score") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()

p_colon_ISG = p9 + p10 + p11 + p12 + plot_layout(ncol = 4)

p13 = FeaturePlot(seu$`Colon|Mock`$SeuratObject, features = "JAK.STAT1", pt.size = 0.1) + labs(title = "", subtitle = "Colon|Mock - JAK STAT score") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()
p14 = FeaturePlot(seu$`Colon|Beta`$SeuratObject, features = "JAK.STAT1", pt.size = 0.1) + labs(title = "", subtitle = "Colon|Beta - JAK STAT score") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()
p15 = FeaturePlot(seu$`Colon|Lambda`$SeuratObject, features = "JAK.STAT1", pt.size = 0.1) + labs(title = "", subtitle = "Colon|Lambda - JAK STAT score") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()
p16 = FeaturePlot(seu$`Colon|Beta_Lambda`$SeuratObject, features = "JAK.STAT1", pt.size = 0.1) + labs(title = "", subtitle = "Colon|Beta+Lambda - JAK STAT score") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()

p_colon_JAK.STAT = p13 + p14 + p15 + p16 + plot_layout(ncol = 4)

p_colon = p_colon_denovo / p_colon_predict / p_colon_ISG / p_colon_JAK.STAT

ggsave(filename = paste0(path, "results/merged_after_kb/TreatmentSeparate/Colon_celltypes.pdf"), plot = p_colon, width = 9, height = 9)

rm(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p_colon_denovo, p_colon_predict, p_colon_ISG, p_colon_JAK.STAT)

#-------------------------------------------------------------------------------
# Plot cell composition change
#-------------------------------------------------------------------------------

prop.table = sapply(seu, function(x) as.matrix(table(x$SeuratObject$predicted.id)/ncol(x$SeuratObject)))
prop.table = melt(prop.table)
prop.table = prop.table[,c(1,3:4)]
colnames(prop.table) = c("CellTypes", "Fraction", "Treatment")
prop.table$Organoid = sapply(strsplit(prop.table$Treatment, "|", fixed = T), function(x) x[1])

prop.table.ileum = prop.table[prop.table$Organoid == "Ileum",]
prop.table.colon = prop.table[prop.table$Organoid == "Colon",]

# Ileum
p_prop_ileum1 = ggplot(prop.table.ileum, aes(x = Treatment, y = Fraction, fill = CellTypes)) + theme_bw(base_size = 8) + 
  labs(x = "", y = "% composition",  fill = "CellTypes", subtitle = "Ileum") +  coord_flip() +
  geom_bar(stat = "identity") + scale_fill_manual(values = c(brewer.pal(8, "Set1"),"black")) +
  theme(legend.position = "right", panel.grid = element_blank())

p_prop_ileum2 = ggplot(prop.table.ileum, aes(x = CellTypes, y = Fraction, fill = Treatment)) + theme_bw(base_size = 8) + 
  labs(x = "", y = "% composition",  fill = "Treatment", subtitle = "Ileum") +  coord_flip() +
  geom_bar(stat = "identity", position='dodge') + facet_wrap(~CellTypes, scales = "free") +
  theme(legend.position = "right", panel.grid = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

p_prop_ileum = p_prop_ileum1 + p_prop_ileum2 + plot_layout(nrow = 2, heights = c(0.3, 1))
ggsave(filename = paste0(path, "results/merged_after_kb/TreatmentSeparate/Ileum_celltypes_proportion_change.pdf"), plot = p_prop_ileum, width = 7, height = 7)

rm(p_prop_ileum1, p_prop_ileum2)

# Colon
p_prop_colon1 = ggplot(prop.table.colon, aes(x = Treatment, y = Fraction, fill = CellTypes)) + theme_bw(base_size = 8) + 
  labs(x = "", y = "% composition",  fill = "CellTypes", subtitle = "Colon") +  coord_flip() +
  geom_bar(stat = "identity") + scale_fill_manual(values = c(brewer.pal(8, "Set1"),"black")) +
  theme(legend.position = "right", panel.grid = element_blank())

p_prop_colon2 = ggplot(prop.table.colon, aes(x = CellTypes, y = Fraction, fill = Treatment)) + theme_bw(base_size = 8) + 
  labs(x = "", y = "% composition",  fill = "Treatment", subtitle = "Colon") +  coord_flip() +
  geom_bar(stat = "identity", position='dodge') + facet_wrap(~CellTypes, scales = "free") +
  theme(legend.position = "right", panel.grid = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

p_prop_colon = p_prop_colon1 + p_prop_colon2 + plot_layout(nrow = 2, heights = c(0.3, 1))
ggsave(filename = paste0(path, "results/merged_after_kb/TreatmentSeparate/Colon_celltypes_proportion_change.pdf"), plot = p_prop_colon, width = 7, height = 7)
rm(p_prop_colon1, p_prop_colon2)

rm(prop.table, prop.table.colon, prop.table.ileum)

#-------------------------------------------------------------------------------
# Comparison between ISG score and Stemness
#-------------------------------------------------------------------------------

ISGvsStem = setNames(vector("list", length(seu)), names(seu))
for(i in names(seu))
{
  ISGvsStem[[i]] = data.frame(ISGscore = seu[[i]]$SeuratObject$ISGscore1,
                  STEMscore = seu[[i]]$SeuratObject$cytotrace,
                  CellTypes = seu[[i]]$SeuratObject$predicted.id,
                  Treatment = i, Organoid = strsplit(i, "|", fixed =T)[[1]][1])
}
rm(i)

df = do.call("rbind", ISGvsStem)

p = ggscatter(df, x = "ISGscore", y = "STEMscore",
              color = "CellTypes", shape = 20, size = 1.3, # Points color, shape and size
              add = "reg.line",  # Add regressin line
              add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
              conf.int = TRUE, # Add confidence interval
              cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
              cor.coeff.args = list(method = "spearman", label.x = 0.2, label.y = 0.9, label.sep = "\n"),
              facet.by = c("Treatment", "Organoid"))
 p =  ggpar(p, ylim = c(-0.01,1))

ggsave(filename = paste0(path, "results/merged_after_kb/TreatmentSeparate/ISGexpression_vs_stemness_analysis.pdf"),
         plot = p, width = 10, height = 15)

#-------------------------------------------------------------------------------
# Progeny analysis
#-------------------------------------------------------------------------------

progeny_plots = lapply(names(seu), function(type){
  
  dat = seu[[type]]$SeuratObject
  
  ## Computing Progeny activity matrix
  CellsClusters <- data.frame(Cell = names(Idents(dat)),
                              CellType = dat$predicted.id,
                              Treatment = dat$treatment,
                              stringsAsFactors = FALSE)
  
  dat <- progeny(dat, scale = FALSE, organism = "Human", top = 500, perm = 1, return_assay = TRUE)
  
  ## We can now directly apply Seurat functions in our Progeny scores.
  ## For instance, we scale the pathway activity scores.
  dat <- Seurat::ScaleData(dat, assay = "progeny")
  
  ## We transform Progeny scores into a data frame to better handling the results
  progeny_scores_df <- as.data.frame(t(GetAssayData(dat, slot = "scale.data", assay = "progeny"))) %>%
    rownames_to_column("Cell") %>% gather(Pathway, Activity,-Cell)
  
  ## We match Progeny scores with the cell clusters.
  progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)
  
  ## We summarize the Progeny scores by cell population
  summarized_progeny_scores_cluster <- progeny_scores_df %>% group_by(Pathway, CellType) %>%
    summarise(avg = mean(Activity), std = sd(Activity))
  
  ## We prepare the data for the plot
  summarized_progeny_scores_cluster <- summarized_progeny_scores_cluster %>%
    dplyr::select(-std) %>%
    spread(Pathway, avg) %>%
    data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  
  ## Progeny plots
  paletteLength = 100
  myColor = colorRampPalette(c("red", "white", "Darkblue"))(paletteLength)
  
  progenyBreaks = c(seq(min(summarized_progeny_scores_cluster), 0, 
                        length.out=ceiling(paletteLength/2) + 1),
                    seq(max(summarized_progeny_scores_cluster)/paletteLength, 
                        max(summarized_progeny_scores_cluster), 
                        length.out=floor(paletteLength/2)))
  
  p = as.ggplot(pheatmap(t(summarized_progeny_scores_cluster[,-1]), clustering_method = "ward.D2", 
                 fontsize=8, treeheight_row = 15,
                 color=myColor, breaks = progenyBreaks, 
                 main = type, angle_col = 90,
                 treeheight_col = 15,  border_color = "white"))
  dev.off()
  
  return(p)
  
})

p = progeny_plots[[1]] + progeny_plots[[2]] + progeny_plots[[3]] + progeny_plots[[4]] + 
    progeny_plots[[5]] + progeny_plots[[6]] + progeny_plots[[7]] + progeny_plots[[8]] + plot_layout(ncol = 4, nrow = 2)

ggsave(filename = paste0(path, "results/merged_after_kb/TreatmentSeparate/Progeny_analysis_per_treatment.pdf"), plot = p, height = 6.5, width = 11)

rm(p)
