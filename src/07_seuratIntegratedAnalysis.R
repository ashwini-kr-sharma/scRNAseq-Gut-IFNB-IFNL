#-------------------------------------------------------------------------------
# Install required libraries and set the root path
#-------------------------------------------------------------------------------

library(Seurat)
library(CytoTRACE)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(reshape2)
library(progeny)
library(dorothea)
library(tidyverse)
library(pheatmap)
library(rio)

path = "/icgc/dkfzlsdf/analysis/B080/sharma/boulantLab/sc_interferon/"
seu = readRDS(paste0(path, "analysis/postFilter_merge_counts/all_SeuratObjects_separate.RDS"))
seu = lapply(seu, function(x) x$SeuratObject)

#-------------------------------------------------------------------------------
# Create output directory 
#-------------------------------------------------------------------------------

# Create output directory 
if( ! dir.exists(paste0(path, "results/merged_after_kb/TreatmentIntegrated"))){
  dir.create(paste0(path, "results/merged_after_kb/TreatmentIntegrated"),recursive = T)
}

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
# Label transfer with Sergio et.al https://doi.org/10.15252/msb.202110232
#-------------------------------------------------------------------------------

# Create output directory 
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

#-------------------------------------------------------------------------------
# Function for Data integration
#-------------------------------------------------------------------------------

SeuratIntegrate = function(tissue)
{
  dat = seu[grep(tissue, names(seu))]
  
  if(tissue == "Ileum"){
    ref = il2d
  }else if(tissue == "Colon"){
    ref = co2d
  }
  
  features <- SelectIntegrationFeatures(object.list = dat)
  anchors  <- FindIntegrationAnchors(object.list = dat, anchor.features = features)
  combined <- IntegrateData(anchorset = anchors)
    
  DefaultAssay(combined) <- "integrated"
  
  # Run the standard workflow for visualization and clustering
  combined <- ScaleData(combined, verbose = FALSE)
  combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
  combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)
  combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30)
  combined <- FindClusters(combined, resolution = 0.4)
  
  combined = labelPredict(reference = ref, query = combined)
  
  # Add ISG module score
  DefaultAssay(combined) <- "RNA"
  combined = AddModuleScore(object = combined, features = list(isg2), name = "ISGscore")
  
  # Add JAK.STAT module score
  combined = AddModuleScore(object = combined, features = list(jak.stat), name = "JAK.STAT")
  
  # Add CytoTRACE stemness scores
  cytodat = CytoTRACE(as.matrix(GetAssayData(combined, slot = "counts")))
  
  combined@meta.data$cytotrace = cytodat$CytoTRACE
  
  # Changing back the default assay
  DefaultAssay(combined) <- "integrated"
  
  return(combined)
}

il.int = SeuratIntegrate(tissue = "Ileum")
co.int = SeuratIntegrate(tissue = "Colon")

saveRDS(object = list(IleumCIntegrated = il.int, ColonIntegrated = co.int), 
        file = paste0(path, "analysis/postFilter_merge_counts/all_SeuratObjects_Integrated.RDS"))

# Plotting integrated data
plotCombined = function(combined, tissue)
{
  p1 <- DimPlot(combined, reduction = "umap", pt.size = 0.1, label = T) + 
    NoLegend() + labs(title = "", subtitle = tissue) + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()
  
  p2 <- FeaturePlot(combined, features = "ISGscore1", pt.size = 0.1) + labs(title = "", subtitle = "ISG score") + 
    theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()
  
  p3 <- DimPlot(combined, reduction = "umap", pt.size = 0.1, cols = brewer.pal(9, "Set1"), group.by = "predicted.id") + 
    labs(title = "") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()
  
  p4 <- DimPlot(combined, reduction = "umap", pt.size = 0.1, cols = brewer.pal(9, "Set1"), split.by = "treatment", group.by = "predicted.id") + 
     labs(title = "") + theme(text = element_text(size = 8)) + FontSize(8) + NoAxes()
  
  a = as.matrix(table(combined$treatment, combined$predicted.id)) 
  b = c(table(combined$treatment))
  m = c(table(combined$predicted.id))
  
  prop_1 = round(a/b * 100, 2)
  prop_1 = melt(prop_1)
  prop_1$Var2 = factor(prop_1$Var2, levels = unique(prop_1$Var2))
  
  prop_2 = round(t(a)/m * 100, 2)
  prop_2 = reshape2::melt(prop_2)
  prop_2$Var1 = factor(prop_2$Var1, levels = unique(prop_2$Var1))
  
  p5 = ggplot(prop_2, aes(x = Var1, y = value, fill = Var2)) + theme_bw(base_size = 8) + 
    labs(x = "Clusters", y = "Percentage (%)", fill = "Treatment") + coord_flip(expand = F) +
    geom_bar(stat = "identity") + scale_fill_manual(values = brewer.pal(8, "Set2")) + 
    theme(legend.position = "bottom", panel.grid = element_blank()) #+ guides(fill = F)
  
  p6 = ggplot(prop_1, aes(x = Var1, y = value, fill = Var2)) + theme_bw(base_size = 8) + 
    labs(x = "", y = "Percentage (%)", fill = "Clusters") +  coord_flip(expand = F) + 
    geom_bar(stat = "identity") + scale_fill_manual(values = brewer.pal(9, "Set1")) + 
    theme(legend.position = "bottom", panel.grid = element_blank()) #+ guides(fill = F)
  
  #p = p + p5 / p6
  p = (p1 | p2 | p3) / p4 / p5 / p6
  
  ggsave(filename = paste0(path, "results/merged_after_kb/TreatmentIntegrated/", tissue, "_cell_types.pdf"), plot = p, width = 10, height = 10)
  
  return(NULL)
}

plotCombined(combined = il.int, tissue = "Ileum")
plotCombined(combined = co.int, tissue = "Colon")

# #-------------------------------------------------------------------------------
# # Function to plot marker gene expression in the integrated data
# #-------------------------------------------------------------------------------
# 
# source(paste0(path,"src/modifiedDotPlot.R"))
# 
# MarkerGenePlot = function(study, type, dat.combined)
# {
#   if(type == "Ileum")
#   {
#     sergio.markers = markers.il
#   }else if(type == "Colon"){
#     sergio.markers = markers.co
#   }
#   
#   p1 = modifiedDotPlot(dat = dat.combined, selgenes = sergio.markers, type = type)
#   p2 = modifiedDotPlot(dat = dat.combined, selgenes = c(isg, jak.stat), type = type)
#   
#   ggsave(filename = paste0(path, "/results/TreatmentIntegrated/", study, "_", type, "_integrated_markergenes_sergio.pdf"), plot = p1, width = 7, height = 4.5)
#   ggsave(filename = paste0(path, "/results/TreatmentIntegrated/", study, "_", type, "_integrated_markergenes_isg_jak_stat.pdf"), plot = p2, width = 7, height = 4.5)
#   
#   return(NULL)
# }
# 
# MarkerGenePlot(study = "HJHG3BGXF", type = "Ileum", dat.combined = il.int)
# MarkerGenePlot(study = "HJHG3BGXF", type = "Colon", dat.combined = co.int)
# 
# #-------------------------------------------------------------------------------
# # Progeny analysis on the integrated data
# #-------------------------------------------------------------------------------
# 
# progenyIntegrated = function(study, type, dat.combined)
# {
#     DefaultAssay(dat.combined) = "RNA"
#     
#     ## Computing Progeny activity matrix
#     CellsClusters <- data.frame(Cell = names(Idents(dat.combined)),
#                                 CellType = dat.combined$seurat_clusters,
#                                 Treatment = dat.combined$treatment,
#                                 stringsAsFactors = FALSE)
#     
#     dat.combined <- progeny(dat.combined, scale = FALSE, organism = "Human", top = 500, perm = 1, return_assay = TRUE)
#     
#     ## We can now directly apply Seurat functions in our Progeny scores.
#     ## For instance, we scale the pathway activity scores.
#     dat.combined <- Seurat::ScaleData(dat.combined, assay = "progeny")
#     
#     ## We transform Progeny scores into a data frame to better handling the results
#     progeny_scores_df <- as.data.frame(t(GetAssayData(dat.combined, slot = "scale.data", assay = "progeny"))) %>%
#       rownames_to_column("Cell") %>% gather(Pathway, Activity,-Cell)
#     
#     ## We match Progeny scores with the cell clusters.
#     progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)
#     
#     ################################
#     # Summarize by Seurat clusters
#     ################################
#     
#     summarized_progeny_scores_cluster <- progeny_scores_df %>% group_by(Pathway, CellType) %>%
#       summarise(avg = mean(Activity), std = sd(Activity))
#     
#     summarized_progeny_scores_cluster <- summarized_progeny_scores_cluster %>%
#       dplyr::select(-std) %>%
#       spread(Pathway, avg) %>%
#       data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
#     
#     ################################
#     # Summarize by Treatment
#     ################################
#     
#     summarized_progeny_scores_treatment <- progeny_scores_df %>% group_by(Pathway, Treatment) %>%
#       summarise(avg = mean(Activity), std = sd(Activity))
#     
#     summarized_progeny_scores_treatment <- summarized_progeny_scores_treatment %>%
#       dplyr::select(-std) %>%
#       spread(Pathway, avg) %>%
#       data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
#     
#     ## Progeny plots
#     paletteLength = 100
#     myColor = colorRampPalette(c("red", "white", "Darkblue"))(paletteLength)
#     
#     progenyBreaks.clusters = c(seq(min(summarized_progeny_scores_cluster), 0, 
#                               length.out=ceiling(paletteLength/2) + 1),
#                               seq(max(summarized_progeny_scores_cluster)/paletteLength, 
#                               max(summarized_progeny_scores_cluster), 
#                               length.out=floor(paletteLength/2)))
#     
#     progenyBreaks.treatment = c(seq(min(summarized_progeny_scores_treatment), 0, 
#                                 length.out=ceiling(paletteLength/2) + 1),
#                                 seq(max(summarized_progeny_scores_treatment)/paletteLength, 
#                                 max(summarized_progeny_scores_treatment), 
#                                 length.out=floor(paletteLength/2)))
#     
#     pheatmap(t(summarized_progeny_scores_cluster[,-1]), clustering_method = "ward.D2", 
#              fontsize=8, treeheight_row = 15,
#              color=myColor, breaks = progenyBreaks.clusters, 
#              main = type, angle_col = 45,
#              treeheight_col = 15,  border_color = "white",
#              filename = paste0(path, "results/TreatmentIntegrated/", study, "_", type, "_Progeny_clusters.pdf"), height = 4, width = 4)
#     
#     pheatmap(t(summarized_progeny_scores_treatment[,-1]), clustering_method = "ward.D2", 
#              fontsize=8, treeheight_row = 15,
#              color=myColor, breaks = progenyBreaks.treatment, 
#              main = type, angle_col = 45,
#              treeheight_col = 15,  border_color = "white",
#              filename = paste0(path, "results/TreatmentIntegrated/", study, "_", type, "_Progeny_treatment.pdf"), height = 4, width = 2.5)
#     
#     return(NULL)
# }
# 
# progenyIntegrated(study = "HJHG3BGXF", type = "Ileum", dat.combined = il.int)
# progenyIntegrated(study = "HJHG3BGXF", type = "Colon", dat.combined = co.int)
# 
