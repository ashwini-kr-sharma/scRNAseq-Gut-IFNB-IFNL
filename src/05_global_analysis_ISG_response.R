library(Seurat)
library(progeny)
library(patchwork)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)

path = "/icgc/dkfzlsdf/analysis/B080/sharma/boulantLab/sc_interferon/"

#-------------------------------------------------------------------------------
# Create directories
#-------------------------------------------------------------------------------

# Create output directory 
if( ! dir.exists(paste0(path, "results/merged_after_kb/GlobalISGresponse"))){
  dir.create(paste0(path, "results/merged_after_kb/GlobalISGresponse"),recursive = T)
}

#-------------------------------------------------------------------------------
# Signature gene lists
#-------------------------------------------------------------------------------

# ISG gene signature | PMC5963606
isg = c("IFNAR1", "IFNAR2", "IFNLR1", "IL10RB", 
        "CXCL10", "DDX60", "EPSTI1", "GBP1", "HERC5", "HERC6", "IFI27", "IFI44", 
        "IFI44L", "IFI6", "IFIT1", "IFIT2", "USP18", "CGAS", "IFIT3", "IFIT5", "ISG15",
        "ISG20","LAMP3", "LY6E", "MX1", "OAS1", "OAS2", "CH25H", "OAS3", "OASL", "PLSCR1", 
        "RSAD2", "RTP4", "SIGLEC1", "SOCS1", "SPATS2L", "ISGF3")

# JAK STAT pathway genes
jak.stat = c("JAK1", "JAK2", "TYK2", 
             "STAT1",  "STAT2", "STAT3", "STAT5", 
             "IRF1", "IRF2", "IRF3", "IRF4", "IRF5", "IRF6", "IRF7", "IRF8", "IRF9")

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
# Check basal expression of INF receptors
#-------------------------------------------------------------------------------

tmp1 = lapply(seu, function(x){
  x = x@assays$RNA@counts
  x = x[which(rownames(x) %in% c(isg, jak.stat)),, drop = F]
  x = t(as.matrix(x))
  x = round(c(apply(x, 2, function(y) sum(y > 1)) / nrow(x)) * 100, 4)
})
tmp1 = unlist(tmp1)
tmp1 = tmp1[tmp1 > 0]

tmp2 = lapply(seu, function(x){
  x = x@assays$RNA@counts
  x = x[which(rownames(x) %in% c(isg, jak.stat)),, drop =F]
  x = t(as.matrix(x))
  x = apply(x, 2, function(y) mean(log2(y[y > 0])))
})
tmp2 = unlist(tmp2)
tmp2 = tmp2[names(tmp2) %in% names(tmp1)]

if(identical(names(tmp1), names(tmp2))){
  tmp = data.frame(Gene = sapply(strsplit(names(tmp1), ".", fixed=T), function(x) x[2]),
                   Condition = sapply(strsplit(names(tmp1), ".", fixed=T), function(x) x[1]),
                   Percent_expressed = tmp1, 
                   Mean_expression = as.numeric(tmp2),
                   stringsAsFactors = F)
  
}

tmp$Gene  = factor(tmp$Gene, levels = c(isg, jak.stat))
tmp = droplevels(tmp)
rm(tmp1, tmp2)

xtext.col = rep("#ff7f00", length(levels(tmp$Gene)))
xtext.col[levels(tmp$Gene) %in% c("IFNAR1", "IFNAR2", "IFNLR1", "IL10RB")] = "#377eb8"
xtext.col[levels(tmp$Gene) %in% jak.stat] = "#4daf4a"

p1 = ggplot(tmp, aes(y = Condition, x = Gene, color = Mean_expression, size = Percent_expressed)) + theme_bw(base_size = 8) + 
labs(x = "", y ="", color = "log2 mean expression", size = "Expressed in % cells") +
geom_point() + scale_colour_gradient(low="pink", high="blue", 
                                     limits = c(0, ceiling(max(tmp$Mean_expression))), 
                                     guide = guide_colourbar(title.position = "top", barheight = 0.5, ticks = F)) +
theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1, color = xtext.col), legend.position = "bottom") +
guides(size = guide_legend(title.position = "top"))

ggsave(filename = paste0(path, "results/merged_after_kb/GlobalISGresponse/ISG_response.pdf"), 
       plot = p1, width = 7, height = 2.5)

#-------------------------------------------------------------------------------
# Merge data
#-------------------------------------------------------------------------------

dat = seu
dat = merge(x = dat$`Ileum|Mock`, y = c(dat$`Ileum|Beta`, dat$`Ileum|Lambda`, dat$`Ileum|Beta_Lambda`,
                                        dat$`Colon|Mock`, dat$`Colon|Beta`, dat$`Colon|Lambda`, dat$`Colon|Beta_Lambda`))
dat[["type"]] = paste(dat@meta.data$site, dat@meta.data$treatment, sep= "|")

# Seurat pre-processing
dat <- NormalizeData(object = dat)
dat <- FindVariableFeatures(object = dat)
dat <- ScaleData(object = dat)
dat <- RunPCA(object = dat, npcs = 15)
dat <- FindNeighbors(object = dat, dim = 1:15)
dat <- FindClusters(object = dat, resolution = 0.5)
dat <- RunUMAP(object = dat, dim = 1:15)

#-------------------------------------------------------------------------------
# Progeny analysis
#-------------------------------------------------------------------------------

## Computing Progeny activity matrix
CellsClusters <- data.frame(Cell = names(Idents(dat)),
                            CellType = dat$type,
                            stringsAsFactors = FALSE)

dat <- progeny(dat, scale = FALSE, organism = "Human", top = 1000, perm = 1, return_assay = TRUE)

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

pheatmap(t(summarized_progeny_scores_cluster[,-1]), clustering_method = "ward.D2", 
         fontsize=8, treeheight_row = 15,
         color=myColor, breaks = progenyBreaks, 
         angle_col = 90,
         treeheight_col = 15,  border_color = "white",
         filename = paste0(path, "results/merged_after_kb/GlobalISGresponse/Global_Progeny_analysis.pdf"), height = 3, width = 2.5)
