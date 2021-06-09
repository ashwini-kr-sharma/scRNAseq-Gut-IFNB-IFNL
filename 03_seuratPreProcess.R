#-------------------------------------------------------------------------------
# Install required libraries and set the root path
#-------------------------------------------------------------------------------

library(DropletUtils)
library(Matrix)
library(tidyverse)
library(Seurat)
library(scDblFinder)
library(scico)
library(reshape2)
library(patchwork)

path = "/icgc/dkfzlsdf/analysis/B080/sharma/boulantLab/sc_interferon/"

#-------------------------------------------------------------------------------
# Run analysis over all samples and replicates
#-------------------------------------------------------------------------------

for (rep in c("HJHG3BGXF", "HJJCYBGXF"))
{
  
  # Reading the corresponding annotation file
  anno = readRDS(paste0(path,"data/", "scINF_",rep,"_anno.RDS"))[,2:5]
  anno = unique(anno)
  
  for(type in paste0("Int",1:8))
  {
    print(paste0("Running analysis for - ", rep, " and type - ", type))
    
    #---------------------------------------------------------------------------
    # Create output directory 
    #---------------------------------------------------------------------------
    
    if( ! dir.exists(paste0(path, "results/merged_after_kb/QC/", rep))){
      dir.create(paste0(path, "results/merged_after_kb/QC/", rep), recursive = T)
    }
    
    if( ! dir.exists(paste0(path, "analysis/", rep, "_counts/Seurat"))){
      dir.create(paste0(path, "analysis/", rep, "_counts/Seurat"), recursive = T)
    }
    
    #---------------------------------------------------------------------------
    # Function to read kb-python output
    #---------------------------------------------------------------------------
    
    read_count_output <- function(dir, name) {
      dir <- normalizePath(dir, mustWork = TRUE)
      m <- readMM(paste0(dir, "/", name, ".mtx"))
      m <- Matrix::t(m)
      m <- as(m, "dgCMatrix")
      # The matrix read has cells in rows
      ge <- ".genes.txt"
      genes <- readLines(file(paste0(dir, "/", name, ge)))
      barcodes <- readLines(file(paste0(dir, "/", name, ".barcodes.txt")))
      colnames(m) <- barcodes
      rownames(m) <- genes
      return(m)
    }
    
    res_mat <- read_count_output(dir = paste0(path, "analysis/", rep, "_counts/kb/", type, "/counts_unfiltered"), name = "cells_x_genes")
    
    #---------------------------------------------------------------------------
    # 1. QC plot - Knee plot function
    #---------------------------------------------------------------------------
    
    knee_plot <- function(bc_rank) {
                  knee_plt <- tibble(rank = bc_rank[["rank"]],
                                     total = bc_rank[["total"]]) %>% distinct() %>% dplyr::filter(total > 0)
                  
                  annot <- tibble(inflection = metadata(bc_rank)[["inflection"]],
                                  rank_cutoff = max(bc_rank$rank[bc_rank$total > metadata(bc_rank)[["inflection"]]]))
                  
                  p <- ggplot(knee_plt, aes(total+1, rank)) + theme_bw(base_size = 8) +
                      geom_line() +
                      geom_hline(aes(yintercept = rank_cutoff), data = annot, linetype = 2) +
                      geom_vline(aes(xintercept = inflection), data = annot, linetype = 2) +
                      scale_x_log10() +
                      scale_y_log10() +
                      annotation_logticks() +
                      labs(y = "Rank", x = "Total UMIs")
                  return(p)
    }
    
    bc_rank <- barcodeRanks(res_mat, lower = 300)
    
    p1 = knee_plot(bc_rank)
    
    #---------------------------------------------------------------------------
    # 2. QC plot - nCount vs nGene per cell
    #---------------------------------------------------------------------------
    
    df2 <- tibble(nCount = colSums(res_mat),
                  nGene = colSums(res_mat > 0))
    
    p2 =  ggplot(df2, aes(nCount+1, nGene+1)) + theme_bw(base_size = 8) +
          geom_bin2d(bins = 50) +
          scale_fill_scico(palette = "devon", direction = -1, end = 0.95) +
          scale_x_log10() + scale_y_log10() + annotation_logticks() +
          geom_vline(xintercept = metadata(bc_rank)[["inflection"]]) +
          geom_hline(yintercept = 500) +
          labs(x = "Total UMI counts", y = "Number of genes detected") +
          theme(legend.position = "bottom",
                legend.text = element_text(angle = 45), 
                legend.text.align = 1) + guides(fill = guide_colourbar(title.position="left", title.vjust = 0.75))
        
    #---------------------------------------------------------------------------
    # 3. QC plot - distribution of nGene per cell
    #---------------------------------------------------------------------------
  
    p3 = ggplot(df2, aes(nGene+1)) + theme_bw(base_size = 8) + scale_x_log10() +
         #geom_vline(xintercept = 100) + 
         labs(y = "Density") +
         geom_vline(xintercept = 500) +
         annotation_logticks() + geom_density()
    
    #---------------------------------------------------------------------------
    # Filtering empty barcodes based on -
    #
    # 1. Total UMI count in each cell > knee plot inflection
    # 2. Number of genes detected in each cell > 300
    #
    # [After step (1 + 2)]
    #
    # 3. Genes detected in at least 1% of non-empty cells 
    #---------------------------------------------------------------------------
    
    res_mat <- res_mat[,colSums(res_mat) > metadata(bc_rank)$inflection & colSums(res_mat > 0) > 500]
    res_mat <- res_mat[rowSums(res_mat > 0) > 0.01 * ncol(res_mat),]
    
    rm(df2, bc_rank)
    
    # #-------------------------------------------------------------------------------
    # # Identify empty barcodes
    # #-------------------------------------------------------------------------------
    # 
    # # Identify empty droplets
    # # Any barcode with total UMI count < 100 will be automatically considered empty and assigned a NA value
    # set.seed(123)
    # e.out = emptyDrops(res_mat, lower = 100)
    # is.cell = e.out$FDR <= 0.01
    # 
    # # Check to ensure that the number of iteration in the MC simulation was sufficient or not
    # # See http://bioconductor.org/books/release/OSCA/droplet-processing.html
    # table(Significant=is.cell, Limited=e.out$Limited)
    # 
    # # Should Monte Carlo simulation iterations be increased (if 0 no, if > 0 yes)
    # mc = table(Significant=is.cell, Limited=e.out$Limited)[1,2]
    # 
    # # Empty drops selected empty cells
    # ed = sum(is.cell, na.rm=TRUE) / length(is.cell) * 100
    # 
    # # Low expression based empty cell (if total UMI count per cell < 100)
    # # Empty drops assigns such cell NA values
    # le = sum(is.na(is.cell)) / length(is.cell) * 100
    # 
    # # Non empty cells
    # ne = sum(!is.cell, na.rm=TRUE) / length(is.cell) * 100
    # 
    # # Total number of barcodes
    # bc = length(is.cell)
    # 
    # # Filter empty barcodes
    # res_mat = res_mat[,which(e.out$FDR <= 0.01)]
    # 
    # # Saving empty droplets QC
    # saveRDS(setNames(c(mc, ed, le, ne, bc), c("monte_carlo_sim", "empty_drops", "low_expr", "non_empty", "barcode_count")),
    #         paste0(path, "analysis/", rep,"_QC/",type,"/emptydroplets/edStats.RDS"))
    # saveRDS(e.out, paste0(path, "analysis/", rep,"_QC/",type,"/emptydroplets/emptyDropsOutput.RDS"))
    # 
    # rm(e.out, is.cell, mc, ed, le, ne, bc)
    # 
  
    #-------------------------------------------------------------------------------
    # Gene name mapping
    #-------------------------------------------------------------------------------
    
    tr2g <- read_tsv(paste0(path, "data/kallistoIndex/t2g.txt"), 
                     col_names = c("transcript", "gene", "gene_symbol")) %>% select(-transcript) %>% distinct()
    rownames(res_mat) <- tr2g$gene_symbol[match(rownames(res_mat), tr2g$gene)]
    rm(tr2g)
    
    #-------------------------------------------------------------------------------
    # Seurat pre-processing
    #-------------------------------------------------------------------------------
    
    # Make a Seurat object, with the following filtering criteria -
    # removing cells with less than 100 genes 
    # removing genes expressed in less than 100 cells
    seu = CreateSeuratObject(res_mat, project = paste0("scINF_", rep, "_", type))
    rm(res_mat)
    
    # Additional study annotation 
    seu[["site"]] = anno$site[anno$type == type]
    seu[["treatment"]] = anno$treatment[anno$type == type]
    seu[["replicate"]] = anno$replicate[anno$type == type]
    
    # Estimate mitochondrial composition per barcode
    seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
    
    # Identify doublets
    dbl = scDblFinder(as.SingleCellExperiment(seu))
    seu[["doublets"]] <- dbl$scDblFinder.class

    dbl = data.frame(table(dbl$scDblFinder.class))
    p4  = ggplot(data=dbl, aes(x=Var1, y=Freq)) + theme_bw(base_size = 8) +
          geom_bar(stat='identity', position="dodge") + labs(x="", y = "Counts",  fill = "") +
          geom_text(aes(label = Freq), position = position_dodge(width=0.9), vjust = 0, size = 3)
    rm(dbl)
    
    #---------------------------------------------------------------------------
    # Final consolidated QC plots
    #---------------------------------------------------------------------------
    
    qcdat = melt(seu@meta.data)
    qcdat = data.frame(Treatment=paste(qcdat$site, qcdat$treatment, sep= "-"),
                       Type = qcdat$variable, value = qcdat$value, stringsAsFactors = F)
    
    p5 = ggplot(qcdat, aes(x = Type, y = value)) + labs(x="", y="") +
      theme_bw(base_size = 8) + geom_violin(color = "grey60", draw_quantiles = 0.5, lwd=0.3) + 
      facet_wrap(~Type, scales = "free") + theme(strip.background = element_blank(), strip.text = element_blank())
    
    p6 = ggplot(seu@meta.data, aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) + 
      theme_bw(base_size = 8) +  labs(x = "Total RNA counts", y = "Number of genes detected") + 
      geom_point(size = 0.2) + scale_shape_manual(values = c(1,3)) +
      scale_x_log10() + scale_y_log10() + annotation_logticks() +
      geom_vline(xintercept = quantile(seu@meta.data$nCount_RNA, 0.99)) +
      geom_hline(yintercept = quantile(seu@meta.data$nFeature_RNA, 0.99)) +
      facet_wrap(~doublets, scales = "free") +  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    
    # Integrated QC plot
    p1 = p1 + labs(subtitle = paste(type, unique(qcdat$Treatment),sep="-"))
    p = (p1 + p2 + p3) / (p4 + p5 + p6)
    
    ggsave(filename = paste0(paste0(path, "results/merged_after_kb/QC/", rep, "/"), 
                             paste(type, unique(qcdat$Treatment), "QCresults.pdf", sep="-")
                             ), plot = p, width = 7.5, height = 6)
    
    rm(p, p1, p2, p3, p4, p5, p6, qcdat)
    
    #---------------------------------------------------------------------------
    # Filtering poor quality cells, i.e
    #
    # 1. With MT content > 20%
    # 2. Cells classified as doublets
    #---------------------------------------------------------------------------
    
    seufilt <- subset(seu, subset = percent.mt < 20)
    seufilt <- subset(seufilt, subset = doublets == "singlet")
    
    # Saving filtered objects
    saveRDS(seufilt, paste0(path, "analysis/", rep,"_counts/Seurat/",type,"_filteredCounts.RDS"))
    rm(seu, seufilt)
  }
  rm(anno)
}
