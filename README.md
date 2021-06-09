# scRNAseq-Gut-IFNB-IFNL
Ileum and Colon organoids were derived from small intestine biopsies from the same patient (i.e. same genetic background). The organoids were then treated with either IFN-B or IFN-L or both followed by single-cell RNA sequencing.

## Experimental design

![experimental design](/ExpDesign_scRNAseq-Gut-IFNB-IFNL.png)

## Scripts

- `01_makeAnnotation.R` - Scripts for generating the file annotations
- `02_kallistoAlignment.sh` - Bash script to run multiple kallisto alignments using file `02_runKB`
- `02_runKB.sh` - Bash script to run kallisto|bustools on the scRNAseq data
- `03_seuratPreProcess.R` - Script to perform standard filtering of scRNAseq data
- `04_techReplicateMerged.R` - Script to merge technical replicates
- `05_global_analysis_ISG_response.R` - Scripts to analyse global expression of ISG's, treating the scRNAseq data as pseudo bulk
- `06_seuratTreatmentSeparateAnalysis.R` - Standard Seurat workflow and Label transfer for individual samples
- `06_seuratTreatmentSeparateAnalysis_onlyISGs.R` - Standard Seurat workflow subsetting for only ISGs
- `07_seuratIntegratedAnalysis.R` - Standard Seurat workflow and Label transfer on integrated samples
- `08_CellChat.R` - Script for analysing cell-cell communication
- `modifiedDotPlot.R` - My own modified version of DotPlot, that I use in some situations.
