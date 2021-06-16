# scRNAseq-Gut-IFNB-IFNL
Ileum and Colon organoids were derived from small intestine biopsies from the same patient (i.e. same genetic background). The organoids were then treated with either IFN-B or IFN-L or both followed by single-cell RNA sequencing.

## Experimental design

![experimental design](/ExpDesign_scRNAseq-Gut-IFNB-IFNL.png)

## Scripts

- `src/01_makeAnnotation.R` - Scripts for generating the file annotations
- `src/02_kallistoAlignment.sh` - Bash script to run multiple kallisto alignments using file `02_runKB`
- `src/02_runKB.sh` - Bash script to run kallisto|bustools on the scRNAseq data
- `src/03_seuratPreProcess.R` - Script to perform standard filtering of scRNAseq data
- `src/04_techReplicateMerged.R` - Script to merge technical replicates
- `src/05_global_analysis_ISG_response.R` - Scripts to analyse global expression of ISG's, treating the scRNAseq data as pseudo bulk
- `src/06_seuratTreatmentSeparateAnalysis.R` - Standard Seurat workflow and Label transfer for individual samples
- `src/06_seuratTreatmentSeparateAnalysis_onlyISGs.R` - Standard Seurat workflow subsetting for only ISGs
- `src/07_seuratIntegratedAnalysis.R` - Standard Seurat workflow and Label transfer on integrated samples
- `src/08_CellChat.R` - Script for analysing cell-cell communication
- `src/modifiedDotPlot.R` - My own modified version of DotPlot, that I use in some situations.
