library(Seurat)

path = "/icgc/dkfzlsdf/analysis/B080/sharma/boulantLab/sc_interferon/"

#-------------------------------------------------------------------------------
# Create directories
#-------------------------------------------------------------------------------

# Create output directory 
if( ! dir.exists(paste0(path, "analysis/postFilter_merge_counts/Seurat"))){
  dir.create(paste0(path, "analysis/postFilter_merge_counts/Seurat"),recursive = T)
}

f1 = list.files(paste0(path, "analysis/HJHG3BGXF_counts/Seurat"), full.names = T)
f2 = list.files(paste0(path, "analysis/HJJCYBGXF_counts/Seurat"), full.names = T)

for(i in paste0("Int", 1:8))
{
  a = readRDS(grep(i, f1, value =T))
  b = readRDS(grep(i, f2, value =T))
  m = merge(a, b)
  
  saveRDS(object = m, file = paste0(path, "analysis/postFilter_merge_counts/Seurat/", i, "_filterCounts.RDS"))
  
  rm(a,b, m)
}
rm(i)
