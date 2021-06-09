path = "/icgc/dkfzlsdf/analysis/B080/sharma/boulantLab/sc_interferon/"

# Information from Megan Stanifer
df = data.frame(rbind( 
                       c("Int1", "Ileum", "Mock"),
                       c("Int2", "Ileum", "Beta"),
                       c("Int3", "Ileum", "Lambda"),
                       c("Int4", "Ileum", "Beta_Lambda"),
                       c("Int5", "Colon", "Mock"),
                       c("Int6", "Colon", "Beta"),
                       c("Int7", "Colon", "Lambda"),
                       c("Int8", "Colon", "Beta_Lambda")
                     ),
                stringsAsFactors = F)
colnames(df) = c("type", "site", "treatment")

# Technical replicate 1
f1  = list.files(paste0(path, "data/sc_Interferon/Int_HJJCYBGXF"), full.names = T)
id1 = gsub("_sequence.txt.gz", "", gsub("HJJCYBGXF_Pool_Interferon_20s001227-1-1_Triana_lane1", "", basename(f1)))
id1 = do.call("rbind", strsplit(id1, "_", fixed = T))
f1  = cbind(files = f1, type = id1[,1])
f1  = merge(x = f1, y = df, all = T)
f1  = data.frame(f1[,c(2, 1, 3,4)], replicate = "HJJCYBGXF", stringsAsFactors = F)
saveRDS(f1, paste0(path, "data/scINF_HJJCYBGXF_anno.RDS"))
rm(id1)
#split(f1, paste(f1$site,f1$treatment,sep="|"))

# Technical replicate 2
f2  = list.files(paste0(path, "/data/sc_Interferon/Int_HJHG3BGXF"), full.names = T)
id2 = gsub("_sequence.txt.gz", "", gsub("HJHG3BGXF_Pool_Interferon_20s001227-2-1_Triana_lane1", "", basename(f2)))
id2 = do.call("rbind", strsplit(id2, "_", fixed = T))
f2  = cbind(files = f2, type = id2[,1])
f2  = merge(x = f2, y = df, all = T)
f2  = data.frame(f2[,c(2, 1, 3,4)], replicate = "HJHG3BGXF", stringsAsFactors = F)
saveRDS(f2, paste0(path, "data/scINF_HJHG3BGXF_anno.RDS"))
rm(id2, df)

#split(f2, paste(f2$site,f2$treatment,sep="|"))
