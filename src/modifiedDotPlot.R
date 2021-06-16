modifiedDotPlot = function(dat, selgenes, type){
  
  mat = dat@assays$RNA@counts
  mat = mat[which(rownames(mat) %in% selgenes),, drop = F]
  mat = data.frame(t(as.matrix(mat)))
  mat = split(mat, paste0("Cluster-",dat$seurat_clusters))
  
  # Percentage cells expressed
  tmp1 = lapply(mat, function(x){
    round(c(apply(x, 2, function(y) sum(y > 1)) / nrow(x)) * 100, 4)
  })
  tmp1 = unlist(tmp1)
  tmp1 = tmp1[tmp1 > 0]
  
  # Average expression
  tmp2 = lapply(mat, function(x){
    apply(x, 2, function(y) mean(log2(y[y > 0])))
  })
  tmp2 = unlist(tmp2)
  tmp2 = tmp2[names(tmp2) %in% names(tmp1)]
  
  # Merging percent expressed and average expression
  if(identical(names(tmp1), names(tmp2))){
    tmp = data.frame(Gene = sapply(strsplit(names(tmp1), ".", fixed=T), function(x) x[2]),
                     Condition = sapply(strsplit(names(tmp1), ".", fixed=T), function(x) x[1]),
                     Percent_expressed = tmp1, 
                     Mean_expression = as.numeric(tmp2),
                     stringsAsFactors = F)}
  
  tmp$Gene = factor(tmp$Gene, levels = selgenes)
  tmp = droplevels(tmp)
  
  xtext.col = rep("#ff7f00", length(levels(tmp$Gene)))
  xtext.col[levels(tmp$Gene) %in% c("IFNAR1", "IFNAR2", "IFNLR1", "IL10RB")] = "#377eb8"
  xtext.col[levels(tmp$Gene) %in% jak.stat] = "#4daf4a"
  
  p = ggplot(tmp, aes(y = Condition, x = Gene, color = Mean_expression, size = Percent_expressed)) + theme_bw(base_size = 8) + 
    labs(x = "", y ="", subtitle = type, color = "log2 mean expression", size = "Expressed in % cells") +
    geom_point() + scale_colour_gradient(low="pink", high="blue", 
                                         limits = c(0, ceiling(max(tmp$Mean_expression))), 
                                         guide = guide_colourbar(title.position = "top", barheight = 0.5, ticks = F)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1, color = xtext.col), legend.position = "bottom") +
    guides(size = guide_legend(title.position = "top"))
  
  return(p)
}
