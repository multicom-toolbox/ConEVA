if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}
png("L5.png", width=2000, height=2000, res=300)
table_jaccard = read.table(file="mat_topL5.txt")
jacc_matrix = data.matrix(table_jaccard)
heatmap.2(jacc_matrix, tracecol = NA, symm = T, dendrogram = "col", Colv = T, scale = 'none', col=colorRampPalette(brewer.pal(8,"YlOrRd"))(length(jacc_matrix)), key.ylab = NA, key.xlab = "Jaccard Similarity", lmat = rbind(c(3,4),c(1,2)), lwid = c(4,1), lhei = c(1,4), density.info="none", margins = c(12, 12), srtCol=90, main="Top-L/5")
