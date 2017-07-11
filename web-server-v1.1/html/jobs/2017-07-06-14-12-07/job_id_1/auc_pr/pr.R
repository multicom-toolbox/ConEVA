if (!require("ggplot2")){
  install.packages("ggplot2", dependencies = TRUE)
  library(ggplot2)
}
d1 = read.table("pr.txt", sep = "	", header = F)
png("auc_pr.png", width=3000, height=1800, res=450)
p1 <- ggplot(data = d1, aes(d1$V1, d1$V2)) +
  geom_line(aes(colour = factor(d1$V3)), lwd = 0.5) +
  xlab("Recall") +
  ylab("Precision") +
  theme_bw() +
  theme(legend.position = "right", legend.title=element_blank())
print(p1)
