png("2L_long.png", width=1600, height=1600, res=300)
infile = "temp.rr"
con <- file(infile, "rt")
seqlen <- readLines(con, 1)
L <- nchar(seqlen)
data <- read.table(infile, skip = 1, header = FALSE)
col_count <- ncol(data)
data$V7 <- abs(data$V1-data$V2)
data$V7[data$V7 > 23] <- 200
data$V7[data$V7 < 12] <- 100
data$V7[data$V7 < 24] <- 150
data$V7[data$V7 == 100] <- 4
data$V7[data$V7 == 150] <- 3
data$V7[data$V7 == 200] <- 2
data$V8[data$V7 == 4] <- "Short-Range"
data$V8[data$V7 == 3] <- "Medium-Range"
data$V8[data$V7 == 2] <- "Long-Range"
data$V9 = as.numeric(data$V6)
par(mar=c(2,2,0.5,0.5))
# Plot for an RR file with multiple sources
native <- subset(data, data$V6 == "PDB")
plot(native$V1, native$V2, col=rgb(0,0,0,0.1), pch=0, xlab = NULL, ylab = NULL, ylim=c(L, 1), xlim=c(1, L))
predic <- subset(data, data$V6 != "PDB")
points(predic$V1, predic$V2, col=predic$V6, xlab = NULL, ylab = NULL, ylim=c(L, 1), xlim=c(1, L), pch = predic$V9)
legend("topright", bty="n", legend=unique(predic$V6), pch = unique(predic$V9), col=unique(predic$V9))
text(0, 0, "Top-2L", cex=2.0, adj = c(0,1))
