png("coverage_long.png", width=2500, height=1400, res=300)
par(mar=c(4.5,4.5,0.2,11), xpd = TRUE)
infile = "coverage_long.txt"
original <- read.table(file=infile, header=F)
sorted <- original[order(original$V1, original$V2),]
y_range <- range(0, original$V3)
sources <- unique(original$V1)
start_flag = 1
src_id = 1
for(src in sources){
	selected <- original$V3[which(original$V1 == src)]
	if(start_flag == 1){
		plot(selected, col=src_id, type="l", lwd = 2, xaxt='n', ylim=y_range, ylab="Coverage", xlab="Top Ranked Contacts")
		start_flag = 0
	}
	else{
		lines(selected, col=src_id, lty=src_id,  lwd = 3)
	}
	src_id = src_id + 1
}
axis(1, at=1:6, labels=c("Top5","L/10","L/5","L/2","L","2L"))
if (length(sources) > 1){
	legend("topright", inset=c(-0.4,0), legend = sources, cex=1.0, lwd = 3, col = seq(1, length(sources), by = 1), lty = seq(1, length(sources), by = 1))
}
