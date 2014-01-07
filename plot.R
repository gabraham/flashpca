

d <- read.table("pcs.txt", header=FALSE, sep="")

colnames(d) <- paste("PC", 1:ncol(d))

#pdf("pcs.pdf")
png("pcs.png", width=1000, height=1000, res=100)
pairs(d[, 1:10])
dev.off()

