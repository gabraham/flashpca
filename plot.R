

d <- read.table("pcs.txt", header=FALSE, sep="")

colnames(d) <- paste("PC", 1:ncol(d))

pdf("pcs.pdf")
pairs(d)
dev.off()

