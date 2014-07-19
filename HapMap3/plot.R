
library(plink2R)
library(grid)

readmat <- function(f)
{
   con <- file(f, "rb")
   p <- readBin(con=con, what="integer", n=1)
   K <- readBin(con=con, what="integer", n=1)
   matrix(readBin(con=con, what="double", n=p * K), nrow=p, ncol=K)
}

dat <- read_plink("data", impute="random")

# Price 2006 standardisation
X <- apply(dat$bed, 2, function(x) {
   u <- mean(x)
   p <- u / 2
   (x - u) / sqrt(p * (1 - p))
})

pu <- colMeans(dat$bed) / 2

# reverse standardisation
R <- sapply(1:ncol(X), function(j) {
   x <- X[,j] - 2 * pu[j]
   x / sqrt(2 * pu[j] * (1 - pu[j]))
})

pr <- prcomp(X, center=FALSE)

max(abs(pr$x))

lf <- list.files(pattern="\\.pca\\.evec$")

fam <- read.table("data.fam", header=FALSE, sep="", stringsAsFactors=FALSE)
m <- read.table("relationships_w_pops_121708.txt", header=TRUE, sep="",
   stringsAsFactors=FALSE)
rownames(m) <- m[,2]
pop <- factor(m[fam[,2], 7])

################################################################################
# Compare the eigenvectors / principal components

# smartpca
d1 <- read.table("data.pca.evec",
   header=FALSE, sep="", stringsAsFactors=FALSE, skip=1, row.names=1)

# flashpca
d2 <- read.table("pcs.txt", header=FALSE, sep="")

# shellfish
d3 <- read.table("shellfish.evecs", header=FALSE, sep="")
   
k <- 10
x1 <- as.matrix(d1[, 1:k])
x2 <- as.matrix(d2[, 1:k])
x3 <- t(as.matrix(d3))[, 1:k]
x4 <- pr$x[, 1:k]

max(abs(x1))
max(abs(x2))
max(abs(x3))
max(abs(x4))

colnames(x1) <- paste("PC", 1:ncol(x1), sep="")
colnames(x2) <- paste("PC", 1:ncol(x2), sep="")
colnames(x3) <- paste("PC", 1:ncol(x3), sep="")
colnames(x4) <- paste("PC", 1:ncol(x4), sep="")
   
panel.cor <- function(x, y, ...)
{
   usr <- par("usr"); on.exit(par(usr))
   par(usr = c(0, 1, 0, 1))
   r <- cor(x, y)
   txt <- format(c(r, 0.123456789), digits=2)[1]
   #cex.cor <- 0.5 / strwidth(txt)
   cex.cor <- 1.5
   text(0.5, 0.5, txt, cex = cex.cor)
}

pdf("hapmap3_comparison.pdf")
for(i in 1:k) {
   z <- cbind(smartpca=x1[,i], flashpca=x2[,i],
      shellfish=x3[,i], R=x4[,i])
   pairs(z, main=paste("PC", i), col=pop, lower.panel=panel.cor)
}
dev.off()

# Figure 1a
pdf("hapmap3.pdf", width=5.5, height=5.5)
plot(x2[, 1:2], col=pop, pch=as.integer(pop))
legend(x=min(x2[,1]), y=1.6,
   legend=levels(pop), col=1:length(levels(pop)),
   lwd=3, lty=1, pch=1:length(levels(pop)), ncol=2,
   box.lty=0)
grid.text("a", x=unit(0.04, "npc"), y=unit(0.96, "npc"),
   gp=gpar(fontsize=40))
dev.off()

# Figure 1b
i <- 1
z <- cbind(smartpca=x1[,i], flashpca=x2[,i],
   shellfish=x3[,i], R=x4[,i])
pdf("hapmap3_pairs.pdf", width=5.5, height=5.5)
pairs(z, main=paste("PC", i), col=pop, lower.panel=panel.cor,
   cex=1.4, cex.axis=1.4)
grid.text("b", x=unit(0.04, "npc"), y=unit(0.96, "npc"),
   gp=gpar(fontsize=40))
dev.off()

d2 <- read.table("pcs.txt", header=FALSE, sep="")
x2 <- as.matrix(d2[, 1:k])
i <- 1
z <- cbind(smartpca=x1[,i], flashpca=x2[,i],
     shellfish=x3[,i], R=x4[,i])
cor(z)

################################################################################
# Compare the eigenvalues (squared principal values)

v1 <- scan("data.eval")[1:k]^2
v2 <- scan("eigenvalues.txt")[1:k]
v3 <- scan("shellfish.evals")[1:k]^2
v4 <- pr$sdev[1:k]^2

eig <- cbind(smartpca=v1, flashpca=v2, shellfish=v3, R=v4)
cor(eig)



