
library(doMC)
library(fields)
library(ggplot2)

registerDoMC(cores=20)

# Your FAM file here
famfile <- ""

# Your phenotype file here
phenofile <- ""

# Number of eigenvectors used
ndim <- 10

lpx <- list.files(pattern="pcsX")
lpy <- list.files(pattern="pcsY")
lpu <- list.files(pattern="eigenvectorsX")
names(lpx) <- lpx
names(lpy) <- lpy
names(lpu) <- lpu
Px <- lapply(lpx, function(x) matrix(scan(x), byrow=TRUE, ncol=ndim))
Py <- lapply(lpy, function(x) matrix(scan(x), byrow=TRUE, ncol=ndim))
nzu <- foreach(x=lpu) %dopar% {
   x <- matrix(scan(x), byrow=TRUE, ncol=ndim)
   colSums(x != 0)
}
R.trn <- sapply(seq(along=Px), function(i) diag(cor(Px[[i]], Py[[i]])))
lambda <- sapply(strsplit(lpx, "_|\\.txt"), function(x) as.numeric(x[2:3]))
R.trn1 <- matrix(R.trn[1,], byrow=TRUE, ncol=3)
nzu1 <- matrix(sapply(nzu, function(x) x[1]), byrow=TRUE, ncol=3)[,1]
lambda1 <- matrix(lambda[1,], byrow=TRUE, ncol=3)[, 1]

png("R_train.png", width=2000, height=1000, res=200)
par(mfrow=c(1, 2), mar=c(4, 4, 3, 1))
matplot(lambda1, R.trn1, log="x")
matplot(nzu1, R.trn1, log="x")
dev.off()

# Your phenotype file here
d <- read.table(phenofile, header=FALSE, sep="", stringsAsFactors=FALSE)

# Your FAM file here
fam <- read.table(famfile, header=FALSE, sep="", stringsAsFactors=FALSE)

rownames(d) <- d[,1]
dtest <- d[fam[,1],]

Y <- as.matrix(dtest[, -(1:2)])

# TODO: scale by original scaling of Y in training data
Y <- scale(Y)

lfy <- list.files(pattern="eigenvectorsY_")

V <- lapply(lfy, function(f) {
   matrix(scan(f), byrow=TRUE, ncol=ndim)
})

Py <- lapply(V, function(v) Y %*% v)

lfx <- list.files(pattern="predX")
Px <- lapply(lfx, function(f) {
   matrix(scan(f), byrow=TRUE, ncol=ndim)
})

r <- sapply(strsplit(lfy, "_|\\.txt"), function(x) x[2:3])
l1 <- as.numeric(r[1,])
l2 <- as.numeric(r[2,])

R <- lapply(seq(along=Px), function(i) {
   cor(Px[[i]], Py[[i]])
})

Rmax <- sapply(R, max)

res <- data.frame(l1=l1, l2=factor(l2), R=Rmax)

g <- ggplot(res, aes(x=l1, y=R, colour=l2))
g <- g + geom_line()
g <- g + geom_point()
print(g)
dev.off()

