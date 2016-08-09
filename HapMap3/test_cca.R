rm(list=ls())

library(plink2R)

################################################################################
# Sparse CCA implementation in R

soft.thresh <- function(x, a)
{
   sign(x) * pmax(abs(x) - a, 0)
}

norm.thresh <- function(x, a)
{
   s <- sqrt(sum(x^2))
   if(s > 0) {
      x <- x / s
      x <- soft.thresh(x, a)
      s <- sqrt(sum(x^2))
      if(s > 0) {
         x <- x / s
      }
   }
   x
}

cca <- function(X, Y, lambdax=0, lambday=0, ndim=10, V=NULL,
   maxiter=100)
{
   k <- ncol(Y)
   p <- ncol(X)
   U <- matrix(0, p, ndim)
   d <- numeric(ndim)
   
   if(is.null(V)) {
      V <- matrix(rnorm(ncol(Y) * ndim), ncol(Y), ndim)
      cat("setting random V\n")
   }

   XY <- crossprod(X, Y)

   for(j in 1:ndim) {
      cat("dim ", j, "\n")
      if(j == 1) {
         XYj <- XY
      } else {
         XYj <- XYj - d[j - 1] * tcrossprod(U[, j - 1], V[, j - 1])
      }

      for(iter in 1:maxiter) {
	 cat(iter, " ")
         U[,j] <- XYj %*% V[, j]
         U[,j] <- norm.thresh(U[,j], lambdax)

         V[,j] <- crossprod(XYj, U[, j])
         V[,j] <- norm.thresh(V[,j], lambday)
      }
      cat("\n")
      #d[j] <- crossprod(X %*% U[,j], Y %*% V[,j])
      d[j] <- crossprod(U[,j], XYj) %*% V[,j]
   }

   list(u=U, d=d, v=V)
}


################################################################################
# Simulate data 

dat <- read_plink("data", impute="none")

scale2 <- function(X)
{
   p <- colSums(X, na.rm=TRUE) / (2 * colSums(!is.na(X)))
   S <- sweep(
      sweep(X, MARGIN=2, STATS=2 * p, FUN="-"),
      MARGIN=2, STATS=sqrt(2 * p * (1 - p)), FUN="/"
   )
   S[is.na(S)] <- 0
   attr(S, "scaled:center") <- 2 * p
   attr(S, "scaled:scale") <- sqrt(2 * p * (1 - p))
   S
}

X <- scale2(dat$bed)
n <- nrow(X)
p <- ncol(X)
k <- 25

B <- matrix(rnorm(p * k), p, k)
Y <- scale(X %*% B + rnorm(n * k)) 

d <- format(data.frame(dat$fam[, 1:2], Y), digits=6)
write.table(d, file="pheno.txt", col.names=FALSE, row.names=FALSE,
   quote=FALSE)

################################################################################
# Test univariate CCA (PLINK-multivariate). Equivalent to regression of SNP on
# all phenotypes.

system(paste0("../flashpca --bfile data --pheno pheno.txt --ucca"))

# sample a subset
w <- sample(ncol(X), 1000)
d.ucca <- read.table("ucca.txt", header=TRUE, sep="", stringsAsFactors=FALSE)
d.ucca <- d.ucca[w, ]

r <- lapply(w, function(j) {
   cat(j, " ")
   s <- summary(lm(X[,j] ~ Y))
   data.frame(SNP=dat$bim[j, 2], R=sqrt(s$r.squared),
      Fstat=s$fstatistic[1],
      P=pf(s$fstatistic["value"], s$fstatistic["numdf"],
	 s$fstatistic["dendf"], lower=FALSE))
})
d.lm <- do.call(rbind, r)

err.tol <- 1e-6
stopifnot(all(d.ucca$SNP == d.lm$SNP))
stopifnot(mean((d.ucca$R - d.lm$R)^2) < err.tol)
stopifnot(mean((d.ucca$Fstat - d.lm$Fstat)^2) < err.tol)
stopifnot(mean((d.ucca$P - d.lm$P)^2) < err.tol)

################################################################################
# Test sparse canonical correlation analysis (SCCA)

l1 <- 2e-2
l2 <- 2e-2

system(paste0(
   "../flashpca --bfile data --pheno pheno.txt --scca",
   " --lambda1 ", l1, " --lambda2 ", l2, " --ndim 10 --verbose",
   " --maxiter 100 --tol 1e-10 --precision 20"))
v0 <- matrix(scan("scca_v0.txt"), byrow=TRUE, nrow=k)
c1 <- cca(X, Y, lambdax=l1, lambday=l2, ndim=10, V=v0)

evecx <- matrix(scan("eigenvectorsX.txt"), byrow=TRUE, ncol=10)
evecy <- matrix(scan("eigenvectorsY.txt"), byrow=TRUE, ncol=10)
eval.obs <- scan("eigenvalues.txt")
px <- X %*% evecx
py <- Y %*% evecy
eval.obs2 <- diag(crossprod(px, py))

pcx <- matrix(scan("pcsX.txt"), byrow=TRUE, ncol=10)
pcy <- matrix(scan("pcsY.txt"), byrow=TRUE, ncol=10)
eval.obs3 <- diag(crossprod(pcx, pcy))

eval.c1 <- with(c1, diag(crossprod(X %*% u, Y %*% v)))

stopifnot(mean((eval.obs - eval.obs2)^2) < err.tol)
stopifnot(mean((c1$d - eval.obs)^2) < err.tol)
stopifnot(mean((c1$d - eval.obs2)^2) < err.tol)


