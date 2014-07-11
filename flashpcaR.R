
library(plink2R)
library(flashpcaR)

set.seed(1)

n <- 500
p <- 1000

#X <- matrix(sample(0:2, n * p, replace=TRUE), n, p)
dat <- read_plink("small", impute="random")
X <- dat$bed

# Price 2006 standardisation
#S <- apply(X, 2, function(x) {
#   u <- mean(x)
#   p <- u / 2
#   (x - u) / sqrt(p * (1 - p))
#})
S <- scale(X, center=TRUE, scale=TRUE)
stand <- "sd"

s <- svd(S / sqrt(nrow(S) - 1))
K <- 20

r1 <- flashpca(X, do_loadings=TRUE, verbose=TRUE, stand=stand, maxiter=100,
   tol=1e-6, transpose=FALSE, ndim=K)
r2 <- flashpca(t(X), do_loadings=TRUE, verbose=TRUE, stand=stand,
   maxiter=100, tol=1e-6, transpose=TRUE, ndim=K)

D <- cbind(r1=r1$values, r2=r2$values, eig=s$d[1:K]^2)
cor(D)

U1 <- cbind(r1$vectors[,1], r2$vectors[,1], s$u[,1])
cor(U1)

V1 <- cbind(r1$loadings[,1], r2$loadings[,1], s$v[,1])
cor(V1)

diag(crossprod(r1$vectors))
diag(crossprod(r2$vectors))
diag(crossprod(r1$loadings))
diag(crossprod(r2$loadings))

