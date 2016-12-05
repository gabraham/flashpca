context("Testing UCCA")

n <- 500
p <- 100
k <- 15
m <- 10
M <- matrix(rnorm(n * m), n, m)
Bx <- matrix(rnorm(m * p), m, p)
By <- matrix(rnorm(m * k), m, k)
X0 <- scale(M %*% Bx + rnorm(n * p))
Y0 <- scale(M %*% By + rnorm(n * k))

data(hm3.chr1)

bedf <- gsub("\\.bed", "",
   system.file("extdata", "data_chr1.bed", package="flashpcaR"))

test.tol <- 1e-4

test <- function(s, X, Y)
{
   p <- ncol(X)

   r2.obs <- s$result[,"R"]^2
   r2.exp <- numeric(p)

   F.obs <- s$result[, "Fstat"]
   F.exp <- numeric(p)
   
   p.obs <- s$result[, "P"]
   p.exp <- numeric(p)

   # expect_equal also tests names
   names(r2.exp) <- names(F.exp) <- names(p.exp) <- colnames(X)

   for(i in 1:p) {
      l <- lm(X[,i] ~ Y)
      l2 <- summary(l)
      r2.exp[i] <- l2$r.squared
      l3 <- anova(l)
      F.exp[i] <- l3[1, "F value"]
      p.exp[i] <- l3[1, "Pr(>F)"]
   }

   expect_equal(r2.exp, r2.obs, tol=test.tol)
   expect_equal(F.exp, F.obs, tol=test.tol)
   expect_equal(p.exp, p.obs, tol=test.tol)
   expect_equal(log(p.exp), log(p.obs), tol=test.tol)
}

test_that("Testing UCCA (binomial) with matrices and PLINK", {
   X <- scale2(hm3.chr1$bed, type="1")
	
   n <- nrow(X)
   p <- ncol(X)
   By <- matrix(rnorm(p * k), p, k)
   Y <- scale(X %*% By + rnorm(n * k))

   s1 <- ucca(X, Y, standx="none", standy="none")
   s2 <- ucca(bedf, Y, standx="binom", standy="none")
   test(s1, X, Y)
   test(s2, X, Y)
})

test_that("Testing UCCA (binomial2) with matrices and PLINK", {
   X <- scale2(hm3.chr1$bed, type="2")
	
   n <- nrow(X)
   p <- ncol(X)
   By <- matrix(rnorm(p * k), p, k)
   Y <- scale(X %*% By + rnorm(n * k))

   s1 <- ucca(X, Y, standx="none", standy="none")
   s2 <- ucca(bedf, Y, standx="binom2", standy="none")
   test(s1, X, Y)
   test(s2, X, Y)
})

#test_that("Testing UCCA (binomial) with PLINK data", {
#
#   s1 <- ucca(bedf, Y0, standx="binom", standy="none")
#   test(s1, Xb, Y0)
#})

#test_that("Testing self-self UCCA with stand='sd'", {
#   # X0 is already unit-variance so scaling within scca is moot but still
#   # useful as a sanity check
#   s1 <- scca(X0, X0, lambda1=l1, lambda2=l2, ndim=ndim,
#      standx="sd", standy="sd", mem="high")
#   s2 <- scca(X0, X0, lambda1=l1, lambda2=l2, ndim=ndim,
#      standx="sd", standy="sd", mem="low")
#   eval <- eigen(crossprod(X0))$val[1:ndim]
#
#   r1 <- diag(cor(s1$Px, s1$Py))
#   r2 <- diag(cor(s2$Px, s2$Py))
#
#   expect_equal(rep(1, length(r1)), r1, tol=test.tol)
#   expect_equal(rep(1, length(r2)), r2, tol=test.tol)
#
#   expect_equal(eval, s1$d, tol=test.tol)
#   expect_equal(eval, s2$d, tol=test.tol)
#})
#
#test_that("Testing self-self UCCA with stand='none'", {
#   # X0 is already unit-variance so not scaling within scca is ok
#   s1 <- scca(X0, X0, lambda1=l1, lambda2=l2, ndim=ndim,
#      standx="none", standy="none", mem="high")
#   s2 <- scca(X0, X0, lambda1=l1, lambda2=l2, ndim=ndim,
#      standx="none", standy="none", mem="low")
#
#   eval <- eigen(crossprod(X0))$val[1:ndim]
#
#   r1 <- diag(cor(s1$Px, s1$Py))
#   r2 <- diag(cor(s2$Px, s2$Py))
#
#   expect_equal(rep(1, length(r1)), r1, tol=test.tol)
#   expect_equal(rep(1, length(r2)), r2, tol=test.tol)
#
#   expect_equal(eval, s1$d, tol=test.tol)
#   expect_equal(eval, s2$d, tol=test.tol)
#})
#
#test_that("Testing self-self UCCA with stand='center'", {
#   # X0 is already zero-mean so scaling within scca is moot but still
#   # useful as a sanity check
#   s1 <- scca(X0, X0, lambda1=l1, lambda2=l2, ndim=ndim,
#      standx="center", standy="center", mem="high")
#   s2 <- scca(X0, X0, lambda1=l1, lambda2=l2, ndim=ndim,
#      standx="center", standy="center", mem="low")
#   eval <- eigen(crossprod(X0))$val[1:ndim]
#
#   r1 <- diag(cor(s1$Px, s1$Py))
#   r2 <- diag(cor(s2$Px, s2$Py))
#
#   expect_equal(rep(1, length(r1)), r1, tol=test.tol)
#   expect_equal(rep(1, length(r2)), r2, tol=test.tol)
#
#   expect_equal(eval, s1$d, tol=test.tol)
#   expect_equal(eval, s2$d, tol=test.tol)
#})
#
#test_that("Testing self-self UCCA (X with X), initialising V0", {
#
#   # Round the continuous values into {0, 1, 2}
#   Xb <- apply(X0, 2, function(x) {
#      round(2 * (x - min(x)) / (max(x) - min(x)))
#   })
#
#   # "binomial" standardization
#   q <- colMeans(Xb) / 2
#   Xbs <- scale(Xb, center=TRUE, scale=sqrt(q * (1 - q)))
#   sv <- svd(Xbs)
#   eval <- sv$d[1:ndim]^2
#   
#   # Essentially power method for eigen-decomposition of XX'
#   s1 <- scca(Xb, Xb, lambda1=l1, lambda2=l2, ndim=ndim,	 
#      standx="binom", standy="binom", V=sv$v[, 1:ndim], mem="high")
#   s2 <- scca(Xb, Xb, lambda1=l1, lambda2=l2, ndim=ndim,	 
#      standx="binom", standy="binom", V=sv$v[, 1:ndim], mem="low")
#
#   r1 <- diag(cor(s1$Px, s1$Py))
#   r2 <- diag(cor(s2$Px, s2$Py))
#
#   expect_equal(rep(1, length(r1)), r1, tol=test.tol)
#   expect_equal(rep(1, length(r2)), r2, tol=test.tol)
#
#   expect_equal(eval, s1$d, tol=test.tol)
#   expect_equal(eval, s2$d, tol=test.tol)
#})
#
