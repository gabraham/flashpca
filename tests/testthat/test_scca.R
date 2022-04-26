context("Testing SCCA")

## It's kind of hard to test that SCCA is working, but we can at least test
## that SCCA of X with X gives 
## - canonical correlations, i.e., diag(cor(Px, Py)), are equal to 1.
## - and that the canonical covariances ``d'' are the same as the eigenvalues
## of X^T X.
##
## We use very small penalisation so as to give very close results to
## classic eigen-decomposition.

data(hm3.chr1)

bedf <- gsub("\\.bed", "",
   system.file("extdata", "data_chr1.bed", package="flashpcaR"))

k <- 50
m <- ncol(hm3.chr1$bed)
X <- scale2(hm3.chr1$bed)
n <- nrow(X)
By <- matrix(rnorm(m * k), m, k)
Y <- scale(X %*% By + rnorm(n * k))
colnames(Y) <- paste0("marker", 1:k)

# very small penalties, to ensure the SCCA of X with X converges to the
# eigen-decomposition
l1 <- 1e-6
l2 <- 1e-6
test.tol <- 1e-4
ndim <- min(n, m, k, 5)

test_that(paste0("Testing self-self SCCA (X with X), l1=", l1, ", l2=", l2), {

   eval <- eigen(tcrossprod(X / sqrt(n - 1)))$val[1:ndim]
   
   # Essentially power method for eigen-decomposition of XX' / n
   s1 <- scca(X, X, lambda1=l1, lambda2=l2, ndim=ndim,	 
      standx="none", standy="none")
   s2 <- scca(bedf, X, lambda1=l1, lambda2=l2, ndim=ndim,	 
      standx="binom2", standy="none")
   r1 <- diag(cor(s1$Px, s1$Py))
   r2 <- diag(cor(s2$Px, s2$Py))

   expect_equal(s1$d, eval, tol=test.tol)
   expect_equal(s1$d, s2$d, tol=test.tol)

   expect_equal(rep(1, length(r1)), r1, tol=test.tol)
   expect_equal(rep(1, length(r2)), r2, tol=test.tol)
})

test_that(paste0("Testing self-self SCCA (X with X), l1=", l1, ", l2=", l2,
   "(no dividing by n)"), {

   eval <- eigen(tcrossprod(X))$val[1:ndim]
   
   # Essentially power method for eigen-decomposition of XX' / n
   s1 <- scca(X, X, lambda1=l1, lambda2=l2, ndim=ndim,	 
      standx="none", standy="none", divisor="none")
   s2 <- scca(bedf, X, lambda1=l1, lambda2=l2, ndim=ndim,	 
      standx="binom2", standy="none", divisor="none")
   r1 <- diag(cor(s1$Px, s1$Py))
   r2 <- diag(cor(s2$Px, s2$Py))

   expect_equal(s1$d, eval, tol=test.tol)
   expect_equal(s1$d, s2$d, tol=test.tol)

   expect_equal(rep(1, length(r1)), r1, tol=test.tol)
   expect_equal(rep(1, length(r2)), r2, tol=test.tol)
})

test_that("Testing self-self SCCA (X with X), initialising V0", {
   eval <- eigen(tcrossprod(X / sqrt(n - 1)))$val[1:ndim]

   Vx <- matrix(rnorm(m * ndim), m, ndim)
   
   # Essentially power method for eigen-decomposition of XX'/n
   s1 <- scca(X, X, lambda1=l1, lambda2=l2, ndim=ndim,	 
      standx="none", standy="none", V=Vx)
   s2 <- scca(bedf, X, lambda1=l1, lambda2=l2, ndim=ndim,	 
      standx="binom2", standy="none", V=Vx)
   r1 <- diag(cor(s1$Px, s1$Py))
   r2 <- diag(cor(s2$Px, s2$Py))

   expect_equal(s1$d, eval, tol=test.tol)
   expect_equal(s1$d, s2$d, tol=test.tol)

   expect_equal(rep(1, length(r1)), r1, tol=test.tol)
   expect_equal(rep(1, length(r2)), r2, tol=test.tol)
})

test_that("Testing SCCA (X with Y)", {

   # These penalties don't have to be tiny
   l1 <- runif(1, 1e-6, 1e-3)
   l2 <- runif(1, 1e-6, 1e-3)

   s1 <- scca(X, Y, lambda1=l1, lambda2=l2, ndim=ndim,	 
      standx="none", standy="none")
   s2 <- scca(bedf, Y, lambda1=l1, lambda2=l2, ndim=ndim,	 
      standx="binom2", standy="none")
   r1 <- diag(cor(s1$Px, s1$Py))
   r2 <- diag(cor(s2$Px, s2$Py))

   expect_equal(s1$d, s2$d, tol=test.tol)

   expect_equal(s1$V[,1], s2$V[,1], tol=test.tol)
   expect_equal(s1$V[,2], s2$V[,2], tol=test.tol)
   expect_equal(s1$U[,1], s2$U[,1], tol=test.tol)
   expect_equal(s1$U[,2], s2$U[,2], tol=test.tol)

   expect_equal(r1, r2, tol=test.tol)
})

test_that("Testing SCCA (X with Y, no divide by n)", {

   # These penalties don't have to be tiny
   l1 <- runif(1, 1e-6, 1e-3)
   l2 <- runif(1, 1e-6, 1e-3)

   s1 <- scca(X, Y, lambda1=l1, lambda2=l2, ndim=ndim,	 
      standx="none", standy="none")
   s2 <- scca(bedf, Y, lambda1=l1, lambda2=l2, ndim=ndim,	 
      standx="binom2", standy="none")

   s3 <- scca(X, Y, lambda1=l1, lambda2=l2, ndim=ndim,	 
      standx="none", standy="none", divisor="none")
   s4 <- scca(bedf, Y, lambda1=l1, lambda2=l2, ndim=ndim,	 
      standx="binom2", standy="none", divisor="none")

   r1 <- diag(cor(s1$Px, s1$Py))
   r2 <- diag(cor(s2$Px, s2$Py))

   r3 <- diag(cor(s3$Px, s3$Py))
   r4 <- diag(cor(s4$Px, s4$Py))

   n <- nrow(X)

   expect_equal(s1$d, s3$d / (n - 1), tol=test.tol)
   expect_equal(s2$d, s4$d / (n - 1), tol=test.tol)

   expect_equal(s1$V[,1], s3$V[,1], tol=test.tol)
   expect_equal(s2$V[,1], s4$V[,1], tol=test.tol)

   expect_equal(s1$V[,2], s3$V[,2], tol=test.tol)
   expect_equal(s2$V[,2], s4$V[,2], tol=test.tol)

   expect_equal(s1$U[,1], s3$U[,1], tol=test.tol)
   expect_equal(s2$U[,1], s4$U[,1], tol=test.tol)

   expect_equal(s1$U[,2], s3$U[,2], tol=test.tol)
   expect_equal(s2$U[,2], s4$U[,2], tol=test.tol)

   expect_equal(r1, r3, tol=test.tol)
   expect_equal(r2, r4, tol=test.tol)
})

test_that("Testing SCCA (X with Y), initialising V0", {

   # These penalties don't have to be tiny
   l1 <- runif(1, 1e-6, 1e-3)
   l2 <- runif(1, 1e-6, 1e-3)

   V0 <- matrix(rnorm(k * ndim), k, ndim)

   s1 <- scca(X, Y, lambda1=l1, lambda2=l2, ndim=ndim,	 
      standx="none", standy="none", V=V0)
   s2 <- scca(bedf, Y, lambda1=l1, lambda2=l2, ndim=ndim,	 
      standx="binom2", standy="none", V=V0)
   r1 <- diag(cor(s1$Px, s1$Py))
   r2 <- diag(cor(s2$Px, s2$Py))

   expect_equal(s1$d, s2$d, tol=test.tol)

   expect_equal(s1$V[,1], s2$V[,1], tol=test.tol)
   expect_equal(s1$V[,2], s2$V[,2], tol=test.tol)
   expect_equal(s1$U[,1], s2$U[,1], tol=test.tol)
   expect_equal(s1$U[,2], s2$U[,2], tol=test.tol)

   expect_equal(r1, r2, tol=test.tol)
})

test_that("Testing input checking", {
   
   # Test incompatible number of rows
   Z <- matrix(rnorm((nrow(X) + 3) * 100), nrow(X) + 3, 100)
   expect_error(scca(X, Z, lambda1=l1, lambda2=l2, ndim=ndim,	 
      standx="none", standy="none"))
   expect_error(scca(bedf, Z, lambda1=l1, lambda2=l2, ndim=ndim,	 
      standx="binom2", standy="none"))

   # Test negative penalties
   expect_error(scca(X, Z, lambda1=l1, lambda2=-1, ndim=ndim,	 
      standx="none", standy="none"))
   expect_error(scca(X, Z, lambda1=-1, lambda2=l2, ndim=ndim,	 
      standx="none", standy="none"))
   expect_error(scca(bedf, Z, lambda1=-1, lambda2=-1, ndim=ndim,	 
      standx="binom2", standy="none"))
})

test_that("Testing cv.scca", {
   lambda1 <- seq(1e-6, 1e-2, length=10)   
   lambda2 <- seq(1e-6, 1e-2, length=10)   
   s1 <- cv.scca(X, Y, lambda1=lambda1, lambda2=lambda2,
      ndim=ndim, standx="none", standy="none", nfolds=3)
   expect_true(all(s1$converged))
   expect_true(all(s1$corr >= 0.0))
   expect_true(all(s1$corr <= 1.0))
   expect_equal(s1$ndim, ndim)
   expect_true(s1$best.lambda1 > 0)
   expect_true(s1$best.lambda2 > 0)
   
   expect_error(
      cv.scca(bedf, Y, lambda1=lambda1, lambda2=lambda2,
       ndim=ndim, standx="none", standy="none", nfolds=3))
})

