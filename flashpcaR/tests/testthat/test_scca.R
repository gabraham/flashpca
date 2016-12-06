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

# very small penalties, to ensure the SCCA of X with X converges to the
# eigen-decomposition
l1 <- 1e-6
l2 <- 1e-6
test.tol <- 1e-4
ndim <- min(n, m, k, 5)

test_that(paste0("Testing self-self SCCA (X with X), l1=", l1, ", l2=", l2), {

   eval <- eigen(tcrossprod(X))$val[1:ndim]
   
   # Essentially power method for eigen-decomposition of XX'
   s1 <- scca(X, X, lambda1=l1, lambda2=l2, ndim=ndim,	 
      standx="none", standy="none", mem="high")
   s2 <- scca(X, X, lambda1=l1, lambda2=l2, ndim=ndim,	 
      standx="none", standy="none", mem="low")
   s3 <- scca(bedf, X, lambda1=l1, lambda2=l2, ndim=ndim,	 
      standx="binom2", standy="none")
   r1 <- diag(cor(s1$Px, s1$Py))
   r2 <- diag(cor(s2$Px, s2$Py))
   r3 <- diag(cor(s3$Px, s3$Py))

   expect_equal(s1$d, eval, tol=test.tol)
   expect_equal(s1$d, s2$d, tol=test.tol)
   expect_equal(s1$d, s3$d, tol=test.tol)

   expect_equal(rep(1, length(r1)), r1, tol=test.tol)
   expect_equal(rep(1, length(r2)), r2, tol=test.tol)
})

test_that("Testing self-self SCCA (X with X), initialising V0", {
   eval <- eigen(tcrossprod(X))$val[1:ndim]

   Vx <- matrix(rnorm(m * ndim), m, ndim)
   
   # Essentially power method for eigen-decomposition of XX'
   s1 <- scca(X, X, lambda1=l1, lambda2=l2, ndim=ndim,	 
      standx="none", standy="none", mem="high", V=Vx)
   s2 <- scca(X, X, lambda1=l1, lambda2=l2, ndim=ndim,	 
      standx="none", standy="none", mem="low", V=Vx)
   s3 <- scca(bedf, X, lambda1=l1, lambda2=l2, ndim=ndim,	 
      standx="binom2", standy="none", V=Vx)
   r1 <- diag(cor(s1$Px, s1$Py))
   r2 <- diag(cor(s2$Px, s2$Py))
   r3 <- diag(cor(s3$Px, s3$Py))

   expect_equal(s1$d, eval, tol=test.tol)
   expect_equal(s1$d, s2$d, tol=test.tol)
   expect_equal(s1$d, s3$d, tol=test.tol)

   expect_equal(rep(1, length(r1)), r1, tol=test.tol)
   expect_equal(rep(1, length(r2)), r2, tol=test.tol)
})

test_that("Testing SCCA (X with Y)", {

   # These penalties don't have to be tiny
   l1 <- runif(1, 1e-6, 1e-3)
   l2 <- runif(1, 1e-6, 1e-3)

   # Essentially power method for eigen-decomposition of XX'
   s1 <- scca(X, Y, lambda1=l1, lambda2=l2, ndim=ndim,	 
      standx="none", standy="none", mem="high")
   s2 <- scca(X, Y, lambda1=l1, lambda2=l2, ndim=ndim,	 
      standx="none", standy="none", mem="low")
   s3 <- scca(bedf, Y, lambda1=l1, lambda2=l2, ndim=ndim,	 
      standx="binom2", standy="none")
   r1 <- diag(cor(s1$Px, s1$Py))
   r2 <- diag(cor(s2$Px, s2$Py))
   r3 <- diag(cor(s3$Px, s3$Py))

   expect_equal(s1$d, s2$d, tol=test.tol)
   expect_equal(s1$d, s3$d, tol=test.tol)

   expect_equal(s1$V[,1], s2$V[,1], tol=test.tol)
   expect_equal(s1$V[,1], s3$V[,1], tol=test.tol)

   expect_equal(s1$V[,2], s2$V[,2], tol=test.tol)
   expect_equal(s1$V[,2], s3$V[,2], tol=test.tol)

   expect_equal(s1$U[,1], s2$U[,1], tol=test.tol)
   expect_equal(s1$U[,1], s3$U[,1], tol=test.tol)

   expect_equal(s1$U[,2], s2$U[,2], tol=test.tol)
   expect_equal(s1$U[,2], s3$U[,2], tol=test.tol)

   expect_equal(r1, r2, tol=test.tol)
   expect_equal(r1, r3, tol=test.tol)
})

test_that("Testing SCCA (X with Y), initialising V0", {

   # These penalties don't have to be tiny
   l1 <- runif(1, 1e-6, 1e-3)
   l2 <- runif(1, 1e-6, 1e-3)

   V0 <- matrix(rnorm(k * ndim), k, ndim)

   # Essentially power method for eigen-decomposition of XX'
   s1 <- scca(X, Y, lambda1=l1, lambda2=l2, ndim=ndim,	 
      standx="none", standy="none", mem="high", V=V0)
   s2 <- scca(X, Y, lambda1=l1, lambda2=l2, ndim=ndim,	 
      standx="none", standy="none", mem="low", V=V0)
   s3 <- scca(bedf, Y, lambda1=l1, lambda2=l2, ndim=ndim,	 
      standx="binom2", standy="none", V=V0)
   r1 <- diag(cor(s1$Px, s1$Py))
   r2 <- diag(cor(s2$Px, s2$Py))
   r3 <- diag(cor(s3$Px, s3$Py))

   expect_equal(s1$d, s2$d, tol=test.tol)
   expect_equal(s1$d, s3$d, tol=test.tol)

   expect_equal(s1$V[,1], s2$V[,1], tol=test.tol)
   expect_equal(s1$V[,1], s3$V[,1], tol=test.tol)

   expect_equal(s1$V[,2], s2$V[,2], tol=test.tol)
   expect_equal(s1$V[,2], s3$V[,2], tol=test.tol)

   expect_equal(s1$U[,1], s2$U[,1], tol=test.tol)
   expect_equal(s1$U[,1], s3$U[,1], tol=test.tol)

   expect_equal(s1$U[,2], s2$U[,2], tol=test.tol)
   expect_equal(s1$U[,2], s3$U[,2], tol=test.tol)

   expect_equal(r1, r2, tol=test.tol)
   expect_equal(r1, r3, tol=test.tol)
})

test_that("Testing input checking", {
   
   # Test incompatible number of rows
   Z <- matrix(rnorm((nrow(X) + 3) * 100), nrow(X) + 3, 100)
   expect_error(scca(X, Z, lambda1=l1, lambda2=l2, ndim=ndim,	 
      standx="none", standy="none", mem="high"))
   expect_error(scca(X, Z, lambda1=l1, lambda2=l2, ndim=ndim,	 
      standx="none", standy="none", mem="low"))
   expect_error(scca(bedf, Z, lambda1=l1, lambda2=l2, ndim=ndim,	 
      standx="binom2", standy="none"))

   # Test negative penalties
   expect_error(scca(X, Z, lambda1=l1, lambda2=-1, ndim=ndim,	 
      standx="none", standy="none", mem="high"))
   expect_error(scca(X, Z, lambda1=-1, lambda2=l2, ndim=ndim,	 
      standx="none", standy="none", mem="low"))
   expect_error(scca(bedf, Z, lambda1=-1, lambda2=-1, ndim=ndim,	 
      standx="binom2", standy="none"))
})

