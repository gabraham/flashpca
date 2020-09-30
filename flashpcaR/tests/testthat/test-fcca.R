
context("Testing FCCA")

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

test_that("Testing fcca", {
   lambda1 <- 1e-3
   lambda2 <- 1e-2
   gamma1 <- 1e-3
   gamma2 <- 1e-2
   ndim <- 3
   s1 <- fcca(X, Y, lambda1=lambda1, lambda2=lambda2,
      gamma1=gamma1, gamma2=gamma2,
      ndim=ndim, standx="none", standy="none")
   expect_true(s1$converged)
   expect_true(is(s1, "fcca"))
   expect_equal(s1$ndim, ndim)
   expect_equal(length(s1$d), ndim)
   expect_equal(ncol(s1$U), ndim)
   expect_equal(ncol(s1$V), ndim)
   expect_equal(ncol(s1$a), ndim)
   expect_equal(ncol(s1$b), ndim)
   expect_equal(ncol(s1$Px), ndim)
   expect_equal(ncol(s1$Py), ndim)
   expect_equal(nrow(s1$Px), nrow(X))
   expect_equal(nrow(s1$Py), nrow(X))
   expect_equal(nrow(s1$U), ncol(X))
   expect_equal(nrow(s1$V), ncol(Y))
   expect_equal(nrow(s1$a), ncol(X))
   expect_equal(nrow(s1$b), ncol(Y))
   expect_true(all(s1$d > 0))

   expect_error(
      fcca(X, Y, lambda1=lambda1, lambda2=lambda2,
	 gamma1=gamma1, gamma2=gamma2,
	 ndim=0, standx="none", standy="none"))
   expect_error(
      fcca(X, NULL, lambda1=lambda1, lambda2=lambda2,
	 gamma1=gamma1, gamma2=gamma2,
	 ndim=ndim, standx="none", standy="none"))
   expect_error(
      fcca(NULL, Y, lambda1=lambda1, lambda2=lambda2,
	 gamma1=gamma1, gamma2=gamma2,
	 ndim=ndim, standx="none", standy="none"))
   expect_error(
      fcca(X, Y, lambda1=lambda1, lambda2=NA,
	 gamma1=gamma1, gamma2=gamma2,
	 ndim=ndim, standx="none", standy="none"))
   expect_error(
      fcca(X, Y, lambda1=-1, lambda2=lambda2,
	 gamma1=gamma1, gamma2=gamma2,
	 ndim=ndim, standx="none", standy="none"))
   expect_error(
      fcca(X, Y, lambda1=lambda1, lambda2=lambda2,
	 gamma1=gamma1, gamma2=gamma2,
	 ndim=ndim, standx="blah", standy="none"))
   expect_error(
      fcca(X, Y, lambda1=lambda1, lambda2=lambda2,
	 gamma1=gamma1, gamma2=gamma2,
	 ndim=ndim, standx="none", standy="foo"))
   expect_error(
      fcca(X, Y, lambda1=lambda1, lambda2=lambda2,
	 gamma1=-1, gamma2=gamma2,
	 ndim=ndim, standx="none", standy="none"))
   expect_error(
      fcca(X, Y, lambda1=lambda1, lambda2=lambda2,
	 gamma1=gamma1, gamma2=NULL,
	 ndim=ndim, standx="none", standy="none"))
})

test_that("Testing cv.fcca", {
   lambda1 <- seq(1e-6, 1e-2, length=3)
   lambda2 <- seq(1e-6, 1e-2, length=3)
   gamma1 <- 10^c(-2, -1)
   gamma2 <- 10^c(-3, -2)
   ndim <- 3
   nfolds <- 2
   s1 <- cv.fcca(X, Y, lambda1=lambda1, lambda2=lambda2,
      gamma1=gamma1, gamma2=gamma2, ndim=ndim, nfolds=nfolds)

   expect_equal(s1$ndim, ndim)
   expect_equal(s1$nfolds, nfolds)
   expect_equal(s1$lambda1, lambda1)
   expect_equal(s1$lambda2, lambda2)
   expect_equal(s1$gamma1, gamma1)
   expect_equal(s1$gamma2, gamma2)

   expect_error(
      cv.fcca(X, NULL, lambda1=lambda1, lambda2=lambda2,
	 gamma1=gamma1, gamma2=gamma2, ndim=ndim, nfolds=nfolds)
   )
   expect_error(
      cv.fcca(NULL, Y, lambda1=lambda1, lambda2=lambda2,
	 gamma1=gamma1, gamma2=gamma2, ndim=ndim, nfolds=nfolds)
   )
   expect_error(
      cv.fcca(X, Y, lambda1=NULL, lambda2=lambda2,
	 gamma1=gamma1, gamma2=gamma2, ndim=ndim, nfolds=nfolds)
   )
   expect_error(
      cv.fcca(X, Y, lambda1=lambda1, lambda2=NULL,
	 gamma1=gamma1, gamma2=gamma2, ndim=ndim, nfolds=nfolds)
   )
   expect_error(
      cv.fcca(X, Y, lambda1=lambda1, lambda2=lambda2,
	 gamma1=NULL, gamma2=gamma2, ndim=ndim, nfolds=nfolds)
   )
   expect_error(
      cv.fcca(X, Y, lambda1=lambda1, lambda2=lambda2,
	 gamma1=gamma1, gamma2=NULL, ndim=ndim, nfolds=nfolds)
   )
   expect_error(
      cv.fcca(X, Y, lambda1=lambda1, lambda2=lambda2,
	 gamma1=gamma1, gamma2=gamma2, ndim=-1, nfolds=nfolds)
   )
   expect_error(
      cv.fcca(X, Y, lambda1=lambda1, lambda2=lambda2,
	 gamma1=gamma1, gamma2=gamma2, ndim=ndim, nfolds=0)
   )

   # Compare the results from cross-validation with the single model.
   # The cross-validation is set up to use the same exact data
   # and parameters as the single model so should be identical (up to
   # reasonable numerical precision etc).
   X2 <- rbind(X, X)
   Y2 <- rbind(Y, Y)

   folds <- rep(1:2, each=nrow(X))

   s1 <- fcca(X, Y, lambda1=lambda1[1], lambda2=lambda2[1],
      gamma1=gamma1[1], gamma2=gamma2[1], ndim=ndim, maxiter=5000)
   s2 <- cv.fcca(X2, Y2, lambda1=lambda1[1], lambda2=lambda2[1],
      gamma1=gamma1[1], gamma2=gamma2[1], ndim=ndim, folds=folds,
      maxiter=5000, return.models=TRUE)

   # Extract the model using the saved indexing
   # The indexing is a bit clunky; the results have one row for each
   # dimension, but because each model contains all dimensions, the
   # model index is repeated several times (once for each dimension).
   #
   # fold 1
   idx1 <- as.numeric(
      s2$result.raw[1, .(idx.i, idx.j, idx.k, idx.m, idx.n)])
   s2.mod1 <- s2$models[[
      idx1[1] ]][[ idx1[2] ]][[ idx1[3] ]][[ idx1[4] ]][[ idx1[5] ]]$model
   # fold 2
   idx2 <- as.numeric(
      s2$result.raw[4, .(idx.i, idx.j, idx.k, idx.m, idx.n)])
   s2.mod2 <- s2$models[[
      idx2[1] ]][[ idx2[2] ]][[ idx2[3] ]][[ idx2[4] ]][[ idx2[5] ]]$model

   expect_equal(s1$ndim, s2.mod1$ndim)
   expect_true(mean(abs(s1$U - s2.mod1$U)) < 1e-4)
   expect_true(mean(abs(s1$V - s2.mod1$V)) < 1e-4)
   expect_equal(diag(cor(s1$U, s2.mod1$U)), rep(1, ndim), tol=1e-6)
   expect_equal(diag(cor(s1$V, s2.mod1$V)), rep(1, ndim), tol=1e-6)
   expect_equal(diag(cor(s1$a, s2.mod1$a)), rep(1, ndim), tol=1e-6)
   expect_equal(diag(cor(s1$b, s2.mod1$b)), rep(1, ndim), tol=1e-6)
   expect_equal(diag(cor(s1$Px, s2.mod1$Px)), rep(1, ndim), tol=1e-6)
   expect_equal(diag(cor(s1$Py, s2.mod1$Py)), rep(1, ndim), tol=1e-6)
   expect_equal(s1$d, s2.mod1$d, tol=1e-5)

   expect_equal(s1$ndim, s2.mod2$ndim)
   expect_true(mean(abs(s1$U - s2.mod2$U)) < 1e-4)
   expect_true(mean(abs(s1$V - s2.mod2$V)) < 1e-4)
   expect_equal(diag(cor(s1$U, s2.mod2$U)), rep(1, ndim), tol=1e-6)
   expect_equal(diag(cor(s1$V, s2.mod2$V)), rep(1, ndim), tol=1e-6)
   expect_equal(diag(cor(s1$a, s2.mod2$a)), rep(1, ndim), tol=1e-6)
   expect_equal(diag(cor(s1$b, s2.mod2$b)), rep(1, ndim), tol=1e-6)
   expect_equal(diag(cor(s1$Px, s2.mod2$Px)), rep(1, ndim), tol=1e-6)
   expect_equal(diag(cor(s1$Py, s2.mod2$Py)), rep(1, ndim), tol=1e-6)
   expect_equal(s1$d, s2.mod2$d, tol=1e-5)
})

