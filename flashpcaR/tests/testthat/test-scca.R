context("Testing SCCA")

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

test_that("Testing fcca", {
   #lambda1 <- seq(1e-6, 1e-2, length=10)
   #lambda2 <- seq(1e-6, 1e-2, length=10)
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


test_that("Testing optim.cv.fcca", {
   lambda1 <- seq(1e-5, 5e-2, length=4)
   lambda2 <- seq(1e-6, 1e-2, length=3)
   gamma1 <- 10^c(-2, -1)
   gamma2 <- 10^c(-3, 0)
   ndim <- 3
   nfolds <- 3
   folds <- sample(1:nfolds, nrow(X), replace=TRUE)

   # Test the grid optimisation, without returning models
   res1 <- optim.cv.fcca(X, Y, ndim=ndim, folds=folds,
      lambda1.grid=lambda1, lambda2.grid=lambda2,
      gamma1.grid=gamma1, gamma2.grid=gamma2,
      method="grid")
   
   expect_equal(nrow(res1$grid.path),
      length(lambda1) * length(lambda2) * length(gamma1) * length(gamma2))
   expect_equal(res1$nfolds, nfolds)
   expect_equal(res1$ndim, ndim)
   expect_null(res1$final.model.cv)
   expect_null(res1$final.model.cv.Px)
   expect_null(res1$final.model.cv.Py)
   expect_is(res1$final.model, "fcca")
   expect_equal(res1$final.model$ndim, ndim)
   expect_true(res1$opt.param["lambda1"] %in% lambda1)
   expect_true(res1$opt.param["lambda2"] %in% lambda2)
   expect_true(res1$opt.param["gamma1"] %in% gamma1)
   expect_true(res1$opt.param["gamma2"] %in% gamma2)

   # Test the grid optimisation, with returning models
   res2 <- optim.cv.fcca(X, Y, ndim=ndim, folds=folds,
      lambda1.grid=lambda1, lambda2.grid=lambda2,
      gamma1.grid=gamma1, gamma2.grid=gamma2,
      method="grid", final.model=TRUE, final.model.cv=TRUE)
   
   expect_equal(nrow(res2$grid.path),
      length(lambda1) * length(lambda2) * length(gamma1) * length(gamma2))
   expect_equal(res2$nfolds, nfolds)
   expect_equal(res2$ndim, ndim)
   expect_is(res2$final.model.cv, "cv.fcca")
   expect_equal(res2$final.model.cv$folds, folds)
   expect_equal(ncol(res2$final.model.cv.Px), ndim)
   expect_equal(ncol(res2$final.model.cv.Py), ndim)
   expect_equal(res2$final.model$ndim, ndim)
   expect_equal(res2$final.model.cv$lambda1, res2$opt.param["lambda1"])
   expect_equal(res2$final.model.cv$lambda2, res2$opt.param["lambda2"])
   expect_equal(res2$final.model.cv$gamma1, res2$opt.param["gamma1"])
   expect_equal(res2$final.model.cv$gamma2, res2$opt.param["gamma2"])
   
   # Test the Bayesian optimisation, if installed
   skip_if_not_installed("mlrMBO")
   skip_if_not_installed("DiceKriging")

   res3 <- optim.cv.fcca(X, Y, ndim=ndim, folds=folds,
      lambda1.grid=lambda1, lambda2.grid=lambda2,
      gamma1.grid=gamma1, gamma2.grid=gamma2,
      method="bopt", final.model=TRUE, final.model.cv=TRUE)
   
   expect_equal(nrow(res3$grid.path),
      length(lambda1) * length(lambda2) * length(gamma1) * length(gamma2))
   expect_equal(res3$nfolds, nfolds)
   expect_equal(res3$ndim, ndim)
   expect_is(res3$final.model.cv, "cv.fcca")
   expect_equal(res3$final.model.cv$folds, folds)
   expect_equal(ncol(res3$final.model.cv.Px), ndim)
   expect_equal(ncol(res3$final.model.cv.Py), ndim)
   expect_equal(res3$final.model$ndim, ndim)
   expect_equal(res3$final.model.cv$lambda1, res3$opt.param["lambda1"])
   expect_equal(res3$final.model.cv$lambda2, res3$opt.param["lambda2"])
   expect_equal(res3$final.model.cv$gamma1, res3$opt.param["gamma1"])
   expect_equal(res3$final.model.cv$gamma2, res3$opt.param["gamma2"])
})

test_that("Testing optim.cv.fcca, more", {
   lambda1 <- seq(1e-5, 5e-2, length=4)
   lambda2 <- seq(1e-6, 1e-2, length=3)
   gamma1 <- 10^c(-2, -1)
   gamma2 <- 10^c(-3, 0)
   ndim <- 3
   nfolds <- 3
   folds <- sample(1:nfolds, nrow(X), replace=TRUE)

   expect_error(optim.cv.fcca(X, Y, ndim=ndim, folds=folds,
      lambda1.grid=lambda1, lambda2.grid=NULL,
      gamma1.grid=gamma1, gamma2.grid=gamma2,
      method="grid"))

   expect_error(optim.cv.fcca(X, Y, ndim=ndim, folds=folds,
      lambda1.grid=NULL, lambda2.grid=lambda2.grid,
      gamma1.grid=gamma1, gamma2.grid=gamma2,
      method="grid"))

   expect_error(optim.cv.fcca(X, Y, ndim=ndim, folds=folds,
      lambda1.grid=c(-1, 2), lambda2.grid=lambda2.grid,
      gamma1.grid=gamma1, gamma2.grid=gamma2,
      method="grid"))

   expect_error(optim.cv.fcca(X, Y, ndim=ndim, folds=folds,
      lambda1.grid=lambda1, lambda2.grid=lambda2.grid,
      gamma1.grid=gamma1, gamma2.grid=gamma2,
      lambda1.bopt=-1, lambda2.bopt=NULL,
      method="bopt"))

   expect_error(optim.cv.fcca(X, Y, ndim=ndim, folds=folds,
      lambda1.grid=lambda1, lambda2.grid=lambda2.grid,
      gamma1.grid=gamma1, gamma2.grid=gamma2,
      lambda1.bopt=0, lambda2.bopt=0,
      gamma1.bopt=0, gamma2.bopt=0,
      method="bopt"))
})

