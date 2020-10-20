
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
   expect_equal(ncol(s1$A), ndim)
   expect_equal(ncol(s1$B), ndim)
   expect_equal(ncol(s1$Px), ndim)
   expect_equal(ncol(s1$Py), ndim)
   expect_equal(nrow(s1$Px), nrow(X))
   expect_equal(nrow(s1$Py), nrow(X))
   expect_equal(nrow(s1$U), ncol(X))
   expect_equal(nrow(s1$V), ncol(Y))
   expect_equal(nrow(s1$A), ncol(X))
   expect_equal(nrow(s1$B), ncol(Y))
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
      maxiter=5000, return_models=TRUE)

   # Extract the model using the saved indexing
   # The indexing is a bit clunky; the results have one row for each
   # dimension, but because each model contains all dimensions, the
   # model index is repeated several times (once for each dimension).
   #
   # fold 1
   idx1 <- as.numeric(
      s2$result_raw[1, c("idx.i", "idx.j", "idx.k", "idx.m", "idx.n")])
   s2.mod1 <- s2$models[[
      idx1[1] ]][[ idx1[2] ]][[ idx1[3] ]][[ idx1[4] ]][[ idx1[5] ]]$model
   # fold 2
   idx2 <- as.numeric(
      s2$result_raw[4, c("idx.i", "idx.j", "idx.k", "idx.m", "idx.n")])
   s2.mod2 <- s2$models[[
      idx2[1] ]][[ idx2[2] ]][[ idx2[3] ]][[ idx2[4] ]][[ idx2[5] ]]$model

   expect_is(s2.mod2, "scca")
   expect_equal(s1$ndim, s2.mod1$ndim)
   expect_true(mean(abs(s1$U - s2.mod1$U)) < 1e-4)
   expect_true(mean(abs(s1$V - s2.mod1$V)) < 1e-4)
   expect_equal(diag(cor(s1$U, s2.mod1$U)), rep(1, ndim), tol=1e-6)
   expect_equal(diag(cor(s1$V, s2.mod1$V)), rep(1, ndim), tol=1e-6)
   expect_equal(diag(cor(s1$A, s2.mod1$A)), rep(1, ndim), tol=1e-6)
   expect_equal(diag(cor(s1$B, s2.mod1$B)), rep(1, ndim), tol=1e-6)
   expect_equal(diag(cor(s1$Px, s2.mod1$Px)), rep(1, ndim), tol=1e-6)
   expect_equal(diag(cor(s1$Py, s2.mod1$Py)), rep(1, ndim), tol=1e-6)
   expect_equal(s1$d, s2.mod1$d, tol=1e-5)

   expect_equal(s1$ndim, s2.mod2$ndim)
   expect_true(mean(abs(s1$U - s2.mod2$U)) < 1e-4)
   expect_true(mean(abs(s1$V - s2.mod2$V)) < 1e-4)
   expect_equal(diag(cor(s1$U, s2.mod2$U)), rep(1, ndim), tol=1e-6)
   expect_equal(diag(cor(s1$V, s2.mod2$V)), rep(1, ndim), tol=1e-6)
   expect_equal(diag(cor(s1$A, s2.mod2$A)), rep(1, ndim), tol=1e-6)
   expect_equal(diag(cor(s1$B, s2.mod2$B)), rep(1, ndim), tol=1e-6)
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
      lambda1_grid=lambda1, lambda2_grid=lambda2,
      gamma1_grid=gamma1, gamma2_grid=gamma2,
      method="grid")
   
   expect_equal(nrow(res1$grid_path),
      length(lambda1) * length(lambda2) * length(gamma1) * length(gamma2))
   expect_equal(res1$nfolds, nfolds)
   expect_equal(res1$ndim, ndim)
   expect_null(res1$final_model_cv)
   expect_null(res1$final_model_cv_Px)
   expect_null(res1$final_model_cv_Py)
   expect_is(res1$final_model, "fcca")
   expect_equal(res1$final_model$ndim, ndim)
   expect_true(res1$opt_param["lambda1"] %in% lambda1)
   expect_true(res1$opt_param["lambda2"] %in% lambda2)
   expect_true(res1$opt_param["gamma1"] %in% gamma1)
   expect_true(res1$opt_param["gamma2"] %in% gamma2)

   # Test the grid optimisation, with returning models
   res2 <- optim.cv.fcca(X, Y, ndim=ndim, folds=folds,
      lambda1_grid=lambda1, lambda2_grid=lambda2,
      gamma1_grid=gamma1, gamma2_grid=gamma2,
      method="grid", final_model=TRUE, final_model_cv=TRUE)
   
   expect_equal(nrow(res2$grid_path),
      length(lambda1) * length(lambda2) * length(gamma1) * length(gamma2))
   expect_equal(res2$nfolds, nfolds)
   expect_equal(res2$ndim, ndim)
   expect_is(res2$final_model_cv, "cv.fcca")
   expect_equal(res2$final_model_cv$folds, folds)
   expect_equal(ncol(res2$final_model_cv_Px), ndim)
   expect_equal(ncol(res2$final_model_cv_Py), ndim)
   expect_equal(res2$final_model$ndim, ndim)
   expect_equal(res2$final_model_cv$lambda1, res2$opt_param["lambda1"])
   expect_equal(res2$final_model_cv$lambda2, res2$opt_param["lambda2"])
   expect_equal(res2$final_model_cv$gamma1, res2$opt_param["gamma1"])
   expect_equal(res2$final_model_cv$gamma2, res2$opt_param["gamma2"])
   
   # Test the Bayesian optimisation, if installed
   skip_if_not_installed("mlrMBO")
   skip_if_not_installed("DiceKriging")

   res3 <- optim.cv.fcca(X, Y, ndim=ndim, folds=folds,
      lambda1_grid=lambda1, lambda2_grid=lambda2,
      gamma1_grid=gamma1, gamma2_grid=gamma2,
      method="bopt", final_model=TRUE, final_model_cv=TRUE)
   
   expect_equal(nrow(res3$grid_path),
      length(lambda1) * length(lambda2) * length(gamma1) * length(gamma2))
   expect_equal(res3$nfolds, nfolds)
   expect_equal(res3$ndim, ndim)
   expect_is(res3$final_model_cv, "cv.fcca")
   expect_equal(res3$final_model_cv$folds, folds)
   expect_equal(ncol(res3$final_model_cv_Px), ndim)
   expect_equal(ncol(res3$final_model_cv_Py), ndim)
   expect_equal(res3$final_model$ndim, ndim)
   expect_equal(res3$final_model_cv$lambda1, res3$opt_param["lambda1"])
   expect_equal(res3$final_model_cv$lambda2, res3$opt_param["lambda2"])
   expect_equal(res3$final_model_cv$gamma1, res3$opt_param["gamma1"])
   expect_equal(res3$final_model_cv$gamma2, res3$opt_param["gamma2"])
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
      lambda1_grid=lambda1, lambda2_grid=NULL,
      gamma1_grid=gamma1, gamma2_grid=gamma2,
      method="grid"))

   expect_error(optim.cv.fcca(X, Y, ndim=ndim, folds=folds,
      lambda1_grid=NULL, lambda2_grid=lambda2_grid,
      gamma1_grid=gamma1, gamma2_grid=gamma2,
      method="grid"))

   expect_error(optim.cv.fcca(X, Y, ndim=ndim, folds=folds,
      lambda1_grid=c(-1, 2), lambda2_grid=lambda2_grid,
      gamma1_grid=gamma1, gamma2_grid=gamma2,
      method="grid"))

   expect_error(optim.cv.fcca(X, Y, ndim=ndim, folds=folds,
      lambda1_grid=lambda1, lambda2_grid=lambda2_grid,
      gamma1_grid=gamma1, gamma2_grid=gamma2,
      lambda1_bopt=-1, lambda2_bopt=NULL,
      method="bopt"))

   expect_error(optim.cv.fcca(X, Y, ndim=ndim, folds=folds,
      lambda1_grid=lambda1, lambda2_grid=lambda2_grid,
      gamma1_grid=gamma1, gamma2_grid=gamma2,
      lambda1_bopt=0, lambda2_bopt=0,
      gamma1_bopt=0, gamma2_bopt=0,
      method="bopt"))
})
