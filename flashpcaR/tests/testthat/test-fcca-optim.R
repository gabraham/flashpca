
context("Testing FCCA optim")

# otherwise .() won't be found during testing
library(data.table)

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

test_that("Testing optim.cv.fcca", {

   skip_on_cran()
   
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

   # Sanity check of test-set predictions
   res.fcca.cv <- list(
      Px=matrix(0, nrow(X), ndim, dimnames=list(NULL, paste0("Px", 1:ndim))),
      Py=matrix(0, nrow(X), ndim, dimnames=list(NULL, paste0("Py", 1:ndim)))
   )

   Xs <- scale(X)
   Ys <- scale(Y)
   for(fold in 1:nfolds)
   {
      m1 <- fcca(Xs[folds != fold,], Ys[folds != fold,],
         ndim=ndim, standx="none", standy="none",
         lambda1=res3$opt_param["lambda1"],
         lambda2=res3$opt_param["lambda2"],
         gamma1=res3$opt_param["gamma1"],
         gamma2=res3$opt_param["gamma2"]
      )
      res.fcca.cv$Px[folds == fold, ] <- Xs[folds == fold,] %*% m1$A
      res.fcca.cv$Py[folds == fold, ] <- Ys[folds == fold,] %*% m1$B
   }

   expect_equivalent(
      diag(cor(res.fcca.cv$Px, res3$final_model_cv_Px)),
      rep(1, ndim), tolerance=1e-3)
   expect_equivalent(
      diag(cor(res.fcca.cv$Py, res3$final_model_cv_Py)),
      rep(1, ndim), tolerance=1e-3)
})

test_that("Testing optim.cv.fcca, more", {

   skip_on_cran()
   
   lambda1_grid <- 1e-3
   lambda2_grid <- seq(1e-6, 1e-2, length=3)
   gamma1_grid <- 10^c(-2, -1)
   gamma2_grid <- 10^c(-3, 0)

   lambda1_bopt <- 1e-3
   lambda2_bopt <- c(0, 1e-1)
   gamma1_bopt <- 10^c(-3, 3)
   gamma2_bopt <- 10^c(-3, 3)

   ndim <- 3
   nfolds <- 3
   folds <- sample(1:nfolds, nrow(X), replace=TRUE)

   # Test the grid optimisation, without returning models
   res1 <- optim.cv.fcca(X, Y, ndim=ndim, folds=folds,
      lambda1_grid=lambda1_grid, lambda2_grid=lambda2_grid,
      gamma1_grid=gamma1_grid, gamma2_grid=gamma2_grid,
      method="grid")
   
   expect_equal(nrow(res1$grid_path),
      length(lambda1_grid) * length(lambda2_grid) 
      * length(gamma1_grid) * length(gamma2_grid))
   expect_equal(res1$nfolds, nfolds)
   expect_equal(res1$ndim, ndim)
   expect_null(res1$final_model_cv)
   expect_null(res1$final_model_cv_Px)
   expect_null(res1$final_model_cv_Py)
   expect_is(res1$final_model, "fcca")
   expect_equal(res1$final_model$ndim, ndim)
   expect_true(res1$opt_param["lambda1"] %in% lambda1_grid)
   expect_true(res1$opt_param["lambda2"] %in% lambda2_grid)
   expect_true(res1$opt_param["gamma1"] %in% gamma1_grid)
   expect_true(res1$opt_param["gamma2"] %in% gamma2_grid)

   # Test the grid optimisation, with returning models
   res2 <- optim.cv.fcca(X, Y, ndim=ndim, folds=folds,
      lambda1_grid=lambda1_grid, lambda2_grid=lambda2_grid,
      gamma1_grid=gamma1_grid, gamma2_grid=gamma2_grid,
      method="grid", final_model=TRUE, final_model_cv=TRUE)
   
   expect_equal(nrow(res2$grid_path),
      length(lambda1_grid) * length(lambda2_grid) 
      * length(gamma1_grid) * length(gamma2_grid))
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
      lambda1_grid=lambda1_grid, lambda2_grid=lambda2_grid,
      gamma1_grid=gamma1_grid, gamma2_grid=gamma2_grid,
      lambda1_bopt=lambda1_bopt, lambda2_bopt=lambda2_bopt,
      gamma1_bopt=gamma1_bopt, gamma2_bopt=gamma2_bopt,
      method="bopt", final_model=TRUE, final_model_cv=TRUE)
   
   expect_equal(nrow(res3$grid_path),
      length(lambda1_grid) * length(lambda2_grid) 
      * length(gamma1_grid) * length(gamma2_grid))
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

   # Sanity check of test-set predictions
   res.fcca.cv <- list(
      Px=matrix(0, nrow(X), ndim, dimnames=list(NULL, paste0("Px", 1:ndim))),
      Py=matrix(0, nrow(X), ndim, dimnames=list(NULL, paste0("Py", 1:ndim)))
   )

   Xs <- scale(X)
   Ys <- scale(Y)
   for(fold in 1:nfolds)
   {
      m1 <- fcca(Xs[folds != fold,], Ys[folds != fold,],
         ndim=ndim, standx="none", standy="none",
         lambda1=res3$opt_param["lambda1"],
         lambda2=res3$opt_param["lambda2"],
         gamma1=res3$opt_param["gamma1"],
         gamma2=res3$opt_param["gamma2"]
      )
      res.fcca.cv$Px[folds == fold, ] <- Xs[folds == fold,] %*% m1$A
      res.fcca.cv$Py[folds == fold, ] <- Ys[folds == fold,] %*% m1$B
   }

   expect_equivalent(
      diag(cor(res.fcca.cv$Px, res3$final_model_cv_Px)),
      rep(1, ndim), tolerance=1e-3)
   expect_equivalent(
      diag(cor(res.fcca.cv$Py, res3$final_model_cv_Py)),
      rep(1, ndim), tolerance=1e-3)

   res4 <- optim.cv.fcca(X, Y, ndim=ndim, folds=folds,
      lambda1_grid=lambda1_grid, lambda2_grid=lambda2_grid,
      gamma1_grid=0, gamma2_grid=0,
      lambda1_bopt=lambda1_bopt, lambda2_bopt=lambda2_bopt,
      gamma1_bopt=0, gamma2_bopt=0,
      method="bopt", final_model=TRUE, final_model_cv=TRUE)

   expect_equivalent(res4$opt_param["gamma1"], 0)
   expect_equivalent(res4$opt_param["gamma2"], 0)

   res5 <- optim.cv.fcca(X, Y, ndim=ndim, folds=folds,
      lambda1_grid=0, lambda2_grid=0,
      gamma1_grid=gamma1_grid, gamma2_grid=gamma2_grid,
      lambda1_bopt=0, lambda2_bopt=0,
      gamma1_bopt=gamma1_bopt, gamma2_bopt=gamma2_bopt,
      method="bopt", final_model=TRUE, final_model_cv=TRUE)

   expect_equivalent(res5$opt_param["lambda1"], 0)
   expect_equivalent(res5$opt_param["lambda2"], 0)

})

test_that("Testing optim.cv.fcca, argument checks", {
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

