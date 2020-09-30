
context("Testing FCCA optim")

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
         lambda1=res3$opt.param["lambda1"],
         lambda2=res3$opt.param["lambda2"],
         gamma1=res3$opt.param["gamma1"],
         gamma2=res3$opt.param["gamma2"]
      )
      res.fcca.cv$Px[folds == fold, ] <- Xs[folds == fold,] %*% m1$a
      res.fcca.cv$Py[folds == fold, ] <- Ys[folds == fold,] %*% m1$b
   }

   expect_equal(diag(cor(res.fcca.cv$Px, res3$final.model.cv.Px)),
      rep(1, ndim), tolerance=1e-3)
   expect_equal(diag(cor(res.fcca.cv$Py, res3$final.model.cv.Py)),
      rep(1, ndim), tolerance=1e-3)

})

test_that("Testing optim.cv.fcca, more", {

   skip_on_cran()
   
   lambda1.grid <- 1e-3
   lambda2.grid <- seq(1e-6, 1e-2, length=3)
   gamma1.grid <- 10^c(-2, -1)
   gamma2.grid <- 10^c(-3, 0)

   lambda1.bopt <- 1e-3
   lambda2.bopt <- c(0, 1e-1)
   gamma1.bopt <- 10^c(-3, 3)
   gamma2.bopt <- 10^c(-3, 3)

   ndim <- 3
   nfolds <- 3
   folds <- sample(1:nfolds, nrow(X), replace=TRUE)

   # Test the grid optimisation, without returning models
   res1 <- optim.cv.fcca(X, Y, ndim=ndim, folds=folds,
      lambda1.grid=lambda1.grid, lambda2.grid=lambda2.grid,
      gamma1.grid=gamma1.grid, gamma2.grid=gamma2.grid,
      method="grid")
   
   expect_equal(nrow(res1$grid.path),
      length(lambda1.grid) * length(lambda2.grid) 
      * length(gamma1.grid) * length(gamma2.grid))
   expect_equal(res1$nfolds, nfolds)
   expect_equal(res1$ndim, ndim)
   expect_null(res1$final.model.cv)
   expect_null(res1$final.model.cv.Px)
   expect_null(res1$final.model.cv.Py)
   expect_is(res1$final.model, "fcca")
   expect_equal(res1$final.model$ndim, ndim)
   expect_true(res1$opt.param["lambda1"] %in% lambda1.grid)
   expect_true(res1$opt.param["lambda2"] %in% lambda2.grid)
   expect_true(res1$opt.param["gamma1"] %in% gamma1.grid)
   expect_true(res1$opt.param["gamma2"] %in% gamma2.grid)

   # Test the grid optimisation, with returning models
   res2 <- optim.cv.fcca(X, Y, ndim=ndim, folds=folds,
      lambda1.grid=lambda1.grid, lambda2.grid=lambda2.grid,
      gamma1.grid=gamma1.grid, gamma2.grid=gamma2.grid,
      method="grid", final.model=TRUE, final.model.cv=TRUE)
   
   expect_equal(nrow(res2$grid.path),
      length(lambda1.grid) * length(lambda2.grid) 
      * length(gamma1.grid) * length(gamma2.grid))
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
      lambda1.grid=lambda1.grid, lambda2.grid=lambda2.grid,
      gamma1.grid=gamma1.grid, gamma2.grid=gamma2.grid,
      lambda1.bopt=lambda1.bopt, lambda2.bopt=lambda2.bopt,
      gamma1.bopt=gamma1.bopt, gamma2.bopt=gamma2.bopt,
      method="bopt", final.model=TRUE, final.model.cv=TRUE)
   
   expect_equal(nrow(res3$grid.path),
      length(lambda1.grid) * length(lambda2.grid) 
      * length(gamma1.grid) * length(gamma2.grid))
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
         lambda1=res3$opt.param["lambda1"],
         lambda2=res3$opt.param["lambda2"],
         gamma1=res3$opt.param["gamma1"],
         gamma2=res3$opt.param["gamma2"]
      )
      res.fcca.cv$Px[folds == fold, ] <- Xs[folds == fold,] %*% m1$a
      res.fcca.cv$Py[folds == fold, ] <- Ys[folds == fold,] %*% m1$b
   }

   expect_equal(diag(cor(res.fcca.cv$Px, res3$final.model.cv.Px)),
      rep(1, ndim), tolerance=1e-3)
   expect_equal(diag(cor(res.fcca.cv$Py, res3$final.model.cv.Py)),
      rep(1, ndim), tolerance=1e-3)

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

