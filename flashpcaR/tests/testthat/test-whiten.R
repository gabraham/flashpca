context("Testing whitening")

data(hm3.chr1)

bedf <- gsub("\\.bed", "",
   system.file("extdata", "data_chr1.bed", package="flashpcaR"))

test_that("Testing rank validation", {
   X <- scale2(hm3.chr1$bed)
   maxdim <- 50
   r <- validate.rank(X, maxdim=maxdim)

   expect_true(all(!is.na(r$mse)))
   expect_true(all(r$mse > 0))
   expect_true(all(r$dim > 0))
   expect_true(all(r$dim <= maxdim))
})

test_that("Testing whitening", {
   X <- scale2(hm3.chr1$bed)
   maxdim <- 50
   r <- validate.rank(X, maxdim=maxdim)
   mx <- r$dim[which.min(r$mse)]

   Xw <- whiten(X, ndim=mx)
   expect_true(all(!is.na(Xw)))
   s <- flashpca(Xw, stand="sd", ndim=maxdim + 10)

   expect_true(all(s$values[1:mx] > 0))
   expect_true(all(s$values[(mx+1):(maxdim+10)] > 0))
})

