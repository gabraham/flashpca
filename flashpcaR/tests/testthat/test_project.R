context("Testing PCA projection onto new samples")

tol <- 1e-5

data(hm3.chr1)

bedf <- gsub("\\.bed", "",
   system.file("extdata", "data_chr1.bed", package="flashpcaR"))

ndim <- 10

test_that("Testing projection", {

   # PCA on all the data
   X1 <- scale2(hm3.chr1$bed, type="2")
   f <- flashpca(X1, ndim=ndim, stand="none", do_loadings=TRUE)
   P1 <- X1 %*% f$loadings / sqrt(ncol(X1))

   refallele <- hm3.chr1$bim[,5]
   names(refallele) <- hm3.chr1$bim[,2]

   # First project the same data used for PCA, onto the PCs   
   pr1 <- project(bedf, loadings=f$loadings, ref_allele=refallele,
      orig_mean=attr(X1, "scaled:center"), orig_sd=attr(X1, "scaled:scale"))
   expect_equal(f$projection, pr1$projection, tol=tol)

   expect_warning(pr2 <- project(hm3.chr1$bed, loadings=f$loadings, ref_allele=refallele,
      orig_mean=attr(X1, "scaled:center"), orig_sd=attr(X1, "scaled:scale")))
   expect_equal(as.numeric(f$projection), as.numeric(pr2$projection), tol=tol)

   # Do PCA on a subset, then project on everyone (not what we
   # would normally do but still good to test we get the right results)
   w <- sample(c(TRUE, FALSE), nrow(X1), prob=c(0.5, 0.5), replace=TRUE)
   X2 <- scale2(hm3.chr1$bed[!w,])

   f2 <- flashpca(X2, ndim=ndim, stand="none", do_loadings=TRUE)
   pr2 <- project(bedf, loadings=f2$loadings, ref_allele=refallele,
      orig_mean=attr(X2, "scaled:center"), orig_sd=attr(X2, "scaled:scale"))
   X1s <- scale(hm3.chr1$bed,
      center=attr(X2, "scaled:center"),
      scale=attr(X2, "scaled:scale"))
   X1s[is.na(X1s)] <- 0
   P2 <- X1s %*% f2$loadings / sqrt(ncol(X2))
   rownames(P2) <- NULL
   dimnames(P2) <- NULL
   expect_equal(P2, pr2$projection, tol=tol)
})

test_that("Testing projection input checking", {
   # PCA on all the data
   X1 <- scale2(hm3.chr1$bed, type="2")
   f <- flashpca(X1, ndim=ndim, stand="none", do_loadings=TRUE)

   refallele <- hm3.chr1$bim[,5]
   names(refallele) <- hm3.chr1$bim[,2]

   # Incorrect ref alleles
   expect_error(
      project(bedf, loadings=f$loadings, ref_allele=sample(refallele),
	 orig_mean=attr(X1, "scaled:center"),
	 orig_sd=attr(X1, "scaled:scale"))
   )

   # Length mismatch for loadings
   expect_error(
      project(bedf, loadings=f$loadings[1:10], ref_allele=refallele,
	 orig_mean=attr(X1, "scaled:center"),
	 orig_sd=attr(X1, "scaled:scale"))
   )

   # Length mismatch for orig_mean
   expect_error(
      project(bedf, loadings=f$loadings, ref_allele=refallele,
	 orig_mean=attr(X1, "scaled:center")[1:10],
	 orig_sd=attr(X1, "scaled:scale"))
   )

   # Length mismatch for orig_sd
   expect_error(
      project(bedf, loadings=f$loadings, ref_allele=refallele,
	 orig_mean=attr(X1, "scaled:center"),
	 orig_sd=attr(X1, "scaled:scale")[1:10])
   )

   # Negative or zero orig_sd
   osd <- attr(X1, "scaled:scale")
   osd[1] <- 0
   osd[2] <- -1
   expect_error(
      project(bedf, loadings=f$loadings, ref_allele=refallele,
	 orig_mean=attr(X1, "scaled:center"), orig_sd=osd)
   )
})

