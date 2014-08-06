
flashpca <- function(X, method=c("eigen", "svd"),
   stand=c("binom", "sd", "center", "none"), transpose=NULL, ndim=10,
   nextra=10, maxiter=50, tol=1e-6, seed=1, kernel=c("linear", "rbf"),
   sigma=NULL, rbf_center=TRUE, rbf_sample=1000, save_kernel=FALSE,
   do_orth=TRUE, verbose=FALSE, num_threads=1, do_loadings=FALSE)
{
   method <- match.arg(method)
   stand <- match.arg(stand)
   kernel <- match.arg(kernel)

   if(method == "eigen") {
      method_i <- 1L
   } else {
      method_i <- 2L
   }

   if(stand == "none") {
      stand_i <- 0L
   } else if(stand == "sd") {
      stand_i <- 1L
   } else if(stand == "binom") {
      stand_i <- 2L
   } else if(stand == "center") {
      stand_i <- 3L
   } 

   if(kernel == "linear") {
      kernel_i <- 1L
   } else {
      kernel_i <- 2L
   }

   # We only transpose if X is transposed, i.e., if X is SNPs by samples
   if(is.null(transpose)) {
      tranpose <- FALSE
   }
   maxdim <- min(dim(X))
   ndim <- min(maxdim, ndim)
   nextra <- min(maxdim - ndim, nextra)
   if(is.null(sigma)) {
      sigma <- 0 # will be estimated from data
   }

   rbf_sample <- min(nrow(X), rbf_sample)

   # otherwise Rcpp will throw an exception
   storage.mode(X) <- "numeric"

   .Call("flashpca", PACKAGE="flashpcaR",
      X, method_i, stand_i, transpose, ndim, nextra, maxiter,
      tol, seed, kernel_i, sigma, rbf_center, rbf_sample,
      save_kernel, do_orth, verbose, num_threads, do_loadings)
}

