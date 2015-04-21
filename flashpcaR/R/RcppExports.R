
flashpca <- function(X, method=c("eigen", "svd"),
   stand=c("binom", "sd", "center", "none"), transpose=NULL, ndim=10,
   nextra=10, maxiter=1e2, tol=1e-6, seed=1, kernel=c("linear", "rbf"),
   sigma=NULL, rbf_center=TRUE, rbf_sample=1000, save_kernel=FALSE,
   do_orth=TRUE, verbose=FALSE, num_threads=1, do_loadings=FALSE,
   mem=c("low", "high"))
{
   method <- match.arg(method)
   stand <- match.arg(stand)
   kernel <- match.arg(kernel)
   mem <- match.arg(mem)

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
      transpose <- FALSE
   }
   maxdim <- min(dim(X))
   ndim <- min(maxdim, ndim)
   nextra <- min(maxdim - ndim, nextra)
   if(is.null(sigma)) {
      sigma <- 0 # will be estimated from data
   }

   rbf_sample <- min(nrow(X), rbf_sample)

   if(mem == "high") {
      mem_i <- 2L
   } else {
      mem_i <- 1L
   }

   # otherwise Rcpp will throw an exception
   storage.mode(X) <- "numeric"

   .Call("flashpca", PACKAGE="flashpcaR",
      X, method_i, stand_i, transpose, ndim, nextra, maxiter,
      tol, seed, kernel_i, sigma, rbf_center, rbf_sample,
      save_kernel, do_orth, verbose, num_threads, do_loadings, mem_i)
}

scca <- function(X, Y, lambda1=0, lambda2=0,
   stand=c("binom", "sd", "center", "none"), ndim=10,
   maxiter=1e3, tol=1e-6, seed=1L, verbose=FALSE, num_threads=1,
   mem=c("low", "high"))
{
   stand <- match.arg(stand)
   mem <- match.arg(mem)

   X <- cbind(X)
   Y <- cbind(Y)

   if(nrow(X) != nrow(Y)) {
      stop("X and Y must have the same number of rows")
   }

   if(!is.numeric(X) || !is.numeric(Y)) {
      stop("X and Y must both be numeric matrices")
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

   maxdim <- min(dim(X))
   ndim <- min(maxdim, ndim)

   if(mem == "high") {
      mem_i <- 2L
   } else {
      mem_i <- 1L
   }

   # otherwise Rcpp will throw an exception
   storage.mode(X) <- "numeric"
   storage.mode(Y) <- "numeric"

   res <- .Call("scca", PACKAGE="flashpcaR",
      X, Y, lambda1, lambda2, ndim, stand_i,
      mem_i, seed, num_threads, maxiter, tol, verbose)
   res
}


