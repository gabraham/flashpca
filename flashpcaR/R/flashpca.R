#' Principal Component Analysis using flashpca
#'
#' @param X A numeric matrix to perform PCA on. X must not have any NAs.
#'
#' @param method A character string indicating which decomposition to use, one of "eigen"
#' (eigen-decomposition of X'X / (n - 1) or X'X, or "svd" (SVD of X).
#'
#' @param stand A character string indicating how to standardize X before PCA,
#' one of "binom" (old Eigenstrat-style), "binom2" (new Eigenstrat-style),
#' "sd" (zero-mean unit-variance), "center" (zero mean), or "none".
#'
#' @param transpose Logical. If X is transposed (samples are on columns rather
#' than rows), set to \code{TRUE}.
#'
#' @param ndim Integer. How many dimensions to return in results.
#' 
#' @param nextra Integer. How many extra dimensions to use for doing PCA.
#' Increasing this will increase accuracy.
#'
#' @param maxiter Integer. How many iterations to use in PCA.
#'
#' @param tol Numeric. Tolerance for determining PCA convergence.
#'
#' @param seed Integer. Seed for random number generator.
#'
#' @param kernel A character string, one of "linear" or "rbf", indicating whether
#' to do standard (linear) PCA or use a radial basis function (RBF) kernel.
#'
#' @param sigma numeric. The sigma parameter for the RBF kernel.
#'
#' @param rbf_center Logical. Whether to center the data in feature space.
#'
#' @param rbf_sample integer. How many observations will be randomly sampled to
#' determine the default RBF sigma (unless sigma is specified).
#'
#' @param save_kernel logical. Save the kernel to disk? (hardcoded as "kernel.txt")
#'
#' @param do_orth logical. Whether to re-orthogonoalize during each PCA step,
#' tends to increase accuracy.
#'
#' @param verbose logical. More verbose output.
#'
#' @param num_threads Integer. How many OpenMP threads to use, if compiled with
#' OpenMP support.
#'
#' @param do_loadings Logical. Whether to compute loadings (matrix V from SVD).
#'
#' @param mem A character string, one of "high" or "low". "High" requires
#' computing the X' X covariance matrix which is memory intensive but fast if
#' using multiple threads. "Low" does not compute X' X and uses less memory,
#' but it tends to be slightly slower.
#'
#' @param check_geno Logical. Whether to explicitly check if the matrix X
#' contains values other than {0, 1, 2}, when stand="binom". This can be
#' be set to FALSE if you are sure your matrix only contains these values
#' (only matters when using stand="binom").
#' 
#' @param return_scale Logical. Whether to return the means and standard
#' deviations used in standardizing the matrix X.
#'
#' @param divide_n Logical. Whether to compute the eigendecomposition
#' of  X' X / (n - 1) or of  X' X.
#'
#'
#' @details
#'    The decomposition is either done on X'X or on X'X/(n-1), depending on
#'    the \code{divide_n} argument. 
#'
#'    The data is standardised by default. \code{stand = "binom"} uses the old Eigensoft
#'    (Price 2006) formulation of
#'	       \deqn{p_j = sum_{i=1}^n X_{i,j} / (2 * n)}
#'	       \deqn{mean_j = 2 * p}
#'	       \deqn{sd_j = sqrt(p * (1 - p)),}
#'    where j is the index for the SNP and i is the index for the sample.
#'    Alternatively, `stand = "binom2"' uses the newer formula, which is
#'    similar to the above except for
#'	       \deqn{sd_j = sqrt(2 * p * (1 - p)).}
#' 
#' @return \code{flashpca} returns a list containing the following components:
#' \describe{  
#'    \item{values}{a numeric vector. The eigenvalues of X' X / (n - 1) or X' X.}
#'    \item{vectors}{a numeric matrix. The eigenvectors 
#'	 of X' X (n - 1) or X' X, i.e., U from SVD.}
#'    \item{projection}{a numeric matrix. Equivalent to X V.}
#'    \item{loadings}{a numeric matrix. The matrix of variable loadings, i.e., V
#'    from SVD.}
#'    \item{scale}{a list of two elements, ``center'' and ''scale'', which was
#'	 used to standardize the input matrix X.}
#' }
#' 
#' @examples
#' 
#' set.seed(123)
#' 
#' ## Toy example
#' n <- 200
#' p <- 500
#' x <- matrix(rnorm(n * p), n, p)
#' ndim <- 20
#' f1 <- flashpca(x, stand="sd", ndim=ndim, method="eigen")
#' f2 <- flashpca(x, stand="sd", ndim=ndim, method="svd")
#' 
#' r <- prcomp(x, center=TRUE, scale.=TRUE)
#' 
#' ## Compare eigenvalues
#' eval <- cbind(r$sdev[1:ndim]^2, f1$values, f2$values)
#' cor(eval)
#' mean((eval[,1] - eval[,2])^2)
#' mean((eval[,1] - eval[,3])^2)
#' 
#' ## Compare eigenvectors
#' diag(cor(r$x[, 1:ndim], f1$projection))
#' diag(cor(r$x[, 1:ndim], f2$projection))
#'
#' \dontrun{
#' ## First get data from
#' ## https://github.com/gabraham/flashpca/blob/devel/HapMap3/data.RData
#' ##
#' ## Missing genotypes have been imputed randomly according to genotype
#' ## proportions in each SNP.
#' load("data.RData")
#' ndim <- 20
#' system.time({
#'    f <- flashpca(hapmap3$bed, stand="center", ndim=ndim, nextra=50)
#' })
#' system.time({
#'    r <- prcomp(hapmap3$bed)
#' })
#' 
#' eval <- cbind(r$sdev[1:ndim]^2, f$values)
#' cor(eval)
#' mean((eval[,1] - eval[,2])^2)
#' 
#' ## Compare eigenvectors
#' diag(cor(r$x[, 1:ndim], f$projection))
#'}
#'
#' @export
flashpca <- function(X, method=c("eigen", "svd"),
   stand=c("binom", "binom2", "sd", "center", "none"),
   transpose=NULL, ndim=10,
   nextra=10, maxiter=1e2, tol=1e-4, seed=1, kernel=c("linear", "rbf"),
   sigma=NULL, rbf_center=TRUE, rbf_sample=1000, save_kernel=FALSE,
   do_orth=TRUE, verbose=FALSE, num_threads=1, do_loadings=FALSE,
   mem=c("low", "high"), check_geno=TRUE, return_scale=TRUE,
   divide_n=TRUE)
{
   method <- match.arg(method)
   stand <- match.arg(stand)
   kernel <- match.arg(kernel)
   mem <- match.arg(mem)
   
   if(is.numeric(X)) {
      if(any(is.na(X))) {
	 stop("X cannot contain any missing values")
      }

      if(stand %in% c("binom", "binom2") && check_geno) {
	 wx <- X %in% 0:2
      	 if(sum(wx) != length(X)) {
      	    stop(
      	       "Your data contains values other than {0, 1, 2}, ",
      	       "stand='binom'/'binom2' can't be used here")
      	 }
      }
   } else if(!is.character(X)) {
      stop("X must be a numeric matrix or a string naming a PLINK fileset")
   }

   if(!is.logical(divide_n)) {
      stop("divide_n must be a logical (TRUE/FALSE)")
   }

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
   } else if(stand == "binom2") {
      stand_i <- 3L
   } else if(stand == "center") {
      stand_i <- 4L
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

   if(is.numeric(X)) {
      maxdim <- min(dim(X))
      ndim <- min(maxdim, ndim)
      nextra <- min(maxdim - ndim, nextra)
      if(is.null(sigma)) {
         sigma <- 0 # will be estimated from data
      }

      rbf_sample <- min(nrow(X), rbf_sample)
   }

   if(mem == "high") {
      mem_i <- 2L
   } else {
      mem_i <- 1L
   }

   # otherwise Rcpp will throw an exception
   if(is.numeric(X)) {
      storage.mode(X) <- "numeric"
   }

   res <- try(
      if(is.character(X)) {
	 flashpca_plink_internal(X, stand_i, ndim, maxiter, tol, verbose)
      } else {
	 flashpca_internal(X, method_i, stand_i, transpose, ndim, nextra, maxiter,
	    tol, seed, kernel_i, sigma, rbf_center, rbf_sample, save_kernel, do_orth,
	    verbose, do_loadings, mem_i, return_scale, num_threads, divide_n)
      }
   )
   if(is(res, "try-error")) {
      NULL
   } else {
      res
   }
}

