#' Principal Component Analysis using flashpca
#'
#' @param X A numeric matrix to perform PCA on. X must not have any NAs.
#'
#' @param method A character string indicating which decomposition to use, one of "eigen"
#' (eigen-decomposition of X'X / (n - 1)) or "svd" (SVD of X).
#'
#' @param stand A character string indicating how to standardize X before PCA,
#' one of "binom" (Eigenstrat-style), "sd" (zero-mean unit-variance), "center"
#' (zero mean), or "none".
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
#' computing the X' X / (n - 1) covariance matrix which is memory intensive but fast if
#' using multiple threads. "Low" does not compute X' X and uses less memory,
#' but it tends to be slightly slower.
#'
#' @return \code{flashpca} returns a list containing the following components:
#'
#' \describe{  
#'    \item{values}{a numeric vector. The eigenvalues of X' X / (n - 1).}
#'    \item{vectors}{a numeric matrix. The eigenvectors of X' X (n - 1), i.e., U from SVD.}
#'    \item{projection}{a numeric matrix. Equivalent to X V.}
#'    \item{loadings}{a numeric matrix. The matrix of variable loadings, i.e., V
#'    from SVD.}
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

   if(!is.numeric(X)) {
      stop("X must both a numeric matrix")
   }

   if(any(is.na(X))) {
      stop("X cannot contain any missing values")
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

#' Performs sparse canonical correlation analysis
#'
#' @param X An n by p numeric matrix
#'
#' @param Y An n by k numeric matrix
#' 
#' @param lambda1 Numeric. Non-negative L1 penalty on canonical vectors of X.
#'
#' @param lambda2 Numeric. Non-negative L1 penalty on canonical vectors of Y.
#"
#' @param stand Character. One of "binom" (zero mean, unit variance
#' where variance is 2*p*(1-p), for SNP data),
#' 	"sd" (zero-mean, unit variance), "center" (zero mean), or "none".
#'
#' @param ndim Integer. Positive number of canonical vectors to compute.
#'
#' @param maxiter Integer. Positive, maximum number of iterations to perform.
#'
#' @param tol Numeric. Tolerance for convergence.
#'
#' @param seed Integer. Random seed for initialisation of CCA iterations.
#'
#' @param verbose Logical.
#'
#' @param num_threads Integer. Number of OpenMP threads to
#' use (only supported under Linux and on Mac if compiled with GCC).
#'
#' @param mem Character. One of "low" or "high", where "low"
#' doesn't explicitly compute X^T Y, and "high" does.
#' This is useful for large X and Y.
#'
#' @return \code{scca} returns a list containing the following components:
#'
#' \describe{  
#'    \item{U}{Top ndim left canonical vectors of X^T Y.}
#'    \item{V}{Top ndim right canonical vectors of X^T Y.}
#'    \item{d}{Top ndim canonical correlations.}
#'    \item{Px}{X times U.}
#'    \item{Py}{Y times V.}
#' }
#'
#' @examples
#'
#' n <- 500
#' p <- 100
#' k <- 50
#' m <- 5
#' M <- matrix(rnorm(n * m), n, m)
#' Bx <- matrix(rnorm(m * p), m, p)
#' By <- matrix(rnorm(m * k), m, k)
#' X <- scale(M %*% Bx + rnorm(n * p))
#' Y <- scale(M %*% By + rnorm(n * k))
#' 
#' s <- scca(X, Y, lambda1=1e-3, lambda2=1e-3, ndim=5)
#'
#' @export
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

   if(any(is.na(X)) || any(is.na(Y))) {
      stop("X and Y cannot contain any missing values")
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


