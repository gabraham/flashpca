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
#' where variance is p * (1 - p), for SNP data), "binom2" (zero mean, unit
#' variance where variance is 2 * p * (1 - p),
#' "sd" (zero-mean, unit Gaussian variance), "center" (zero mean), or "none". Note
#' that if you use "binom" for non-SNP data you may get garbage. Note that
#' currently the same standardization is applied to both X and Y. If you
#' require different standardizations, it's best to standardize yourself
#' and then choose stand="none".
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
#' @param check_geno Logical. Whether to explicitly check if the matrices
#' X and Y contain values other than {0, 1, 2}, when stand="binom". This can
#' be set to FALSE if you are sure your matrices only contain these values
#' (only matters when using stand="binom").
#' 
#' @param V Numeric. A vector to initialise "v" in SCCA iterations. By
#' default, it will be a vector of normally distributed variates.
#'
#' @return \code{scca} returns a list containing the following components:
#'
#' \describe{  
#'    \item{U}{Top ndim left canonical vectors of X^T Y.}
#'    \item{V}{Top ndim right canonical vectors of X^T Y.}
#'    \item{d}{Top ndim canonical covariances, i.e., diag(cov(X U, Y V)). 
#'	 Note that we don't divide by n-1.}
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
#' s <- scca(X, Y, lambda1=1e-2, lambda2=1e-2, ndim=5, stand="sd")
#'
#' ## The canonical correlations
#' diag(cor(s$Px, s$Py))
#'
#' @export
scca <- function(X, Y, lambda1=0, lambda2=0,
   stand=c("binom", "binom2", "sd", "center", "none"),
   ndim=10, maxiter=1e3, tol=1e-4, seed=1L, verbose=FALSE, num_threads=1,
   mem=c("low", "high"), check_geno=TRUE, V=NULL)
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

   if(stand == "binom" && check_geno) {
      wx <- X %in% 0:2
      wy <- Y %in% 0:2
      if(sum(wx) != length(X) || sum(wy) != length(Y)) {
	 stop(
	    "Your data contains values other than {0, 1, 2}, stand='binom' can't be used here")
      }
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

   if(!is.null(V)) {
      V <- cbind(V)
      if(nrow(V) < ncol(Y) || ncol(V) != ndim) {
         stop("dimensions of V must be (nrow(Y) x (ndim))")
      }
   }

   res <- try(
      if(is.null(V)) {
	 scca_internal(X, Y, lambda1, lambda2, ndim, stand_i, mem_i, seed, maxiter, tol,
	    verbose, num_threads, FALSE, matrix(0))
      } else {
	 scca_internal(X, Y, lambda1, lambda2, ndim, stand_i, mem_i, seed, maxiter, tol,
	    verbose, num_threads, TRUE, V)
      }
   )
   class(res) <- "scca"
   if(is(res, "try-error")) {
      NULL
   } else {
      res
   }
}

#' @param x A flashpca object to be printed
#' @export 
print.scca <- function(x, ...)
{
   cat("scca object; ndim=", length(x$d), "\n")
   invisible(x)
}

