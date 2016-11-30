#' Principal Component Analysis using flashpca
#'
#' @param X A numeric matrix to perform PCA on. X must not have any NAs.
#'
#' @param stand A character string indicating how to standardize X before PCA,
#' one of "binom" (old Eigenstrat-style), "binom2" (new Eigenstrat-style),
#' "sd" (zero-mean unit-variance), "center" (zero mean), or "none".
#'
#' @param ndim Integer. How many dimensions to return in results.
#' 
#' @param maxiter Integer. How many iterations to use in PCA.
#'
#' @param tol Numeric. Tolerance for determining PCA convergence.
#'
#' @param seed Integer. Seed for random number generator.
#'
#' @param verbose logical. More verbose output.
#'
#' @param do_loadings Logical. Whether to compute loadings (matrix V from SVD).
#'
#' @param check_geno Logical. Whether to explicitly check if the matrix X
#' contains values other than {0, 1, 2}, when stand="binom". This can be
#' be set to FALSE if you are sure your matrix only contains these values
#' (only matters when using stand="binom").
#' 
#' @param return_scale Logical. Whether to return the means and standard
#' deviations used in standardizing the matrix X.
#'
#' @details
#'    The decomposition is of X X' / m, where m is the number of SNPs.
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
#'    \item{values}{a numeric vector. The eigenvalues of X X' / m.}
#'    \item{vectors}{a numeric matrix. The eigenvectors of X X' / m.}
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
#' f1 <- flashpca(x, stand="sd", ndim=ndim)
#' 
#' r <- prcomp(x, center=TRUE, scale.=TRUE)
#' 
#' ## Compare eigenvalues
#' eval <- cbind(r$sdev[1:ndim]^2, f1$values)
#' cor(eval)
#' mean((eval[,1] - eval[,2])^2)
#' 
#' ## Compare eigenvectors
#' cor(r$x[, 1:ndim], f1$projection)
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
#'    f <- flashpca(hapmap3$bed, stand="center", ndim=ndim)
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
flashpca <- function(X, ndim=10,
   stand=c("binom2", "binom", "sd", "center", "none"),
   maxiter=1e2, tol=1e-4, seed=1, verbose=FALSE,
   do_loadings=FALSE, check_geno=TRUE, return_scale=TRUE)
{
   stand <- match.arg(stand)
   
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

   if(is.numeric(X)) {
      maxdim <- min(dim(X))
      ndim <- min(maxdim, ndim)
   }

   # otherwise Rcpp will throw an exception
   if(is.numeric(X)) {
      storage.mode(X) <- "numeric"
   }

   res <- try(
      if(is.character(X)) {
	 flashpca_plink_internal(X, stand_i, ndim, maxiter,
	    tol, seed, verbose, do_loadings, return_scale)
      } else {
	 flashpca_internal(X, stand_i, ndim, maxiter,
	    tol, seed, verbose, do_loadings, return_scale)
      }
   )
   if(is(res, "try-error")) {
      NULL
   } else {
      res
   }
}

