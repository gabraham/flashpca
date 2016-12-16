#' Project new data onto existing principal components
#'
#' @param X A numeric matrix to project onto the PCs, or a
#' character string pointing to a PLINK dataset. 
#' 
#' @param loadings A numeric matrix of right
#' eigenvectors (SNPs on rows, ndim dimensions on columns).
#'
#' @param orig_mean A numeric vector of the original SNP means (averages)
#' used to standardise the genotypes before PCA.
#'
#' @param orig_sd A numeric vector of the original SNP standard deviations
#' used to standardise the genotypes before PCA.
#'
#' @param ref_alleles A character vector of the reference alleles for each
#' SNP. This must match the reference alleles used for PCA.
#'
#' @param divisor A character string indicating whether to divide the
#' eigenvalues by number of columns of X ("p"), the number of
#' rows of X minus 1 ("n1") or none ("none").
#' 
#' @param block_size Integer. Block size for PCA on PLINK files.
#' 
#' @param verbose Logical. Verbose output.
#' 
#' @param check_geno Logical. Whether to explicitly check if the matrix X
#' contains values other than {0, 1, 2}, when stand="binom". This can be
#' be set to FALSE if you are sure your matrix only contains these values
#' (only matters when using stand="binom").
#'
#' @param check_bim Logical. Whether to check that the number of rows in 
#' the PLINK bim file (if X is a character string) matches the number of
#' rows in the loadings.
#' 
#' @details
#'
#' \code{project} projects data Z onto existing principal components (PCs), i.e.,
#' \deqn{P_z = Z V} where V is the matrix of SNPs loadings. The projected data
#' must match the original data used to derive the PCs in terms of
#' \describe{
#'   \item{}{Same SNPs, in the same order}
#'   \item{}{Same alleles and reference alleles}
#'   \item{}{Same divisor specified (although not strictly necessary}.
#'}
#'
#' @return \code{project} returns a list containing the following components:
#' \describe{  
#'    \item{projection}{A numeric matrix. The projections of the new data
#'   onto the principal components.}
#' }
#'
#' @importFrom utils read.table
#' @export
project <- function(X, loadings, orig_mean=NULL, orig_sd=NULL,
   ref_alleles=NULL, divisor="p", block_size=1000, verbose=FALSE,
   check_geno=TRUE, check_bim=TRUE)
{
   divisor <- match.arg(divisor)

   if(is.null(orig_mean)) {
      stop("The vector of means used for standardising",
	 " the data must be provided via 'orig_mean'")
   } else if(is.null(orig_sd)) {
      stop("The vector of standard deviations used for standardising",
	 " the data must be provided via 'orig_sd'")
   }

   if(is.numeric(X)) {
      if(any(is.na(X))) {
	 warning("X contains missing values, will be mean imputed")
      }

      if(nrow(loadings) != ncol(X)) {
	 stop("The number of rows in X and number ",
	    "of columns of the loadings don't match")
      } else if(length(orig_mean) != ncol(X)) {
	 stop("The number of rows in X and length ",
	    "of orig_mean don't match")
      } else if(length(orig_sd) != ncol(X)) {
	 stop("The number of rows in X and length ",
	    "of orig_sd don't match")
      } else if(any(orig_sd <= 0)) {
	 stop("orig_sd cannot be zero or negative")
      }

      X <- scale(X, center=orig_mean, scale=orig_sd)
      for(j in 1:ncol(X)) {
	 w <- is.na(X[,j])
	 if(any(w)) {
	    X[w, j] <- 0
	 }
      }
      n <- nrow(X)
   } else if(is.character(X)) {
      #if(!stand %in% c("binom", "binom2")) {
      #   stop("When using PLINK data, ",
      #      "you must use stand='binom' or 'binom2'")
      #}
      if(check_bim) {
	 bim <- read.table(paste0(X, ".bim"), header=FALSE, sep="",
	    stringsAsFactors=FALSE)
	 if(nrow(bim) != nrow(loadings)) {
	    stop("The number of rows in ", X, ".bim",
	       " and the number of columns in the loadings don't match")
	 } else if(any(bim[,2] != names(ref_alleles))) {
	    stop("The SNP names in ", X, ".bim",
	       "do not match the names of the ref_alleles vector")
	 } else if(any(bim[,5] != ref_alleles)) {
	    stop("The reference alleles in ", X, ".bim",
	       "do not match the ref_alleles vector")
	 } else if(length(orig_mean) != nrow(bim)){
	    stop("The number of rows in ", X, ".bim",
	       " and the length of orig_mean don't match")
	 } else if(length(orig_sd) != nrow(bim)){
	    stop("The number of rows in ", X, ".bim",
	       " and the length of orig_sd don't match")
	 } else if(any(orig_sd <= 0)) {
	    stop("orig_sd cannot be zero or negative")
	 }
	 rm(bim)
	 fam <- read.table(paste0(X, ".fam"), header=FALSE, sep="",
	    stringsAsFactors=FALSE)
	 n <- nrow(fam)
	 rm(fam)
      }
   } else {
      stop("X must be a numeric matrix or a string naming a PLINK fileset")
   }

   divisors <- c(
      "p"=2,
      "n1"=1,
      "none"=0
   )
   div <- divisors[divisor]

   divisors_val <- c(
      "p"=nrow(loadings),
      "n1"=n,
      "none"=1
   )
   div_val <- divisors_val[divisor]

   # If the matrix is integer, Rcpp will throw an exception
   if(is.numeric(X)) {
      storage.mode(X) <- "numeric"
   }

   res <- try(
      if(is.character(X)) {
	 project_plink_internal(X, loadings, ref_alleles,
	    orig_mean, orig_sd, block_size, div, verbose)
      } else {
	 list(projection=X %*% loadings / sqrt(div_val))
      }
   )
   class(res) <- "flashpca.projection"
   if(is(res, "try-error")) {
      NULL
   } else {
      res
   }
}

