
#' Cross-validation for principal component canonical correlation analysis (PCCA)
#'
#' @param X An n by p numeric matrix. Note: PLINK data is currently not supported.
#'
#' @param Y An n by k numeric matrix.
#'
#' @param ndim Integer. Positive number of canonical vectors to compute.
#'
#' @param kx Numeric. Number of principal components to compute for X.
#'
#' @param ky Numeric. Number of principal components to compute for Y.
#'
#' @param nfolds Integer. Number of cross-validation folds.
#'
#' @param folds Integer. The fold identifiers. Overrides the `nfolds'
#' parameter.
#'
#' @param standx Character. One of "sd" (empricial standard deviation) or "none".
#'
#' @param standy Character. One of "sd" (empricial standard deviation) or "none".
#'
#' @param svd.tol Numeric. Tolerance under which to truncate singular values of X and Y.
#'
#' @importFrom stats cancor
#'
#' @export
cv.pcca <- function(X, Y, ndim=1, kx=3, ky=3, nfolds=3, folds=NULL,
   standx=c("sd", "none"), standy=c("sd", "none"), svd.tol=1e-12)
{
   standx <- match.arg(standx)
   standy <- match.arg(standy)

   if(is.null(X) || !is.numeric(X) || ncol(X) <= 1) {
      stop("X must be a numeric matrix")
   }

   if(is.null(Y) || !is.numeric(Y) || ncol(Y) <= 1) {
      stop("Y must be a numeric matrix")
   }

   if(nrow(X) != nrow(Y)) {
      stop("X and Y must have the same number of rows")
   }

   if(!is.numeric(kx) || kx < 1 || kx > ncol(X)) {
      stop("kx must be >0 and < ncol(X)")
   }

   if(!is.numeric(ky) || ky < 1 || ky > ncol(X)) {
      stop("ky must be >0 and < ncol(X)")
   }

   if(standx == "sd") {
      X <- scale(X)
   }
   if(standy == "sd") {
      Y <- scale(Y)
   }

   if(!is.null(folds)) {
      nfolds <- max(folds)
   } else {
      if(nfolds < 2 || nfolds >= nrow(X)) {
	 stop("nfolds must be >1 and < nrow(X)") 
      }
      folds <- sample(nfolds, nrow(X), replace=TRUE)
   }

   ndim.max <- min(kx, ky)

   Px.tst <- matrix(0, nrow(X), ndim.max)
   Py.tst <- matrix(0, nrow(X), ndim.max)
   colnames(Px.tst) <- paste0("Px", 1:ndim.max)
   colnames(Py.tst) <- paste0("Py", 1:ndim.max)

   for(fold in 1:nfolds) {
      fx <- svd(X[folds != fold,])
      fy <- svd(Y[folds != fold,])

      kxf <- min(kx, sum(fx$d > svd.tol))
      kyf <- min(ky, sum(fy$d > svd.tol))

      ndim <- min(kx, ky)
   
      # On the training data
      rc <- cancor(fx$u[, 1:kxf], fy$u[, 1:kyf], xcenter=FALSE, ycenter=FALSE)
   
      Ux.tst <- X[folds == fold, ] %*% fx$v[, 1:kxf] %*% diag(1 / fx$d[1:kxf])
      Uy.tst <- Y[folds == fold, ] %*% fy$v[, 1:kyf] %*% diag(1 / fy$d[1:kyf])

      Px.tst[folds == fold, ] <- Ux.tst %*% rc$xcoef[, 1:ndim]
      Py.tst[folds == fold, ] <- Uy.tst %*% rc$ycoef[, 1:ndim]
   }
   r.tst <- diag(cor(Px.tst, Py.tst))

   res <- list(ndim=ndim, nfolds=nfolds, folds=folds,
      kx=kx, ky=ky, Px=Px.tst, Py=Py.tst, r=r.tst)
   class(res) <- "cv.pcca"
   res
}

#' Principal component canonical correlation analysis (PCCA)
#'
#' @param X An n by p numeric matrix. Note: PLINK data is currently not supported.
#'
#' @param Y An n by k numeric matrix.
#'
#' @param ndim Integer. Positive number of canonical vectors to compute.
#'
#' @param kx Numeric. Number of principal components to compute for X.
#'
#' @param ky Numeric. Number of principal components to compute for Y.
#'
#' @param nfolds Integer. Number of cross-validation folds.
#'
#' @param folds Integer. The fold identifiers. Overrides the `nfolds'
#' parameter.
#'
#' @param standx Character. One of "sd" (empricial standard deviation) or "none".
#'
#' @param standy Character. One of "sd" (empricial standard deviation) or "none".
#'
#' @param svd.tol Numeric. Tolerance under which to truncate singular values of X and Y.
#'
#' @importFrom stats cancor
#'
#' @export
pcca <- function(X, Y, ndim=1, kx=3, ky=3,
   standx=c("sd", "none"), standy=c("sd", "none"), svd.tol=1e-12)
{
   standx <- match.arg(standx)
   standy <- match.arg(standy)

   if(is.null(X) || !is.numeric(X) || ncol(X) <= 1) {
      stop("X must be a numeric matrix")
   }

   if(is.null(Y) || !is.numeric(Y) || ncol(Y) <= 1) {
      stop("Y must be a numeric matrix")
   }

   if(nrow(X) != nrow(Y)) {
      stop("X and Y must have the same number of rows")
   }

   if(!is.numeric(kx) || kx < 1 || kx > ncol(X)) {
      stop("kx must be >0 and < ncol(X)")
   }

   if(!is.numeric(ky) || ky < 1 || ky > ncol(X)) {
      stop("ky must be >0 and < ncol(X)")
   }

   if(standx == "sd") {
      X <- scale(X)
   }
   if(standy == "sd") {
      Y <- scale(Y)
   }

   fx <- svd(X)
   fy <- svd(Y)
   kx <- min(kx, sum(fx$d > svd.tol))
   ky <- min(ky, sum(fy$d > svd.tol))

   ndim <- min(kx, ky)
   
   rc <- cancor(fx$u[, 1:kx], fy$u[, 1:ky], xcenter=FALSE, ycenter=FALSE)
   
   Px <- fx$u[, 1:kx] %*% rc$xcoef[, 1:ndim]
   Py <- fy$u[, 1:ky] %*% rc$ycoef[, 1:ndim]
   colnames(Px) <- paste0("Px", 1:ndim)
   colnames(Py) <- paste0("Py", 1:ndim)
   
   r <- diag(cor(Px, Py))

   res <- list(ndim=ndim,
      kx=kx, ky=ky, Px=Px, Py=Py, r=r)
   class(res) <- "pcca"
   res
}

#' Prints a cv.pcca object
#'
#' @param x A cv.pcca object to be printed
#' @param ... Ignored
#' @export 
print.cv.pcca <- function(x, ...)
{
   cat(paste0(
      "cv.pcca object; ndim=", x$ndim, "; ", x$nfolds,
      "-fold cross-validation\n"))
   invisible(x)
}

#' Prints a pcca object
#'
#' @param x A pcca object to be printed
#' @param ... Ignored
#' @export 
print.pcca <- function(x, ...)
{
   cat(paste0(
      "pcca object; ndim=", x$ndim, "\n"))
   invisible(x)
}
