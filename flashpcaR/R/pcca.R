
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
#' @param check_sign Logical. Whether to check and correct possible sign flips 
#' of the singular vectors within the cross-validation.
#'
#' @param final_model Logical. Whether to return a \code{pcca} model trained
#' on all the data.
#'
#' @param svd_tol Numeric. Tolerance under which to truncate singular values of X and Y.
#'
#' @importFrom stats cancor
#'
#' @export
cv.pcca <- function(X, Y, ndim=NULL, kx=3, ky=3, nfolds=3, folds=NULL,
   standx=c("sd", "none"), standy=c("sd", "none"), check_sign=TRUE,
   final_model=FALSE, svd_tol=1e-12, verbose=FALSE)
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

   if(is.null(ndim)) {
      ndim.max <- min(kx, ky)
   } else {
      ndim.max <- min(ndim, kx, ky)
   }

   Px.tst <- matrix(0, nrow(X), ndim.max)
   Py.tst <- matrix(0, nrow(X), ndim.max)
   colnames(Px.tst) <- paste0("Px", 1:ndim.max)
   colnames(Py.tst) <- paste0("Py", 1:ndim.max)

   A.list <- vector("list", nfolds)
   B.list <- vector("list", nfolds)

   for(fold in 1:nfolds) {
      Xtrn <- X[folds != fold,]
      Ytrn <- Y[folds != fold,]
      fx <- svd(Xtrn)
      fy <- svd(Ytrn)

      kxf <- min(kx, sum(fx$d > svd_tol))
      kyf <- min(ky, sum(fy$d > svd_tol))

      ndimf <- min(kxf, kyf, ndim.max)
   
      # On the training data
      rc <- cancor(fx$u[, 1:kxf], fy$u[, 1:kyf], xcenter=FALSE, ycenter=FALSE)
   
      A.list[[fold]] <- fx$v[, 1:kxf] %*% diag(1 / fx$d[1:kxf]) %*% rc$xcoef[, 1:ndimf]
      B.list[[fold]] <- fy$v[, 1:kyf] %*% diag(1 / fy$d[1:kyf]) %*% rc$ycoef[, 1:ndimf]

      if(check_sign) {
	 if(verbose) {
	    cat("checking sign\n")
	 }
	 sc <- check.eig.sign(A.list[[fold]], B.list[[fold]], Xtrn, Ytrn)
	 A.list[[fold]] <- sc$A
	 B.list[[fold]] <- sc$B
      }
      Px.tst[folds == fold, ] <- X[folds == fold,] %*% A.list[[fold]]
      Py.tst[folds == fold, ] <- Y[folds == fold,] %*% B.list[[fold]]
   }
   r.tst <- diag(cor(Px.tst, Py.tst))

   mod <- NULL
   if(final_model) {
      mod <- pcca(X, Y, ndim=ndim,
	 kx=kx, ky=ky, standx=standx, standy=standy, check_sign=check_sign,
	 svd_tol=svd_tol)
   }

   res <- list(ndim=ndim.max, nfolds=nfolds, folds=folds,
      kx=kx, ky=ky, final_model=mod,
      final_model_cv_Px=Px.tst, final_model_cv_Py=Py.tst,
      final_model_cv_r=r.tst, final_model_cv_A=A.list, final_model_cv_B=B.list)
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
#' @param check_sign Logical. Whether to check and correct possible sign flips 
#' of the singular vectors within the cross-validation.
#'
#' @param svd_tol Numeric. Tolerance under which to truncate singular values of X and Y.
#'
#' @importFrom stats cancor
#'
#' @export
pcca <- function(X, Y, ndim=NULL, kx=3, ky=3,
   standx=c("sd", "none"), standy=c("sd", "none"), check_sign=TRUE,
   svd_tol=1e-12)
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
   kx <- min(kx, sum(fx$d > svd_tol))
   ky <- min(ky, sum(fy$d > svd_tol))

   if(is.null(ndim)) {
      ndim <- min(kx, ky)
   } else {
      ndim <- min(ndim, kx, ky)
   }
   
   rc <- cancor(fx$u[, 1:kx], fy$u[, 1:ky], xcenter=FALSE, ycenter=FALSE)
   
   A <- fx$v[, 1:kx] %*% diag(1 / fx$d[1:kx]) %*% rc$xcoef[, 1:ndim]
   B <- fy$v[, 1:ky] %*% diag(1 / fy$d[1:ky]) %*% rc$ycoef[, 1:ndim]

   #Px <- fx$u[, 1:kx] %*% rc$xcoef[, 1:ndim]
   #Py <- fy$u[, 1:ky] %*% rc$ycoef[, 1:ndim]
   
   if(check_sign) {
      sc <- check.eig.sign(A, B, X, Y)
      A <- sc$A
      B <- sc$B
   }

   Px <- X %*% A
   Py <- Y %*% B
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

