#' Modified version of algorithm from
#' Bro, Acar, and Kolda ''Resolving the sign ambiguity in the singular value
#' decomposition'' (2008), J Chemometrics (22) 135-140.
#'
#' @param A Numeric matrix. Left canonical vectors (loadings).
#' @param B Numeric matrix. Right canonical vectors (loadings).
#' @param X Numeric matrix.
#' @param Y Numeric matrix.
#'
check.eig.sign <- function(A, B, X, Y) {
   if(nrow(A) != ncol(X)) {
      stop("nrow(A) must be identical to ncol(X)")
   }
   if(nrow(B) != ncol(Y)) {
      stop("nrow(A) must be identical to ncol(X)")
   }
   if(ncol(A) != ncol(B)) {
      stop("A and B must have the same number of columns")
   }
   if(nrow(X) != nrow(Y)) {
      stop("X and Y must have the same number of rows")
   }

   #mx <- crossprod(A, t(X))
   #my <- crossprod(B, t(Y))
   #sx <- rowSums(sign(mx) * mx^2)
   #sy <- rowSums(sign(my) * my^2)
   mx <- X %*% A
   my <- Y %*% B
   sx <- colSums(sign(mx) * mx^2)
   sy <- colSums(sign(my) * my^2)
   ndim <- ncol(A)

   for(j in 1:ndim) {
      if(!is.na(sx[j]) && !is.na(sy[j])) {
	 if(sign(sx[j]) != sign(sy[j])) {
      	    if(abs(sx[j]) < abs(sy[j])) {
      	       sx[j] <- -sx[j]
      	    } else {
      	       sy[j] <- -sy[j]
      	    }
      	 }
      	 A[,j] <- sign(sx[j]) * A[,j]
      	 B[,j] <- sign(sy[j]) * B[,j]
      }
   }

   list(A=A, B=B, sx=sx, sy=sy)
}

