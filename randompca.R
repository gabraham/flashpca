
normalize <- function(x)
{
   r <- apply(x, 2, function(m) sqrt(crossprod(m)))
   sweep(x, 2, r, FUN="/")
}

randompca <- function(X, R=NULL, ndim=10, nextra=10, maxiter=500,
   whiten=FALSE, center=TRUE, scale=FALSE, tol=1e-6)
{
   M <- scale(X, center=center, scale=scale)
   total_dim <- ndim + nextra

   if(is.null(R)) {
      R <- matrix(rnorm(total_dim * ncol(M)), ncol=total_dim)
   }
   Y <- M %*% R
   Y <- normalize(Y)
   for(iter in 1:maxiter) {
      cat("iter", iter)
      Yn <- M %*% crossprod(M, Y)
      Yn <- normalize(Yn)
      diff <- mean((Y - Yn)^2)
      cat(" ", diff, "\n")
      Y <- Yn
      if(diff < tol) {
	 break
      }
   }
   
   q <- qr(Y)
   Q <- qr.Q(q)
   
   B <- crossprod(Q, M)
   s <- svd(B)
   
   U <- (Q %*% s$u)[, 1:ndim]
   D <- s$d[1:ndim]
   V <- s$v[, 1:ndim]

   X2 <- M %*% V

   # ZCA whitening
   Z <- U %*% tcrossprod(diag(1 / D), U)
   W <- Z %*% X

   list(U=U, D=D, V=V, X=X2, W=W)
}

