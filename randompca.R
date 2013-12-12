
randompca <- function(X, ndim=10, nextra=10, maxiter=5,
   center=TRUE, scale=FALSE)
{
   M <- scale(X, center=center, scale=scale)
   total_dim <- ndim + nextra

   R <- matrix(rnorm(total_dim * ncol(M)), ncol=total_dim)
   Y <- M %*% R
   for(iter in 1:maxiter) {
      Y <- M %*% crossprod(M, Y)
      r <- apply(Y, 2, function(x) sqrt(crossprod(x)))
      Y <- sweep(Y, 2, r, FUN="/")
   }
   
   q <- qr(Y)
   Q <- qr.Q(q)
   
   B <- crossprod(Q, M)
   s <- svd(B)
   
   U <- (Q %*% s$u)[, 1:ndim]
   D <- s$d[1:ndim]
   V <- s$v[, 1:ndim]

   X2 <- M %*% V

   list(U=U, D=D, V=V, X=X2)
}

