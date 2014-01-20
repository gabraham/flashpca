
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
   D <- s$d[1:ndim] / sqrt(nrow(X) - 1)
   V <- s$v[, 1:ndim]

   X2 <- M %*% V

   # ZCA whitening
   Z <- U %*% tcrossprod(diag(1 / D), U)
   W <- Z %*% X

   list(U=U, D=D, V=V, X=X2, W=W)
}

randompca2 <- function(X, R=NULL, ndim=10, nextra=10, maxiter=500,
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

   BBT <- tcrossprod(B)
   e <- eigen(BBT)

   U <- (Q %*% e$vectors)[, 1:ndim]
   D <- sqrt(e$values[1:ndim]) / sqrt(nrow(X) - 1)
   UD <- U %*% diag(D)

   # ZCA whitening
   #Z <- U %*% tcrossprod(diag(1 / D), U)
   #W <- Z %*% X

   list(U=U, D=D, X=UD)#, W=W)
}

n <- 100
p <- 5000
X <- matrix(rnorm(n * p), n, p)

system.time({
   r1 <- randompca(X, nextra=100)
})

system.time({
   r2 <- randompca2(X, nextra=100)
})

diag(cor(r1$X, r2$X))

