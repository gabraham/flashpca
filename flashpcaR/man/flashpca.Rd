\name{flashpca}
\alias{flashpca}

\title{Principal Component Analysis using flashpca} 

\description{
   Performs a principal component analysis on the matrix X and returns the
   eigenvalues and eigenvectors of X' X / (n - 1).
}

\details{
}

\usage{
   flashpca(X, ...)

   \method{flashpca}{default}(X, method = c("eigen", "svd"),
      stand = c("binom", "sd", "center", "none"),
      transpose = NULL, ndim = 10, nextra = 10,
      maxiter = 50, tol = 1e-06, seed = 1,
      kernel = c("linear", "rbf"), sigma = NULL, rbf_center = TRUE,
      rbf_sample = 1000, save_kernel = FALSE, do_orth = TRUE,
      verbose = FALSE, num_threads = 1, do_loadings = FALSE,
      mem = c("high", "low"))
   
} 

\arguments{

\item{X}{a numeric matrix to perform PCA on. X must not have any NAs.}
\item{method}{a character string indicating which decomposition to use, one of "eigen"
(eigen-decomposition of X'X / (n - 1)) or "svd" (SVD of X).}
\item{stand}{a character string indicating how to standardize X before PCA,
one of "binom" (Eigenstrat-style), "sd" (zero-mean unit-variance), "center"
(zero mean), or "none".}
\item{transpose}{logical. If X is transposed (samples are on columns rather
than rows), set to \code{TRUE}.}
\item{ndim}{integer. How many dimensions to return in results.}
\item{nextra}{integer. How many extra dimensions to use for doing PCA.
Increasing this will increase accuracy.}
\item{maxiter}{integer. How many iterations to use in PCA.}
\item{tol}{numeric. Tolerance for determining PCA convergence.}
\item{seed}{integer. Seed for random number generator.}
\item{kernel}{a character string, one of "linear" or "rbf", indicating whether
to do standard (linear) PCA or use a radial basis function (RBF) kernel.}
\item{sigma}{numeric. The sigma parameter for the RBF kernel.}
\item{rbf_center}{logical. Whether to center the data in feature space.}
\item{rbf_sample}{integer. How many observations will be randomly sampled to
determine the default RBF sigma (unless sigma is specified).}
\item{save_kernel}{logical. Save the kernel to disk? (hardcoded as "kernel.txt")}
\item{do_orth}{logical. Whether to re-orthogonoalize during each PCA step,
tends to increase accuracy.}
\item{verbose}{logical. More verbose output.}
\item{num_threads}{integer. How many OpenMP threads to use, if compiled with
OpenMP support.}
\item{do_loadings}{logical. Whether to compute loadings (matrix V from SVD).}
\item{mem}{a character string, one of "high" or "low". "High" requires
computing the X' X / (n - 1) covariance matrix which is memory intensive but fast if
using multiple threads. "Low" does not compute X' X and uses less memory,
but it tends to be slightly slower.}

} 

\value{
\code{flashpca} returns a list containing the following components:
  
  \item{values}{a numeric vector. The eigenvalues of X' X / (n - 1).}
  \item{labels}{a numeric matrix. The eigenvectors of X' X (n - 1), i.e., U from SVD.}
  \item{projection}{a numeric matrix. Equivalent to X V.} 
  \item{loadings}{a numeric matrix. The matrix of variable loadings, i.e., V
  from SVD.}

}

\examples{
set.seed(123)

## Toy example
n <- 200
p <- 500
x <- matrix(rnorm(n * p), n, p)
k <- 20
f1 <- flashpca(x, stand="sd", ndim=k, method="eigen")
f2 <- flashpca(x, stand="sd", ndim=k, method="svd")

r <- prcomp(x, center=TRUE, scale.=TRUE)

## Compare eigenvalues
eval <- cbind(r$sdev[1:k]^2, f1$values, f2$values)
cor(eval)
mean((eval[,1] - eval[,2])^2)
mean((eval[,1] - eval[,3])^2)

## Compare eigenvectors
diag(cor(r$x[, 1:k], f1$projection))
diag(cor(r$x[, 1:k], f2$projection))

\dontrun{
## First get data from
## https://github.com/gabraham/flashpca/blob/devel/HapMap3/data.RData
##
## Missing genotypes have been imputed randomly according to genotype
## proportions in each SNP.
load("data.RData")
k <- 20
system.time({
   f <- flashpca(hapmap3$bed, stand="center", ndim=k, nextra=50)
})
system.time({
   r <- prcomp(hapmap3$bed)
})

eval <- cbind(r$sdev[1:k]^2, f$values)
cor(eval)
mean((eval[,1] - eval[,2])^2)

## Compare eigenvectors
diag(cor(r$x[, 1:k], f$projection))
}

}

