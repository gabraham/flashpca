#' Standardise a genotype matrix X.
#'
#' @param X A genotype matrix in dosage encoding 0, 1, 2.
#' @param type A character string indicating which type of standard deviation
#' to use.
#' @param impute Logical. Whether to impute missing values to zero (=mean
#' after standardisation).
#' @details{
#' type 1 (old Eigensoft way, Price 2006)
#'    mean = 2p (p is MAF)
#'    sd = sqrt(p (1 - p))
#'
#' type 2: (new Eigensoft way)
#'    mean = 2p (p is MAF)
#'    sd = sqrt(2 * p (1 - p))
#'
#' Missing genotypes ('NA') are imputed to the mean (=0).
#' }
#'
#' @export
scale2 <- function(X, type=c("2", "1"), impute=TRUE)
{
   type <- match.arg(type)
   mult <- ifelse(type == "1", 1, 2)

   sum2 <- nrow(X) - colSums(apply(X, 2, is.na))
   p <- colSums(X, na.rm=TRUE) / (2 * sum2)
   xsd <- sqrt(mult * p * (1 - p))
   names(p) <- names(xsd) <- colnames(X)

   s <- sweep(
      sweep(X, MARGIN=2, STATS=2 * p, FUN="-"),
         MARGIN=2, STATS=xsd, FUN="/"
   )
   if(impute) {
      s[is.na(s)] <- 0
   }
   attr(s, "scaled:center") <- 2 * p
   attr(s, "scaled:scale") <- xsd
   s
}

