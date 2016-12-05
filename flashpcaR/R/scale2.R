
# Standardise genotype matrix x using 
# either
# type 1 (old Eigensoft way, Price 2006)
#    mean = 2p (p is MAF)
#    sd = sqrt(p (1 - p))
#
# type 2: (new Eigensoft way)
#    mean = 2p (p is MAF)
#    sd = sqrt(2 * p (1 - p))
#
# Missing genotypes ('NA') are imputed to the mean (=0)
#
scale2 <- function(x, type=c("2", "1"))
{
   type <- match.arg(type)
   mult <- ifelse(type == "1", 1, 2)

   sum2 <- nrow(x) - colSums(apply(x, 2, is.na))
   p <- colSums(x, na.rm=TRUE) / (2 * sum2)
   xsd <- sqrt(mult * p * (1 - p))
   names(p) <- names(xsd) <- colnames(x)

   s <- sweep(
      sweep(x, MARGIN=2, STATS=2 * p, FUN="-"),
         MARGIN=2, STATS=xsd, FUN="/"
   )
   s[is.na(s)] <- 0
   attr(s, "scaled:center") <- 2 * p
   attr(s, "scaled:scale") <- xsd
   s
}

