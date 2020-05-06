library(VIM)
library(missForest)
library(data.table)

###############################################################################
# Based on https://github.com/selbouhaddani/OmicsPLS/blob/master/vignettes/OmicsPLS_vignette.pdf
#

# Download the gene expression data from ArrayExpress, if hasn't been already
f <- "~/Downloads/E-TABM-1036.processed.1.zip"
if(!file.exists(f)) {
   download.file(
      paste0("http://www.ebi.ac.uk/arrayexpress/files/E-TABM-1036/",
         basename(f)), f)
}
rna0 <- data.table::fread(unzip(f, "test.tab"))
rna1 <- t(as.data.frame.matrix(rna0[-1, -1, with=F]))
rna2 <- matrix(as.numeric(rna1), nrow=nrow(rna1))
dimnames(rna2) <- list(colnames(rna0)[-1],unlist(rna0[-1, 1, with=FALSE]))

# There shouldn't be any missingness in the gene expression data
stopifnot(sum(is.na(rna2)) == 0)

# Order rows according to the participant ID
rna2 <- rna2[order(row.names(rna2)), ] 

filter_rna <- function(rna=rna, prop=0.75)
{
   # calculate the maximum of gene expression for each gene and take the top
   maxGE <- apply(rna, 2, max)
   propGEmax <- quantile(maxGE, prop)
   
   # take the IQR of each gene and take the top genes
   IQRGE <- apply(rna, 2, IQR, na.rm=TRUE)
   propGEIQR <- quantile(IQRGE, prop)
   
   # selected genes/probes are the intersection of the two previous sets
   intersect(which(maxGE > propGEmax), which(IQRGE > propGEIQR))
}
rna3 <- rna2[,filter_rna(rna2, prop=0.5)]

# Download the metabolite data, if hasn't been already
f <- "~/Downloads/msb201093-sup-0002.zip"
if(!file.exists(f)) {
   download.file(
      "https://www.embopress.org/action/downloadSupplement?doi=10.1038%2Fmsb.2010.93&file=msb201093-sup-0002.zip",
      f)
}
metab0 <- read.table(unzip(f, "metabonomic_data.txt"), header=TRUE) 
metab1 <- t(metab0[, -1])
colnames(metab1) <- metab0$Metabolite

# Visualise the missingness in the metabolites
try(VIM::aggr(metab1, col=c("navyblue", "red"), numbers=TRUE,
   sortVars=FALSE, labels=names(data), cex.axis=0.7, gap=3,
   ylab=c("Histogram of missing data", "Pattern")))

# Remove individuals with 100% missing metabolites
NAs_in_metab1 <- which(rowMeans(is.na(metab1)) == 1.0)
table(NAs_in_metab1)
metab2 <- metab1[-NAs_in_metab1, ]
rna4 <- rna3[-NAs_in_metab1, ]

# Impute the missing metabolites
metab2.imp <- missForest::missForest(metab2, verbose=TRUE)

# Center the variables
metab <- scale(metab2.imp$ximp, scale=FALSE)
rna <- scale(rna4, scale=FALSE)


# Look at the structure of the metabolites before imputation
try(
   gplots::heatmap.2(cor(metab1, use="pair"), Rowv=FALSE, Colv=FALSE,
      trace="n", breaks=seq(-1, 1, length.out=25), col=gplots::bluered),
   silent=TRUE)

# Look at the structure of the metabolites after imputation
try(
   gplots::heatmap.2(cor(metab, use="pair"), Rowv=FALSE, Colv=FALSE,
      trace="n", breaks=seq(-1, 1, length.out=25), col=gplots::bluered),
   silent=TRUE)

# Look at distribution of the top 95% of varying genes
boxplot(rna[, filter_rna(rna, 0.95)])

# Look at distributions of metabolites
boxplot(metab)

fwrite(round(rna, 6), file="rna.txt.gz", sep="\t")
fwrite(round(metab, 6), file="metab.txt.gz", sep="\t")



