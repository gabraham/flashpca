
rm(list=ls())

library(data.table)

devtools::load_all("~/Code/flashpca/flashpcaR")

dir.expr <- "/storage/abyss/abrahamg-u/Data/ArrayExpress/E-TABM-1140"
dir.meth <- "/storage/abyss/abrahamg-u/Data/ArrayExpress/E-MTAB-1866"

d.expr <- fread(paste0(dir.expr,
   "/MuTHER_Fat_normalized_31032010_uncondensed_Ids.txt.gz"), skip=2)
d.expr.col <- fread(paste0(dir.expr,
   "/MuTHER_Fat_normalized_31032010_uncondensed_Ids.txt.gz"), nrows=2)

# The column names indicate technical replicates
d.expr.id <- sapply(strsplit(colnames(d.expr.col), "_"), head, n=1)
table(table(d.expr.id))
setnames(d.expr, d.expr.id)

d.meth <- fread(paste0(dir.meth,
   "/MuTHER_Fat_450K_norm_AE_030913.txt.gz"), skip=2)
d.meth.col <- fread(paste0(dir.meth,
   "/MuTHER_Fat_450K_norm_AE_030913.txt.gz"), nrows=2)
table(table(colnames(d.meth.col)))
setnames(d.meth, colnames(d.meth.col)) # hack, weird header format

# Illumina 450k array marker information
d.meth.info.all <- fread(paste0(dir.meth, "/",
   "HumanMethylation450_15017482_v1-2.csv.gz"), skip=7, fill=TRUE)

exrprb.id <- d.expr[, `Hybridization REF`]
cpg.id <- d.meth[, `Hybridization REF`]

d.expr.sinfo <- fread(paste0(dir.expr,
   "/E-TABM-1140.sdrf.txt"))
d.meth.sinfo <- fread(paste0(dir.meth,
   "/E-MTAB-1866.sdrf.txt"))

d.meth.id <- colnames(d.meth)
indiv <- intersect(d.expr.id, d.meth.id)
indiv <- indiv[indiv != "Hybridization REF"]
table(table(indiv))

d.sinfo <- d.meth.sinfo[,
   .(`Factor Value[individual]`,
     `Characteristics[co-twin]`,
     `Characteristics[twin zygosity]`,
     `Characteristics[sex]`)]
setnames(d.sinfo, c("indiv.id", "co.twin", "twin.zygosity", "sex"))
d.sinfo <- d.sinfo[data.table(indiv.id=indiv), on="indiv.id", ]

expr.idx <- match(indiv, d.expr.id)
meth.idx <- match(indiv, d.meth.id)

expr <- t(as.matrix(d.expr[, ..expr.idx]))
colnames(expr) <- exrprb.id

# These are B-values
meth <- t(as.matrix(d.meth[, ..meth.idx]))
colnames(meth) <- cpg.id

d.meth.info <- d.meth.info.all[
   data.table(IlmnID=colnames(meth)), on="IlmnID", nomatch=0]

# Remove non-autosomal CpGs
d.meth.info <- d.meth.info[
   Chromosome_36 %in% as.character(1:22), ]
d.meth.info[, Chromosome_36 := as.integer(Chromosome_36)]
d.meth.info[, Coordinate_36 := as.integer(Coordinate_36)]
setorder(d.meth.info, Chromosome_36, Coordinate_36)

# Exclude CpGs overlapping with SNPs
d.meth.snp <- fread(paste0(dir.meth, "/",
   "humanmethylation450_dbsnp137.snpupdate.table.v2.sorted.txt.gz"))
d.meth.info[d.meth.snp, on=c("IlmnID"="TargetID"), snp_overlap := i.SNP]

d.meth.info[, table(!is.na(snp_overlap))]

# Beta values
meth <- meth[, d.meth.info[is.na(snp_overlap), IlmnID]]

# missingness by sample 
summary(rowMeans(is.na(meth)))

# missingness by marker
summary(colMeans(is.na(meth)))

meth.b.imp <- meth
# For now, impute missing values to mean
for(j in 1:ncol(meth)) {
   m <- mean(meth[, j], na.rm=TRUE)
   meth.b.imp[is.na(meth[,j]), j] <- m
}

sum(is.na(meth.b.imp))

# M-values (ignoring beta offset)
meth.m.imp <- log2(meth.b.imp / (1 - meth.b.imp))

stopifnot(mode(expr) == "numeric")
stopifnot(mode(meth) == "numeric")
stopifnot(nrow(expr) == nrow(meth))

save(expr, file="twinsuk_fat_expr.RData")
save(meth.b.imp, meth.m.imp,
   file="twinsuk_fat_meth.RData")
save(d.meth.info, file="twinsuk_fat_meth_info.RData")

# Sample information
fwrite(d.sinfo, file="d.sinfo.txt.gz", sep="\t")


