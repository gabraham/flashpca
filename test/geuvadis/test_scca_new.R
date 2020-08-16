
rm(list=ls())
graphics.off()

options(error=dump.frames)

#library(flashpcaR)
library(mixOmics)
library(PMA)
library(data.table)
library(doMC)

library(devtools)

load_all("~/Code/flashpca/flashpcaR")

registerDoMC(cores=10)

set.seed(27211)

#load("exp_gene_sample.rda")
#load("snps_v.rda")
load("exp_snps.rda")

maf <- colMeans(snps) / 2
snps <- snps[, maf > 0.05]
s_var <- apply(snps, 2, var)
snps_v <- snps[, s_var > 0.54]

#X <- scale(snps_v[, sample(ncol(snps_v), 1000)])
#Y <- scale(exp_gene[, sample(ncol(exp_gene), 5000)])
X <- scale(snps_v[, sample(ncol(snps_v), 5000)])
Y <- scale(exp_gene)


################################################################################
# PCCA

n <- nrow(X)
f1 <- svd(X)
f2 <- svd(Y)
mx <- f1$d > 1e-10
my <- f2$d > 1e-10
kx <- 5 # number of SNP PCs
ky <- 30 # number of gene expression PCs

summary(colMeans(f1$u))
summary(colMeans(f2$u))
r4 <- cancor(f1$u[, 1:kx], f2$u[, 1:ky], xcenter=FALSE, ycenter=FALSE)
r4$Px <- f1$u[, 1:kx] %*% r4$xcoef
r4$Py <- f2$u[, 1:ky] %*% r4$ycoef
colnames(r4$Py) <- paste0("Py", 1:ncol(r4$Py))

################################################################################

lambda1 <- seq(1e-4, 0.002, length=15)
lambda2 <- seq(1e-4, 0.005, length=10)
gamma1 <- 10^(-4:4)
gamma2 <- 10^(-4:4)

nfolds <- 10
folds <- sample(nfolds, nrow(X), replace=TRUE)

n <- nrow(X) 
f1 <- svd(X) # note: svd of X, not X/sqrt(n-1)
f2 <- svd(Y)
mx <- f1$d > 1e-10
my <- f2$d > 1e-10

g1w <- g2w <- 0.1

# Whitened data (ridge penalty)
Xw <- with(f1,
   tcrossprod(u[,mx] %*% diag(d[mx] / sqrt(d[mx]^2 + (n - 1) * g1w)), v[,mx]))
Yw <- with(f2,
   tcrossprod(u[,my] %*% diag(d[my] / sqrt(d[my]^2 + (n - 1) * g2w)), v[,my]))
sx.invsqrt <- with(f1,
   tcrossprod(
      v[,mx] %*% diag(sqrt(n - 1) / sqrt(d[mx]^2 + (n - 1) * g1w)), v[,mx]))
sy.invsqrt <- with(f2,
   tcrossprod(
      v[,my] %*% diag(sqrt(n - 1) / sqrt(d[my]^2 + (n - 1) * g2w)), v[,my]))

# SCCA on whitened data
# divided internally by sqrt(n - 1), but we don't want that
Xw1 <- Xw * sqrt(n - 1)
Yw1 <- Yw * sqrt(n - 1)
s3 <- cv.scca(Xw1, Yw1, folds=folds, lambda1=lambda1, lambda2=lambda2,
   parallel=TRUE, standx="none", standy="none", ndim=1, verbose=FALSE)
s4 <- scca(Xw1, Yw1, lambda1=s3$best.lambda1, lambda2=s3$best.lambda2,
   standx="none", standy="none", ndim=3, verbose=FALSE)

s5 <- scca.ridge(X, Y, gamma1=g1w, gamma2=g2w,
   lambda1=s3$best.lambda1, lambda2=s3$best.lambda2,
   ndim=5, verbose=FALSE)

a4 <- sx.invsqrt %*% s4$U
rownames(a4) <- colnames(X)
Px4 <- X %*% a4

b4 <- sy.invsqrt %*% s4$V
rownames(b4) <- colnames(Y)
Py4 <- Y %*% b4


################################################################################

system.time({
   s6 <- cv.scca.ridge(X, Y, ndim=3, gamma1=gamma1, gamma2=gamma2,
      lambda1=lambda1, lambda2=lambda2, folds=folds, verbose=FALSE)
})
r6m <- s6$result[, list(r.trn.mean=mean(r.trn), r.trn.se=sqrt(var(r.trn) / .N),
   r.tst.mean=mean(r.tst), r.tst.se=sqrt(var(r.tst) / .N)),
   by=.(dim, gamma1, gamma2, lambda1, lambda2)]
r6m[, gamma1_label := paste0("gamma1:", gamma1)]
r6m[, gamma2_label := paste0("gamma2:", gamma2)]

# Sum of squared correlations across the dimensions, for each model
# on the penalty grid. I.e., the best overall model taking all dimensions into account.
r6m.sum <- r6m[, list(avg.sq.cor=mean(r.tst.mean^2, na.rm=TRUE)),
   by=.(gamma1, gamma2, lambda1, lambda2)]
r6m.sum.best <- r6m.sum[which.max(avg.sq.cor), ]

# Which is best for each dimension (may be different for each one)
r6m.best <- r6m[, .SD[which.max(r.tst.mean),], by=dim]
r6m.best.d1 <- r6m.best[dim == 1, ]

g1 <- ggplot(r6m.sum, aes(x=lambda1, y=lambda2, fill=avg.sq.cor))
g1 <- g1 + geom_raster()
g1 <- g1 + facet_grid(gamma1 ~ gamma2)
g1 <- g1 + theme_bw()
g1 <- g1 + scale_fill_viridis_c("Average squared\ncorrelation")
g1 <- g1 + scale_x_continuous(guide=guide_axis(check.overlap=TRUE))
g1 <- g1 + scale_y_continuous(guide=guide_axis(check.overlap=TRUE))
g1 <- g1 + geom_point(data=r6m.sum.best, shape="+", size=15)
#g1 <- g1 + geom_point(data=r6m.best.d1,
#	mapping=aes(x=lambda1, y=lambda2), fill="black",
#	shape="x", size=15)

ggsave(g1, file="geuvadis_scca_ridge_avgcorrs_cv.png", width=12, height=11.5)


#r6m.best[, r.tst.1se := r.tst.mean - r.tst.se]

# Select all the models with test performance within 1 stderr
# of the best model
#r6m.1se <- r6m[r.tst.mean >= r6m.best$r.tst.1se, ]

# Out of all canadidate models within 1 stderr, select the most consistent
# model, i.e., least difference between train and test performance
#r6m.1se.best <- r6m.1se[which.min((r.trn.mean - r.tst.mean)^2), ]

s7 <- scca.ridge(X, Y, gamma1=r6m.best.d1$gamma1, gamma2=r6m.best.d1$gamma2,
   lambda1=r6m.best.d1$lambda1, lambda2=r6m.best.d1$lambda2, ndim=5)
s8 <- scca.ridge(X, Y, gamma1=r6m.sum.best$gamma1, gamma2=r6m.sum.best$gamma2,
   lambda1=r6m.sum.best$lambda1, lambda2=r6m.sum.best$lambda2, ndim=5)

################################################################################
# manual cross-validation to test if one model really better than another
res.cv <- lapply(1:nfolds, function(fold) {
   Xtrn <- X[folds != fold,]
   Xtst <- X[folds == fold,]
   Ytrn <- Y[folds != fold,]
   Ytst <- Y[folds == fold,]
   sf1 <- scca.ridge(Xtrn, Ytrn, gamma1=g1w, gamma2=g1w,
      lambda1=s3$best.lambda1, lambda2=s3$best.lambda2, ndim=3)
   sf2 <- scca.ridge(Xtrn, Ytrn,
      gamma1=r6m.sum.best$gamma1, gamma2=r6m.sum.best$gamma2,
      lambda1=r6m.sum.best$lambda1, lambda2=r6m.sum.best$lambda2, ndim=3)
   sf3 <- scca.ridge(Xtrn, Ytrn,
      gamma1=r6m.best$gamma1[1], gamma2=r6m.best$gamma2[1],
      lambda1=r6m.best$lambda1[1], lambda2=r6m.best$lambda2[1], ndim=3)
   Px1 <- Xtst %*% sf1$a
   Py1 <- Ytst %*% sf1$b
   Px2 <- Xtst %*% sf2$a
   Py2 <- Ytst %*% sf2$b
   Px3 <- Xtst %*% sf3$a
   Py3 <- Ytst %*% sf3$b
   rbind(
      data.table(model=1, fold=fold, rbind(diag(cor(Px1, Py1)))),
      data.table(model=2, fold=fold, rbind(diag(cor(Px2, Py2)))),
      data.table(model=3, fold=fold, rbind(diag(cor(Px3, Py3))))
   )
})
res.cv <- do.call(rbind, res.cv)

# Average test correlation for each model, across the dimensions
res.cv.avg <- res.cv[, list(mean(V1), mean(V2), mean(V3)), by=model]

#res.cv.tst1 <- cv.scca.ridge(X, Y,
#   gamma1=g1w, gamma2=g2w,
#   lambda1=s3$best.lambda1, lambda2=s3$best.lambda2, folds=folds)
#res.cv.tst2 <- cv.scca.ridge(X, Y,
#   gamma1=r6m.1se.best$gamma1, gamma2=r6m.1se.best$gamma2,
#   lambda1=r6m.1se.best$lambda1, lambda2=r6m.1se.best$lambda2, folds=folds)
#res.cv.tst3 <- cv.scca.ridge(X, Y,
#   gamma1=r6m.best$gamma1, gamma2=r6m.best$gamma2,
#   lambda1=r6m.best$lambda1, lambda2=r6m.best$lambda2, folds=folds)
#
#mean(abs(res.cv.tst1$result[, r.tst] - res.cv[model == 1, V1]))
#mean(abs(res.cv.tst2$result[, r.tst] - res.cv[model == 2, V1]))
#mean(abs(res.cv.tst3$result[, r.tst] - res.cv[model == 3, V1]))
#
## fixed gamma model
#m1 <- r6m[gamma1 == g1w & gamma2 == g2w & lambda1 == s3$best.lambda1 & lambda2 == s3$best.lambda2,]
#m2 <- res.cv[model == 1, mean(V1)]
#abs(m1$r.tst.mean - m2)
#
## Test for significance of difference between the two models
## TODO: maybe Fisher z-transform/atahn for correlations
#z.diff <- (m1$r.tst.mean - r6m.best$r.tst.mean) / sqrt(m1$r.tst.se^2 +
#   r6m.best$r.tst.se^2)
#z.diff.pval <- pnorm(abs(z.diff), lower=FALSE) * 2
#
#r6m[gamma1 == 1 & gamma2 == 1e-4 & lambda1 == 1e-4 & (abs(lambda2 -
#0.0006444444) < 1e-7), ]

warning("need to check that the results from cv.scca.ridge match exactly what we get in the manual cv")
warning("and if they do match, is the best model selected")
warning("and if the best model is selected, is it meaningfully better")
warning("and should be use reconstruction error rather than correlation to select models")
warning("reconstruction error of what? Brown uses reconstruction of gene expression",
   "but that's not what we're trying to model, we're doing SVD of the X'Y correlation matrix")
warning("since optimising for correlation is only useful for 1st dimension")

pdf("geuvadis_test_scca.pdf", width=9, height=12)
par(mfrow=c(3, 2), mar=c(4, 3, 3, 1))
plot(r4$Py[,1:2], col=factor(sample_info$pop), main="PCCA")
legend("topright", lwd=2, col=1:4, legend=levels(factor(sample_info$pop)))
plot(s4$Py[,1:2], col=factor(sample_info$pop), main="scca on whitened, fixed gamma")
legend("topright", lwd=2, col=1:4, legend=levels(factor(sample_info$pop)))
plot(s5$Py[,1:2], col=factor(sample_info$pop), main="scca.ridge, fixed gamma" )
legend("topright", lwd=2, col=1:4, legend=levels(factor(sample_info$pop)))
plot(s7$Py[,1:2], col=factor(sample_info$pop), main="scca.ridge, best CV gamma d1")
legend("topright", lwd=2, col=1:4, legend=levels(factor(sample_info$pop)))
plot(s8$Py[,1:2], col=factor(sample_info$pop), main="scca.ridge, best CV gamma sum")
legend("topright", lwd=2, col=1:4, legend=levels(factor(sample_info$pop)))
dev.off()

pdf("geuvadis_test_scca_pairs.pdf")
pairs(s5$Py, gap=0, col=factor(sample_info$pop), main="scca.ridge, fixed gamma")
pairs(s8$Py, gap=0, col=factor(sample_info$pop), main="scca.ridge, CV gamma sum")
dev.off()

pdf("geuvadis_test_scca_xy_pairs.pdf")
par(mfrow=c(2, 2))
plot(s5$Px[,1], s5$Py[,1], col=factor(sample_info$pop), main="scca.ridge, fixed gamma")
plot(s5$Px[,2], s5$Py[,2], col=factor(sample_info$pop), main="scca.ridge, fixed gamma")
plot(s8$Px[,1], s8$Py[,1], col=factor(sample_info$pop), main="scca.ridge, CV gamma sum")
plot(s8$Px[,2], s8$Py[,2], col=factor(sample_info$pop), main="scca.ridge, CV gamma sum")
dev.off()

################################################################################
# Seems like optimising CCA for 1st dimension separates out a subgroup of YRI that
# is also visible in PCA
# So optimising for 1st dimension hurts separating out FINs, but helps separate out
# more of the YRI...
#y1 <- s8$Py[,1]
#y2 <- s8$Py[,2]
#grp1 <- names(y1)[y2 > -1 & y1 > 1.5]
#grp2 <- names(y1)[y2 < -1 & y1 > 1.5]
#nm <- c(grp1, grp2)
#grp <- nm %in% grp1
#f1 <- flashpca(X[nm,], stand="none", check_geno=FALSE)
#png("geuvadis_test_genotypes_yri_pca.png")
#pairs(f1$projection[, 1:5], gap=0, col=factor(grp))
#dev.off()

