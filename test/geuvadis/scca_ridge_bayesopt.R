
rm(list=ls())
graphics.off()

options(error=dump.frames)

#library(flashpcaR)
library(mixOmics)
library(PMA)
library(data.table)
library(doMC)
library(mlrMBO)
library(ggplot2)
library(boot)
library(biomaRt)
library(GGally)

library(devtools)

load_all("~/Code/flashpca/flashpcaR")

registerDoMC(cores=20)

set.seed(2721)

load("exp_snps.rda")

# These SNPs have been LD pruned already, and filtered by MAF
maf <- colMeans(snps) / 2
summary(maf)

X <- scale(snps[, sample(ncol(snps), 5000)])
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
colnames(r4$Px) <- paste0("Px", 1:ncol(r4$Px))

################################################################################
# FCCA, cross-validation
nfolds <- 5
folds <- sample(nfolds, nrow(X), replace=TRUE)

ndim <- 3

obj.fun <- smoof::makeSingleObjectiveFunction(
   name="scca",
   fn=function(x) {
      r <- cv.scca.ridge(X, Y, ndim=ndim, lambda1=x[1], lambda2=x[2],
	 gamma1=x[3], gamma2=x[4], folds=folds, verbose=FALSE)
      m <- r$result$avg.sq.cor
      cat("x:", x, "m:", m, "\n")
      ifelse(is.nan(m) || !is.finite(m), runif(1, 0, 1e-2), m)
   },
   par.set=makeParamSet(
      makeNumericParam("x1", lower=1e-4, upper=1e-1),
      makeNumericParam("x2", lower=1e-4, upper=1e-1),
      makeNumericParam("x3", lower=-3, upper=4, trafo=function(x) 10^x),
      makeNumericParam("x4", lower=-3, upper=4, trafo=function(x) 10^x)
   ),
   minimize=FALSE,
   noisy=TRUE
)

# Setup initial 'warmup' results for mlrMBO
run.time1 <- system.time({
   des.cv <- cv.scca.ridge(X, Y, ndim=ndim,
      lambda1=seq(1e-4, 0.002, length=3),
      lambda2=seq(1e-4, 0.002, length=3),
      gamma1=10^(-2:0), gamma2=10^(-2:0), folds=folds)
})
des <- copy(des.cv$result)
setcolorder(des, c(3, 4, 1, 2))
des[, gamma1 := log10(gamma1)]
des[, gamma2 := log10(gamma2)]
setnames(des, c("x1", "x2", "x3", "x4", "y"))
des <- as.data.frame(des)

ctrl <- makeMBOControl()
ctrl <- setMBOControlTermination(ctrl, iters=50L)
ctrl <- setMBOControlInfill(ctrl, crit=makeMBOInfillCritEI())
surr.km <- makeMBOLearner(control=ctrl, obj.fun)

run.time2 <- system.time({
   run <- mbo(obj.fun, design=des, learner=surr.km,
      control=ctrl, show.info=TRUE)
})
run.path <- as.data.table(run$opt.path)
setnames(run.path, c("x1", "x2", "x3", "x4"),
   c("lambda1", "lambda2", "log10.gamma1", "log10.gamma2"))

save(run, file="geuvadis_mbo_run.RData")

# sanity check
s1 <- cv.scca.ridge(X, Y, ndim=ndim,
  lambda1=run$x$x1, lambda2=run$x$x2,
  gamma1=10^(run$x$x3), gamma2=10^(run$x$x4),
  folds=folds, verbose=FALSE)
c(s1$result$avg.sq.cor, run$y)

pdf("scca_ridge_bayesopt_runs.pdf")
plot(run)
dev.off()

system.time({
   s3 <- scca.ridge(X, Y, ndim=ndim,
      lambda1=run$x$x1, lambda2=run$x$x2,
   gamma1=10^(run$x$x3), gamma2=10^(run$x$x4), verbose=TRUE)
})

# sanity check, should be faster by initialising V
system.time({
   s4 <- scca.ridge(X, Y, ndim=ndim,
      lambda1=run$x$x1, lambda2=run$x$x2,
      gamma1=10^(run$x$x3), gamma2=10^(run$x$x4), V=s3$V, verbose=TRUE)
})

mean(abs(s3$a - s4$a))
mean(abs(s3$b - s4$b))

################################################################################
# Bootstrap results

boot.fun <- function(original, idx)
{
   s <- scca.ridge(X[idx,], Y[idx,], ndim=ndim,
      lambda1=run$x$x1, lambda2=run$x$x2,
      gamma1=10^(run$x$x3), gamma2=10^(run$x$x4),
      V=s3$V)

   # boot() can't handle returning of matrices, so we squeeze the matrices
   # as a single vector. Avoids having to re-run bootstrap for each
   # of a and b and their dimensions.
   # Note that as.numeric operates _by_column_
   #c(as.numeric(s$a), as.numeric(s$b))
   as.numeric(s$b)
}

# If this is too big, the results will be VERY large
nreps <- 1000
res.boot <- boot(data=data.frame(idx=1:nrow(X)),
   statistic=boot.fun, R=nreps, parallel="multicore", ncpus=30)

ncolx <- ncol(X)
ncoly <- ncol(Y)
#var.types <- rep(c("a", "b"), c(ncolx * ndim, ncoly * ndim))
#var.names <- rep(c(colnames(X), colnames(Y)), each=ndim)
var.types <- rep("b", ncoly * ndim)
var.names <- rep(colnames(Y), ndim)

res.boot.ci <- foreach(j=seq(along=res.boot$t0), .combine="rbind") %dopar% {
   bc <- try(boot.ci(res.boot, index=j, conf=0.95, type=c("norm", "perc")))
   #if(var.types[j] == "a") {
   #   dim <- floor((j - 1) / ncolx) + 1
   #} else {
   #   dim <- floor((j - ncolx * ndim - 1) / ncoly) + 1
   #}
   dim <- floor((j - 1) / ncoly) + 1
   if(is(bc, "try-error")) {
      r <- data.table(idx=j, var.type=var.types[j],
	 var.name=var.names[j], dim=dim,
	 boot.method=NA, t0=res.boot$t0[j],
	 lower=NA, upper=NA)
      return(r)
   }
   m <- t(sapply(bc[-(1:3)], function(x) tail(c(x), n=2)))
   colnames(m) <- c("lower", "upper")
   r <- data.table(idx=j, var.type=var.types[j],
      var.name=var.names[j], dim=dim,
      boot.method=rownames(m), t0=res.boot$t0[j], m)
   r
}

# Good explanation about the Normal bootstrap
# at https://stats.stackexchange.com/questions/20701/computing-p-value-using-bootstrap-with-r
#
# The below is the same as boot::norm.ci, with the addition of a two-sided p-value.

# check whether bootstrap samples are approximately normal
png("geuvadis_expr_bootstrap_checks.png", width=900, height=900)
par(mfrow=c(5, 5))
for(j in sample(ncol(res.boot$t), 25)) {
   qqnorm(res.boot$t[,j])
   qqline(res.boot$t[,j])
}
dev.off()

res.boot.mean <- colMeans(res.boot$t)
res.boot.se <- sqrt(apply(res.boot$t, 2, var))
res.boot.bias <- res.boot.mean - res.boot$t0
res.boot.ci.norm <- cbind(
   res.boot$t0 - res.boot.bias + qnorm(0.025) * res.boot.se,
   res.boot$t0 - res.boot.bias + qnorm(0.975) * res.boot.se
)
res.boot.ci.norm.p <- 2 * pnorm(
   abs((2 * res.boot$t0 - res.boot.mean) / res.boot.se), lower=FALSE)
res.boot.ci.norm.p.adj <- p.adjust(res.boot.ci.norm.p, "fdr")
names(res.boot.ci.norm.p.adj) <- colnames(Y)
res.boot.ci[boot.method == "normal",
   c("mean", "p.value", "res.boot.ci.norm.p.adj") := list(
      res.boot.mean, res.boot.ci.norm.p, res.boot.ci.norm.p.adj)]

res.boot.b.signif <- res.boot.ci[boot.method == "normal" & p.value.adj < 0.05,]
res.boot.b.signif[, gene.ens := gsub("\\.[0-9]*$", "", var.name)]

mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)

bm.res <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "gene_biotype", "external_gene_name"),
  filter="ensembl_gene_id",
  values=res.boot.b.signif$gene.ens,
  uniqueRows=TRUE)
setDT(bm.res)
res.boot.b.signif[bm.res,
   on=c("gene.ens"="ensembl_gene_id"),
   symbol := i.external_gene_name]
res.boot.b.signif[, t0_abs := abs(t0)]
setorder(res.boot.b.signif, -t0_abs)

top.k <- 5
res.boot.b.signif.top <- res.boot.b.signif[,
   head(.SD[, .(var.name, symbol, t0)], n=top.k), by=dim]
exp_gene_top <- data.table(pop=sample_info$pop,
   exp_gene[, res.boot.b.signif.top$var.name])
setnames(exp_gene_top, c("pop", res.boot.b.signif.top$symbol))
exp_gene_top.m <- melt(exp_gene_top, id.vars="pop")

# check whether a gene appears in more than one dimension
table(colSums(res.boot.b.signif.top[, table(dim, symbol)]))

exp_gene_top.m[res.boot.b.signif.top, on=c("variable"="symbol"),
   dim := i.dim]
exp_gene_top.m[, dim_label := paste0("Dimension ", dim)]

g2 <- ggplot(exp_gene_top.m[dim %in% 1:2,], aes(x=pop, y=value, colour=pop))
g2 <- g2 + geom_violin()
g2 <- g2 + geom_point(alpha=0.1)
g2 <- g2 + facet_wrap(dim_label ~ variable, scales="free_y", ncol=top.k)
g2 <- g2 + theme_bw()
g2 <- g2 + scale_x_discrete("Population")
g2 <- g2 + scale_y_continuous("Gene expression level")
g2 <- g2 + scale_colour_viridis_d(name="Population")
g2 <- g2 + theme(legend.position="none")

ggsave(g2, file="geuvadis_scca_ridge_bayesopt_top_genes_by_pop.png", width=8)

################################################################################


pdf("scca_ridge_bayesopt_final_by_pop.pdf", width=8, height=8)
pairs(cbind(s3$Px[,1:3], s3$Py[, 1:3]),
   col=factor(sample_info$pop), main="SCCA bayes optim", gap=0)
pairs(cbind(r4$Px[,1:3], r4$Py[, 1:3]),
   col=factor(sample_info$pop), main="PCCA", gap=0)
dev.off()

df1 <- data.table(
   pop=sample_info$pop, method="FCCA", s3$Py[,1:3])
df2 <- data.table(
   pop=sample_info$pop, method="PCCA", r4$Py[,1:3])
df <- rbind(df1, df2)
df[, method_label := factor(method, levels=c("PCCA", "FCCA"))]


g1 <- ggplot(df, aes(x=Py1, y=Py2, colour=pop))
g1 <- g1 + geom_point(alpha=0.5)
g1 <- g1 + facet_wrap(. ~ method_label, scales="free")
g1 <- g1 + theme_bw()
g1 <- g1 + scale_colour_viridis_d(name="Population")
g1 <- g1 + scale_x_continuous("Canonical coordinate 1")
g1 <- g1 + scale_y_continuous("Canonical coordinate 2")
g1 <- g1 + guides(colour=guide_legend(override.aes=list(alpha=1)))

ggsave(g1, file="scca_ridge_bayesopt_final_by_pop_v2.pdf",
   width=7, height=3.5)

g2 <- ggscatmat(df[method == "FCCA",], columns=3:5, color="pop")
ggsave(g2, file="scca_ridge_bayesopt_final_by_pop_alldims_v2.pdf")

warning("test-set R^2 are not decreasing with dimension, necessarily")
warning("do we need to sort the canonical directions somehow?")

save.image(file="geuvadis_final_results.RData")

