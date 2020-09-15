
library(data.table)
library(ggplot2)

# The gene expression is from Illumina HT-12 chip
#url <- "https://sapac.support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/humanht-12/v3/humanht-12_v3_0_r3_11283641_a_txt.zip"
#tmpdir <- tempdir()
#f <- paste0(tmpdir, basename(url))
#if(!file.exists(f)) {
#   download.file(url, f)
#}
#annot <- fread(unzip(f, "HumanHT-12_V3_0_R3_11283641_A.txt"),
#   skip=8, fill=TRUE)
#save(annot, file="annot.rda")

load("annot.rda")

load("twinsuk_fat_meth_info.RData")
load("twinsuk_LL_data.RData")
load("twinsuk_LL_fcca_boot_ci.RData")
load("twinsuk_LL_fcca_boot.RData")

res.boot.ci.expr <- res.boot.ci[var.type == "a", ]
res.boot.ci.meth <- res.boot.ci[var.type == "b", ]

setnames(res.boot.ci.expr, "var.name", "Probe_Id")
res.boot.ci.expr[annot, on="Probe_Id", Symbol := i.Symbol]

res.boot.ci.meth[d.meth.info, on=c("var.name"="IlmnID"),
   UCSC_RefGene_Name := i.UCSC_RefGene_Name]   

# Select the markers with CI not overlapping with zero
res.boot.expr.signif <- res.boot.ci.expr[
   sign(lower) == sign(upper) & boot.method == "percent",]
res.boot.meth.signif <- res.boot.ci.meth[
   sign(lower) == sign(upper) & boot.method == "percent",]

res.boot.expr.signif[, t0_abs := abs(t0)]
res.boot.expr.signif[, t.boot.median_abs := abs(t.boot.median)]
setorder(res.boot.expr.signif, -t.boot.median_abs)

res.boot.meth.signif[, t0_abs := abs(t0)]
res.boot.meth.signif[, t.boot.median_abs := abs(t.boot.median)]
setorder(res.boot.meth.signif, -t.boot.median_abs)


stop()

################################################################################
# Plot the gene loadings
res.boot.ci.perc <- copy(res.boot.ci[boot.method == "percent",])
res.boot.ci.perc[, significant := sign(lower) == sign(upper)]

res.boot.ci.perc[an.mod, on="Probe_Id", AN.module := Replicated.Module]
res.boot.ci.perc[is.na(AN.module), AN.module := "(not classified)"]

# Plot the 'loadings'
res.boot.ci.perc.wide <- dcast(
   res.boot.ci.perc,
   AN.module + Probe_Id + Symbol ~ dim, value.var=c("t.boot.median", "significant"))
res.boot.ci.perc.wide[, significant_12 := significant_1 | significant_2]
res.boot.ci.perc.wide[, significant_22 := significant_2 | significant_3]
res.boot.ci.perc.wide.signif <- res.boot.ci.perc.wide[significant_12 == TRUE, ]
#res.boot.ci.perc.wide.signif[,
#   c("t.boot.median_1_abs", "t.boot.median_2_abs", "t.boot.median_3")

g1 <- ggplot(res.boot.ci.perc.wide.signif,
   aes(x=t.boot.median_1, y=t.boot.median_2, colour=AN.module))
g1 <- g1 + geom_hline(yintercept=0, linetype="dashed", alpha=0.5)
g1 <- g1 + geom_vline(xintercept=0, linetype="dashed", alpha=0.5)
g1 <- g1 + geom_point(alpha=0.7)
g1 <- g1 + theme_bw()

ggsave(g1, file="twinsuk_expr_bootstrap_weights.png", width=7, height=6)






