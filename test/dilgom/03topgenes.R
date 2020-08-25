
library(data.table)
library(ggplot2)

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

load("rna.RData")
load("metab.RData")

load("dilgom_scca_ridge_boot_ci.RData")
load("dilgom_scca_ridge_boot.RData")

setnames(res.boot.ci, "var.name", "Probe_Id")
res.boot.ci[annot, on="Probe_Id", Symbol := i.Symbol]


# Select the genes with CI not overlapping with zero
res.boot.a.signif <- res.boot.ci[sign(lower) == sign(upper)
    & boot.method == "percent",]

res.boot.a.signif[, t0_abs := abs(t0)]
res.boot.a.signif[, t.boot.median_abs := abs(t.boot.median)]

setorder(res.boot.a.signif, -t.boot.median_abs)


#pdf("dilgom_scca_ridge_bayesopt_final.pdf", width=8, height=8)
#pairs(cbind(s3$Px[,1:3], s3$Py[, 1:3]), main="SCCA bayes optim", gap=0)
#dev.off()
#
#dx1 <- data.table(Probe_Id=rownames(s3$a), s3$a)
#dy1 <- data.table(Metabolite=rownames(s3$b), s3$b)
#
#
#
#dx1[annot, on="Probe_Id", Symbol := i.Symbol]

LL.mod <- data.table(
   Symbol=c("CPA3", "ENPP3", "FCER1A", "GATA2",
      "HDC", "HS.132563", "MS4A2", "MS4A3", "MS4A3", "SLC45A3", "SPRYD5"),
   Probe_Id=c("ILMN_1766551", "ILMN_1749131", "ILMN_1688423",
   "ILMN_2102670", "ILMN_1792323", "ILMN_1899034", "ILMN_1806721",
   "ILMN_1695530", "ILMN_1751625", "ILMN_1726114", "ILMN_1753648"))

LL.mod[, in.data := Probe_Id %in% colnames(rna)]

res.boot.a.signif[, in.LL.module := FALSE]
res.boot.a.signif[LL.mod, on="Probe_Id", in.LL.module := TRUE]

#setnames(dx1, paste0("V", 1:3), paste0("Component ", 1:3))
#setnames(dy1, paste0("V", 1:3), paste0("Component ", 1:3))
#
#dx1[order(abs(`Component 1`), decreasing=TRUE)[1:20], ]
#dx1[order(abs(`Component 2`), decreasing=TRUE)[1:20], ]
#dx1[order(abs(`Component 3`), decreasing=TRUE)[1:20], ]
#
#dy1[order(abs(`Component 1`), decreasing=TRUE)[1:20], ]
#dy1[order(abs(`Component 2`), decreasing=TRUE)[1:20], ]
#dy1[order(abs(`Component 3`), decreasing=TRUE)[1:20], ]

an.mod <- fread("Module_Genes_DILGOM.csv", skip=2)[, -5]
setnames(an.mod, c("Probe_Id", "Symbol", "Replicated.Module",
   "CoreModuleGenesDILGOM.YFS"))

res.boot.b.signif[an.mod, on="Probe_Id", AN.module := Replicated.Module]

#dx1d <- melt(dx1, measure.vars=grep("^Component", colnames(dx1), value=TRUE))
#dx1d[, pos := 1:.N, by=variable]
#dx1d[, type := "gene"]
#dx1d[, name := Probe_Id]
#dx1d[, col := ifelse(in.LL.module, "red", "black"), by=variable]
#
#g1 <- ggplot(dx1d, aes(sample=value))
## Need to hack the colour so as to not split the plot into group
#g1 <- g1 + stat_qq(
#   colour=dx1d$col[dx1d[, order(value), by=variable]$V1]) + stat_qq_line()
#g1 <- g1 + facet_grid(. ~ variable)
#g1 <- g1 + theme_bw()
#g1 <- g1 + geom_hline(yintercept=0, linetype="dashed", alpha=0.5)
#ggsave(g1, file="dilgom_scca_ridge_expr.png", width=10, height=4)
#
#setnames(dy1, paste0("V", 1:3), paste0("Component ", 1:3))
#dy1d <- melt(dy1, measure.vars=grep("^Component", colnames(dy1), value=TRUE))
#g2 <- ggplot(dy1d, aes(sample=value))
## Need to hack the colour so as to not split the plot into group
#g2 <- g2 + stat_qq()
##   colour=dx1d$col[dx1d[, order(value), by=variable]$V1]) 
#g2 <- g2 + stat_qq_line()
#g2 <- g2 + facet_grid(. ~ variable)
#g2 <- g2 + theme_bw()
#g2 <- g2 + geom_hline(yintercept=0, linetype="dashed", alpha=0.5)
#ggsave(g2, file="dilgom_scca_ridge_metab.png", width=10, height=4)
#
#stop()
#
#Pxy <- cbind(r2$Px, r2$Py)
#
#warning("these correlations are way too high, overfitting?")
#g2 <- ggpairs(data.frame(Pxy), mapping=ggplot2::aes(alpha=0.5))
#g2 <- g2 + theme_bw()
#ggsave(g2, file="dilgom_cca_Pxy.png", width=9, height=9)
