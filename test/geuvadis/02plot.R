
library(ggplot2)
library(data.table)
library(pROC)

load("geuvadis_data.RData")
load("geuvadis_pcca_model.RData")
load("geuvadis_fcca_model.RData")
load("geuvadis_rcca_model.RData")
load("geuvadis_scca_model.RData")
load("geuvadis_pma_model.RData")


df1 <- data.table(
   pop=sample_info$pop, method="FCCA", run.fcca$final_model$Py[,1:3])
df2 <- data.table(
   pop=sample_info$pop, method="PCCA", run.pcca$final_model$Py[,1:3])
df3 <- data.table(
   pop=sample_info$pop, method="RCCA", run.rcca$final_model$Py[,1:3])
df4 <- data.table(
   pop=sample_info$pop, method="SCCA", run.scca$final_model$Py[,1:3])
df5 <- data.table(
   pop=sample_info$pop, method="PMA", run.pma$Py[,1:3])
df.final <- rbind(df1, df2, df3, df4, df5)
df.final[, method_label :=
   factor(method, levels=c("PCCA", "RCCA", "SCCA", "FCCA", "PMA"))]

df.final2 <- df.final[method_label != "PMA",]

g1 <- ggplot(df.final2, aes(x=Py1, y=Py2, colour=pop))
g1 <- g1 + geom_point(alpha=0.5)
g1 <- g1 + facet_wrap(. ~ method_label, scales="free")
g1 <- g1 + theme_bw()
g1 <- g1 + scale_colour_viridis_d(name="Population")
g1 <- g1 + scale_x_continuous("Canonical coordinate 1")
g1 <- g1 + scale_y_continuous("Canonical coordinate 2")
g1 <- g1 + guides(colour=guide_legend(override.aes=list(alpha=1)))
ggsave(g1, file="geuvadis_all_models_final_by_pop.pdf", width=7, height=6)

df.final3 <- df.final[method_label %in% c("PCCA", "SCCA", "FCCA"),]
ggsave(g1 %+% df.final3,
   file="geuvadis_all_models_final_by_pop_v2.pdf", width=8, height=3)

# Cross-validated predictions
df1 <- data.table(
   pop=sample_info$pop, method="FCCA", run.fcca$final_model_cv_Py[,1:3])
df2 <- data.table(
   pop=sample_info$pop, method="PCCA", run.pcca$final_model_cv_Py[,1:3])
df3 <- data.table(
   pop=sample_info$pop, method="RCCA", run.rcca$final_model_cv_Py[,1:3])
df4 <- data.table(
   pop=sample_info$pop, method="SCCA", run.scca$final_model_cv_Py[,1:3])
df.cv <- rbind(df1, df2, df3, df4)
df.cv[, method_label :=
   factor(method, levels=c("PCCA", "RCCA", "SCCA", "FCCA"))]

# The coordinates from PCCA/SCCA/FCCA may be anti-correlaetd (even though may
# have been flipped internally)

g1 <- ggplot(df.cv, aes(x=Py1, y=Py2, colour=pop))
g1 <- g1 + geom_point(alpha=0.5)
g1 <- g1 + facet_wrap(. ~ method_label, scales="free")
g1 <- g1 + theme_bw()
g1 <- g1 + scale_colour_viridis_d(name="Population")
g1 <- g1 + scale_x_continuous("Canonical coordinate 1")
g1 <- g1 + scale_y_continuous("Canonical coordinate 2")
g1 <- g1 + guides(colour=guide_legend(override.aes=list(alpha=1)))
ggsave(g1, file="geuvadis_all_models_cv_by_pop.pdf", width=7, height=6)

df.cv2 <- df.cv[method_label %in% c("PCCA", "SCCA", "FCCA"), ]
ggsave(g1 %+% df.cv2, file="geuvadis_all_models_cv_by_pop_v2.pdf",
   width=8, height=3)

# Measure how well FIN is separated from the rest in Py2, for each method
df.final[, auc(factor(pop == "FIN") ~ Py2, quiet=TRUE), by=method]
df.final[, auc(factor(pop == "FIN") ~ Py3, quiet=TRUE), by=method]
df.cv[, auc(factor(pop == "FIN") ~ Py2, quiet=TRUE), by=method]
df.cv[, auc(factor(pop == "FIN") ~ Py3, quiet=TRUE), by=method]

