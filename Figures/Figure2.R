library(rtracklayer)
library(Rsamtools)
library(grid)
library(GenomicAlignments)
library(ggplot2)
library(GGally)
library(edgeR)
library(stringr)
library(EDASeq)
library(dplyr)
library(matrixStats)
library(gridExtra)
library(reshape2)
library(scales)
library(bigPint)

source("../functions.R")

load("../data/soybean_ir.rda")
data <- soybean_ir
load("../data/soybean_ir_noFilt_metrics.rda")
metrics <- soybean_ir_noFilt_metrics[["N_P"]]

# Filter, normalize, and standardize the data so each gene has mean=0 and stdev=1
res <- filterStandardizeSB(data)
# Fitered data standardized
filts <- res[["filts"]]
# Non-filtered data standardized
datas <- res[["datas"]]
# Hierarchical clustering object
hc <- res[["hc"]]
# Full data standardized
fulls <- rbind(datas, filts)

fulls <- fulls[,c(7,1:6)]
metrics <- list(metrics)
names(metrics) <- "N_P"

nC=4
colList = scales::hue_pal()(nC+1)
colList <- colList[c(3, 2, 5, 1)]

metricUnList = metrics[["N_P"]]
geneList = metricUnList[which(metricUnList$FDR < 0.05), ]$ID
fullGL = fulls[fulls$ID %in% geneList,]
ret <- plotClusters(data=datas[,c(7,1:6)], geneList = geneList, clusterAllData = TRUE, yAxisLabel = "Standardized count", colList = colList, saveFile = FALSE, vxAxis = TRUE, lineAlpha = 1, lineSize = 0.3)
plot(ret[["N_P_4"]])
