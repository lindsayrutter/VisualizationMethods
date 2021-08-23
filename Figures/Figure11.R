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

# Number of clusters
nC = 4
# Number of samples
nCol = 6
colList = scales::hue_pal()(nC+1)
colList <- colList[c(4, 3, 2, 5, 1)]
# Hierarchical clustering
k = cutree(hc, k=nC)

yMin = min(datas[,1:nCol])
yMax = max(datas[,1:nCol])

# Create background boxplot data
boxDat <- melt(fulls, id.vars="ID")
colnames(boxDat) <- c("ID", "Sample", "Count")

logSoy = soybean_ir
logSoy[,-1] <- log(soybean_ir[,-1]+1)

# Plot scatterplot matrix for Cluster 1 significant genes
plotClusterSM(1)
