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
library(data.table)

source("functions.R")

# Get data
load("data/LK_data.RData")
data = as.data.frame(MA.subsetA$M)
rownames(data) = as.character(MA.subsetA$genes$EnsemblGeneID)
setDT(data, keep.rownames = TRUE)[]
colnames(data) = c("ID","K.R1L1","L.R1L2","K.R1L3","L.R1L4","L.R1L6","K.R1L7","L.R1L8","K.R2L2","L.R2L3","K.R2L6")
data = as.data.frame(data)
data = data[,c(1,2,4,7,9,11,3,5,6,8,10)]
# Obtain R1 values
data <- data[,c(1:4,7:9)]
colnames(data) <- c("ID", "K.1", "K.2", "K.3", "L.1", "L.2", "L.3")

logSoy = data
logSoy[,-1] <- log(data[,-1]+1)

# Get raw data metrics
rowNames = data[,1]
dataSel = data[,-1]
rownames(dataSel) = rowNames
group <- factor(c(1,1,1,2,2,2))
y <- DGEList(counts=dataSel,group=group)
design <- model.matrix(~group)
y <- estimateDisp(y, design)
fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
ret = data.frame(ID=rownames(topTags(lrt, n = nrow(y[[1]]))[[1]]), topTags(lrt, n = nrow(y[[1]]))[[1]])
ret$ID = as.character(ret$ID)
ret = as.data.frame(ret)
metricList = list()
metricList[["K_L"]] = ret
metrics <- metricList[["K_L"]]
sigMets = metrics[which(metrics$FDR<0.001),]
sigK_Raw <- sigMets[which(sigMets$logFC<0),]
sigL_Raw <- sigMets[which(sigMets$logFC>0),]

# Get TMM data metrics 
rowNames = data[,1]
dataSel = data[,-1]
rownames(dataSel) = rowNames
group <- factor(c(1,1,1,2,2,2))
y <- DGEList(counts=dataSel,group=group)
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y, design)
fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
ret = data.frame(ID=rownames(topTags(lrt, n = nrow(y[[1]]))[[1]]), topTags(lrt, n = nrow(y[[1]]))[[1]])
ret$ID = as.character(ret$ID)
ret = as.data.frame(ret)
metricList = list()
metricList[["K_L"]] = ret
metrics <- metricList[["K_L"]]

sigMets = metrics[which(metrics$FDR<0.001),]
sigK_TMM <- sigMets[which(sigMets$logFC<0),]
sigL_TMM <- sigMets[which(sigMets$logFC>0),]
addDEG <- sigL_TMM[which(!sigL_TMM$ID %in% sigL_Raw$ID),]

# Filter, normalize, and standardize the data so each gene has mean=0 and stdev=1
res <- filterStandardizeKL(data)
# Fitered data standardized
fulls <- res[["fulls"]]
# Non-filtered data standardized
datas <- res[["datas"]]

# Set the color
totalColor = scales::seq_gradient_pal("purple", "purple4", "Lab")(seq(0,1,length.out=8))[4]

# Combine the filtered and remaining data
fulls <- datas
boxDat <- melt(fulls, id.vars="ID")
colnames(boxDat) <- c("ID", "Sample", "Count")

# Indices of the NAN rows
nID <- which(is.nan(datas$K.1))
# Set these filtered values that have all same values for samples to 0
datas[nID,1:6] <- 0

yMin = min(datas[,1:6])
yMax = max(datas[,1:6])

x = as.data.frame(datas[which(datas$ID %in% addDEG$ID),])
x$cluster = "color"
x$cluster2 = factor(x$cluster)
xNames = rownames(x)
xSig = x
xSigNames = rownames(xSig)
nGenes = nrow(xSig)

dendo = xSig[,1:6]
rownames(dendo) = NULL
d = dist(as.matrix(dendo))
hc = hclust(d, method="ward.D")

nC = 8
colList = scales::seq_gradient_pal("maroon1", "maroon4", "Lab")(seq(0,1,length.out=nC))
k = cutree(hc, k=nC)

plot_clusters = lapply(1:nC, function(j){
  i = rev(order(table(k)))[j]
  x = as.data.frame(xSig[,1:6][which(k==i),])
  nGenes = nrow(x)
  x$cluster = "color"
  x$cluster2 = factor(x$cluster)
  xNames = rownames(x)
  x$ID = xNames
  xSigNames = rownames(x)
  #saveRDS(xSigNames, file=paste0(getwd(), "/", outDir, "/Sig_", nC, "_", j, ".Rds"))
  
  pcpDat <- melt(x[,c(1:6,9)], id.vars="ID")
  colnames(pcpDat) <- c("ID", "Sample", "Count")
  boxDat$Sample <- as.character(boxDat$Sample)
  pcpDat$Sample <- as.character(pcpDat$Sample)
  
  p = ggplot(boxDat, aes_string(x = 'Sample', y = 'Count')) + geom_boxplot() + geom_line(data=pcpDat, aes_string(x = 'Sample', y = 'Count', group = 'ID'), colour = colList[j], alpha=0.3) + ylab("Standardized Count") + ggtitle(paste("Cluster ", j, " Genes (n=", format(nGenes, big.mark=",", scientific=FALSE), ")",sep="")) + theme_gray() + theme(plot.title = element_text(hjust = 0.5, size=8), axis.text=element_text(size=8), axis.title=element_text(size=8))
  p
})

# Plot all 8 clusters on a grid
do.call("grid.arrange", c(plot_clusters, ncol=ceiling(nC/4)))

