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
library(dplyr)
library(cowplot)

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

# Get data metrics
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

# Determine DEGs that are liver-specific and kidney-specific
sigMets = metrics[which(metrics$FDR<0.001),]
sigL <- sigMets[which(sigMets$logFC<0),]
sigK <- sigMets[which(sigMets$logFC>0),]

logSoy = data
logSoy[,-1] <- log(data[,-1]+1)

# Filter, normalize, and standardize the data so each gene has mean=0 and stdev=1
res <- filterStandardizeKL(data)
# Fitered data standardized
fulls <- res[["fulls"]]
# Non-filtered data standardized
datas <- res[["datas"]]

# Choose colors for plots
colPurple = scales::seq_gradient_pal("purple4", "purple", "Lab")(seq(0,1,length.out=9))
colOrange = scales::seq_gradient_pal("orangered4", "darkorange2", "Lab")(seq(0,1,length.out=9))
colList = c(colPurple[5], colOrange[5])

Type = c("Kidney", "Liver")
yMin = min(datas[,1:6])
yMax = max(datas[,1:6])
sigIDs = list(sigL$ID, sigK$ID)

# Create background boxplot data
boxDat <- melt(fulls, id.vars="ID")
colnames(boxDat) <- c("ID", "Sample", "Count")

plot_clustersSigRaw = lapply(1:2, function(i){ 
  x = as.data.frame(datas[which(datas$ID %in% sigIDs[[i]]),])
  x$cluster = "color"
  x$cluster2 = factor(x$cluster)
  xNames = rownames(x)
  metricFDR = metrics[which(as.character(metrics$ID) %in% xNames),]
  sigID = metricFDR[metricFDR$FDR<0.001,]$ID
  xSig = x[which(rownames(x) %in% sigID),]
  xSigNames = rownames(xSig)
  nGenes = nrow(xSig)
  
  xSig$ID = xSigNames
  pcpDat <- melt(xSig[,c(1:7)], id.vars="ID")
  colnames(pcpDat) <- c("ID", "Sample", "Count")
  pcpDat$Sample <- as.character(pcpDat$Sample)
  
  pSig = ggplot(boxDat, aes_string(x = 'Sample', y = 'Count')) + geom_boxplot() +
    geom_line(data=pcpDat, aes_string(x = 'Sample', y = 'Count', group = 'ID'), colour = colList[i], alpha=0.05) + 
    ylab("Standardized Count") +
    ggtitle(paste("Significant ",  Type[i] ," Genes (n=", format(nGenes, big.mark=",", scientific=FALSE), ")",sep="")) + theme_gray() +
    theme(plot.title = element_text(hjust = 0.5, size=12), axis.text=element_text(size=12), axis.title=element_text(size=12))
  pSig
})

# Repeat for kidney
#######################################################

y <- DGEList(counts=dataSel,group=group)
# Only difference is we know have calcNormFactors
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

# Determine DEGs that are liver-specific and kidney-specific
sigMets = metrics[which(metrics$FDR<0.001),]
sigL <- sigMets[which(sigMets$logFC<0),]
sigK <- sigMets[which(sigMets$logFC>0),]

logSoy = data
logSoy[,-1] <- log(data[,-1]+1)

# Filter, normalize, and standardize the data so each gene has mean=0 and stdev=1
res <- filterStandardizeKL(data)
# Fitered data standardized
fulls <- res[["fulls"]]
# Non-filtered data standardized
datas <- res[["datas"]]

# Choose colors for plots
colPurple = scales::seq_gradient_pal("purple4", "purple", "Lab")(seq(0,1,length.out=9))
colOrange = scales::seq_gradient_pal("orangered4", "darkorange2", "Lab")(seq(0,1,length.out=9))
colList = c(colPurple[5], colOrange[5])

Type = c("Kidney", "Liver")
yMin = min(datas[,1:6])
yMax = max(datas[,1:6])
sigIDs = list(sigL$ID, sigK$ID)

# Create background boxplot data
boxDat <- melt(fulls, id.vars="ID")
colnames(boxDat) <- c("ID", "Sample", "Count")

plot_clustersSigTMM = lapply(1:2, function(i){ 
  x = as.data.frame(datas[which(datas$ID %in% sigIDs[[i]]),])
  x$cluster = "color"
  x$cluster2 = factor(x$cluster)
  xNames = rownames(x)
  metricFDR = metrics[which(as.character(metrics$ID) %in% xNames),]
  sigID = metricFDR[metricFDR$FDR<0.001,]$ID
  xSig = x[which(rownames(x) %in% sigID),]
  xSigNames = rownames(xSig)
  nGenes = nrow(xSig)
  
  xSig$ID = xSigNames
  pcpDat <- melt(xSig[,c(1:7)], id.vars="ID")
  colnames(pcpDat) <- c("ID", "Sample", "Count")
  pcpDat$Sample <- as.character(pcpDat$Sample)
  
  pSig = ggplot(boxDat, aes_string(x = 'Sample', y = 'Count')) + geom_boxplot() +
    geom_line(data=pcpDat, aes_string(x = 'Sample', y = 'Count', group = 'ID'), colour = colList[i], alpha=0.05) + 
    ylab("Standardized Count") +
    ggtitle(paste("Significant ",  Type[i] ," Genes (n=", format(nGenes, big.mark=",", scientific=FALSE), ")",sep="")) + theme_gray() +
    theme(plot.title = element_text(hjust = 0.5, size=12), axis.text=element_text(size=12), axis.title=element_text(size=12))
  pSig
})

plot1 <- plot_clustersSigRaw[[1]]
plot2 <- plot_clustersSigRaw[[2]]
plot3 <- plot_clustersSigTMM[[1]]
plot4 <- plot_clustersSigTMM[[2]]

plot_grid(plot1, plot2, plot3, plot4, ncol = 2, labels = c("A", "A", "B", "B"), label_size = 12)
