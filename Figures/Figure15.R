library(edgeR)
library(DESeq2)
library(bigPint)
library(data.table)
library(ggplot2)

load("../data/LK_data.RData")
data = as.data.frame(MA.subsetA$M)
rownames(data) = as.character(MA.subsetA$genes$EnsemblGeneID)
setDT(data, keep.rownames = TRUE)[]
colnames(data) = c("ID","K.R1L1","L.R1L2","K.R1L3","L.R1L4","L.R1L6","K.R1L7","L.R1L8","K.R2L2","L.R2L3","K.R2L6")
data = as.data.frame(data)
data = data[,c(1,2,4,7,9,11,3,5,6,8,10)]

# Obtain R1 values
data <- data[,c(1:4,7:9)]
colnames(data) <- c("ID", "K.1", "K.2", "K.3", "L.1", "L.2", "L.3")

# Get DEGs 
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

# Set all genes to have FDR=1 so they are not overplotted on scatterplot matrix
metricList2 = metricList
metricList2[["K_L"]]$FDR = 1

# Prepare data for plotting
logData <- data
rownames(logData) <- logData$ID
logData[,2:ncol(logData)] <- log(data[,2:ncol(logData)]+1)

ret <- plotSM(data=logData, dataMetrics = metricList2, option = "allPoints", threshVar = "FDR", threshVal=1e-3, pointSize = 0.1, saveFile = FALSE)
ret[["K_L"]] + theme_gray()
