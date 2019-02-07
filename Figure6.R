library(bigPint)
library(dplyr)
library(tidyr)
library(ggplot2)
library(limma)
library(gridExtra)

load("data/soybean_cn.rda")
data <- soybean_cn
load("data/sbCNSwitchedMetrics.rda")
data <- soybean_cn
metrics <- sbCNSwitchedMetrics

# Focus on treatment groups S1 and S2
data <- data[,c(1:7)]
dataSwitched <- data[,c(1,2,3,5,4,6,7)]
colnames(dataSwitched) <- c("ID","S1.1","S1.2","S1.3","S2.1","S2.2","S2.3")

# Produce four subplots below:

# 1. Boxplot of non-switched data
boxSel = data[,-1] %>% gather(Sample,Count)
boxSel$group = c(rep("S1",nrow(boxSel)/2), rep("S2", nrow(boxSel)/2))
ggplot(boxSel, aes(x=Sample, y=Count, fill=group)) + geom_boxplot() + scale_fill_manual(values=c("limegreen","magenta")) + theme(axis.text=element_text(size=12), axis.title=element_text(size=12), legend.position="none")

# 2. Boxplot of switched data
boxSel = dataSwitched[,-1] %>% gather(Sample,Count)
boxSel$group = c(rep("S1",nrow(boxSel)/2), rep("S2", nrow(boxSel)/2))
ggplot(boxSel, aes(x=Sample, y=Count, fill=group)) + geom_boxplot() + scale_fill_manual(values=c("limegreen","magenta")) + theme(axis.text=element_text(size=12), axis.title=element_text(size=12), legend.position="none")

# 3. MDS plot of non-switched data
colors <- c(rep(c("limegreen"), 3), rep(c("magenta"), 3))
plotMDS(data[,2:ncol(data)], col=colors, xlab="Dim 1", ylab="Dim 2", labels=colnames(data[,2:ncol(data)]), xlim=c(-4, 4), ylim=c(-4, 4))

# 4. MDS plot of switched data
plotMDS(dataSwitched[,2:ncol(dataSwitched)], col=colors, xlab="Dim 1", ylab="Dim 2", labels=colnames(dataSwitched[,2:ncol(dataSwitched)]), xlim=c(-4, 4), ylim=c(-4, 4))
