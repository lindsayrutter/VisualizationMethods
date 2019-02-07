library(bigPint)
library(dplyr)
library(tidyr)
library(ggplot2)
library(limma)
library(GGally)
library(gridExtra)
library(cowplot)

load("data/soybean_cn.rda")
data <- soybean_cn
load("data/soybean_cn_metrics.rda")
metrics <- soybean_cn_metrics

# Focus on treatment groups S1 and S2
data <- data[,c(1:7)]

metricList = list()
metricList[["S1_S2"]] = metrics[["S1_S2"]]

p1 = plotSM(data, metricList, option="allPoints", threshVar="FDR", threshVal=0.05, pointSize=0.1, saveFile = FALSE)
p2 = plotPCP(data, metricList, threshVar="FDR", threshVal=0.05, lineSize=0.1, saveFile = FALSE)

plot1 <- ggmatrix_gtable(p1[["S1_S2"]])
plot2 <- p2[["S1_S2"]] + theme(axis.text=element_text(size=10), axis.title =element_text(size=10))
plot_grid(plot1, plot2, nrow = 2, rel_heights = c(0.7, 0.3), labels = c("A", "B"), label_size = 9)

