library(bigPint)
library(dplyr)
library(tidyr)
library(ggplot2)
library(limma)

load("data/soybean_cn.rda")
data <- soybean_cn
load("data/sbCNSwitchedMetrics.rda")
metrics <- sbCNSwitchedMetrics

# Focus on treatment groups S1 and S2
data <- data[,c(1:7)]
data <- data[,c(1,2,3,5,4,6,7)]
colnames(data) <- c("ID","S1.1","S1.2","S1.3","S2.1","S2.2","S2.3")

p = plotSM(data, option="allPoints", pointSize=0.25, saveFile = FALSE)

# Plot figure
p[["S1_S2"]]
