library(bigPint)
library(dplyr)
library(tidyr)
library(ggplot2)
library(limma)

load("data/soybean_cn.rda")
data <- soybean_cn

# Focus on treatment groups S1 and S2
data <- data[,c(1:7)]
p = plotSM(data, option="allPoints", pointSize=0.25, saveFile = FALSE)

# Plot figure
p[["S1_S2"]]
