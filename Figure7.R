library(bigPint)
library(ggplot2)
load("data/soybean_ir.rda")

soybean_ir[,-1] <- log(soybean_ir[,-1]+1)
ret = plotSM(soybean_ir, option ="allPoints", saveFile = FALSE, pointSize = 0.1)
ret[["N_P"]] + theme_gray()
