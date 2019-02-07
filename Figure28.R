library(bigPint)
library(EDASeq)
library(dplyr)

load("data/kidneyLiver.rda")
dat <- data
dat <- dat[,1:7]

# Standardize in this application
RowSD = function(x) {
    sqrt(rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1))
}
dat_Rownames <- dat$ID
dat = dat[,-1]
rownames(dat) <- dat_Rownames
dat <- betweenLaneNormalization(as.matrix(dat), which="full", round=FALSE)
dat = as.data.frame(dat)
dat = mutate(dat, mean = (K.1+K.2+K.3+L.1+L.2+L.3)/6, stdev = RowSD(cbind(K.1,K.2,K.3,L.1,L.2,L.3)))
rownames(dat)=dat_Rownames
dat$ID <- dat_Rownames
datqps <- t(apply(as.matrix(dat[,1:6]), 1, scale))
datqps <- as.data.frame(datqps)
colnames(datqps) <- colnames(dat[,1:6])
datqps$ID <- rownames(datqps)
nID <- which(is.nan(datqps$K.1))
datqps[nID,1:6] <- 0
dat <- datqps[,c(7,1:6)]

load("data/add1_metrics.rda")

# Below is the code to examine litre plots. After the application opens, 1) set the "Treatment Pairs" option to "K" and "L", 2) set the "Metrics" to "FDR", and 3) set the "Metric order" option to "Increasing"

app <- plotLitreApp(data = dat, dataMetrics = metrics, pointColor = '#FF34B3')
if (interactive()) {
    shiny::runApp(app, port = 1234, launch.browser = TRUE)
}
