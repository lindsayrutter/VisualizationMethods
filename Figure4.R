library(devtools)
library(EDASeq)
library(DESeq)
library(RCurl)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(GGally)
library(bigPint)
library(gridExtra)
library(cowplot)
install_github("drisso/yeastRNASeqRisso2011")
library(yeastRNASeqRisso2011)

source("functions.R")

# Read in three .rda files
githubURL <- "https://github.com/drisso/yeastRNASeqRisso2011/blob/master/data/"
load(url(paste0(githubURL, "geneLevelCounts.rda?raw=true")))
load(url(paste0(githubURL, "laneInfo.rda?raw=true")))
load(url(paste0(githubURL, "geneInfo.rda?raw=true")))

# Prepare data frame based on parameters from .rda files
data(yeastGC)
colnames(laneInfo)[2] <- "conditions"
means <- rowMeans(geneLevelCounts)
filter <- means >= 10
geneLevelCounts <- geneLevelCounts[filter,]
sub <- intersect(rownames(geneLevelCounts), names(yeastGC))
mat <- as.matrix(geneLevelCounts[sub, ])
data <- newSeqExpressionSet(mat, phenoData=laneInfo, featureData=AnnotatedDataFrame(data.frame(gc=yeastGC[sub])))

# Run within and between normalization on data frame
dataWithin <- withinLaneNormalization(data, "gc", which="full", offset=FALSE)
dataNorm <- betweenLaneNormalization(dataWithin, which="median")
counts <- as(dataWithin,"CountDataSet")
dataWithin <- as.data.frame(counts@assayData$counts)
counts <- as(dataNorm,"CountDataSet")
dataBetween <- as.data.frame(counts@assayData$counts)
dataWithin <- formatYeastDF(dataWithin)
dataWithin <- select(dataWithin, ID, Y1.1, Y1.2, Y4.1, Y4.2)
dataBetween <- formatYeastDF(dataBetween)
dataBetween <- select(dataBetween, ID, Y1.1, Y1.2, Y4.1, Y4.2)

# Use bigPint function to create scatterplot matrices showing within and between normalization
retWithin <- plotSM(dataWithin, option="allPoints", pointSize = 0.1, saveFile = FALSE)
retWithin <- retWithin[["Y1_Y4"]] + theme_gray()
retBetween <- plotSM(dataBetween, option="allPoints", pointSize = 0.1, saveFile = FALSE)
retBetween <- retBetween[["Y1_Y4"]] + theme_gray()

# Arrange two scatterplot matrices into grid
plot1 <- plot_grid(ggmatrix_gtable(retWithin), labels=c("A"), ncol = 1, nrow = 1, label_size=12) + theme(plot.background = element_rect(size=0.1,linetype="solid",color="black"))
plot2 <- plot_grid(ggmatrix_gtable(retBetween), labels=c("B"), ncol = 1, nrow = 1, label_size=12) + theme(plot.background = element_rect(size=0.1,linetype="solid",color="black"))
grid.arrange(plot1, plot2, ncol=2)
