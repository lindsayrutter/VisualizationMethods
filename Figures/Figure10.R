library(bigPint)

load("../data/soybean_ir.rda")
data <- soybean_ir
data <- data[,1:7]
data[,-1] <- log(data[,-1]+1)
datCol <- colnames(data)[-which(colnames(data) %in% "ID")]
myPairs <- unique(sapply(datCol, function(x) unlist(strsplit(x,"[.]"))[1]))
load("../data/cluster1_metrics.rda")

# Below is the code to examine litre plots of the four clusters. After the application opens, 1) set the "Treatment Pairs" option to "N" and "P", 2) set the "Metrics" to "FDR", and 3) set the "Metric order" option to "Increasing"

# Cluster 1
app <- plotLitreApp(data = data, dataMetrics = metrics, pointColor = '#00BF7D')
if (interactive()) {
  shiny::runApp(app, port = 1234, launch.browser = TRUE)
}

# Cluster 2
load("data/cluster2_metrics.rda")
app <- plotLitreApp(data = data, dataMetrics = metrics, pointColor = '#A3A500')
if (interactive()) {
  shiny::runApp(app, port = 1234, launch.browser = TRUE)
}

# Cluster 3
load("data/cluster3_metrics.rda")
app <- plotLitreApp(data = data, dataMetrics = metrics, pointColor = '#E76BF3')
if (interactive()) {
  shiny::runApp(app, port = 1234, launch.browser = TRUE)
}

# Cluster 4
load("data/cluster4_metrics.rda")
app <- plotLitreApp(data = data, dataMetrics = metrics, pointColor = '#F8766D')
if (interactive()) {
  shiny::runApp(app, port = 1234, launch.browser = TRUE)
}
