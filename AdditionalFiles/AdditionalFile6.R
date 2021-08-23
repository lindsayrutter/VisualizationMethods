library(bigPint)

load("../data/kidneyLiver.rda")
dat <- data
dat <- dat[,1:7]
dat[,-1] <- log(dat[,-1]+1)
load("../data/orig1_metrics.rda")

# Below is the code to examine litre plots. After the application opens, 1) set the "Treatment Pairs" option to "K" and "L", 2) set the "Metrics" to "FDR", and 3) set the "Metric order" option to "Increasing"

app <- plotLitreApp(data = dat, dataMetrics = metrics, pointColor = '#FF8C00')
if (interactive()) {
    shiny::runApp(app, port = 1234, launch.browser = TRUE)
}
