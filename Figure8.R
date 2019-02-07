library(ggplot2)

# Create dataset
time=c(rep("30", 2), rep("60", 2), rep("120", 2))
Expression=rep(c("Induced", "Repressed") , 3)
value=c(52, 43, 1421, 1682, 57, 4)
data=data.frame(time,Expression,value)
levels(data$time) <- c("30", "60", "120")

# Plot
ggplot(data, aes(fill=Expression, y=value, x=time)) + geom_bar(position="dodge", stat="identity") + geom_text(aes(label=value), size=3, position=position_dodge(width=0.9), vjust=-0.25) + xlab("Minutes after onset of iron stress") + ylab("Number of DEGs")
