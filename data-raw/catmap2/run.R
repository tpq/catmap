# Name the working directory where your file is located:
setwd("./data-raw/catmap2")

source("catmap2.R")
catmapdata <- read.table("catmapdata.txt", header = TRUE, sep = "\t")

# Build cm.obj from file location
cm.obj <- catmap(catmapdata, 0.95, TRUE, TRUE)

# Call catmap functions
catmap.forest(cm.obj, TRUE, TRUE)
catmap.sense(cm.obj, TRUE, TRUE, TRUE, TRUE)
catmap.cumulative(cm.obj, TRUE, TRUE, TRUE, TRUE)
catmap.funnel(cm.obj, TRUE)
