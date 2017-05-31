# Name the working directory where your file is located:
setwd("./data-raw/catmap163")

library(catmap)
data(catmapdata)

# Build cm.obj from file location
cm.obj <- catmap(catmapdata, 0.95, TRUE)

# Call catmap functions
catmap.forest(cm.obj, TRUE, TRUE)
catmap.sense(cm.obj, TRUE, TRUE, TRUE)
catmap.cumulative(cm.obj, TRUE, TRUE, TRUE)
catmap.funnel(cm.obj, TRUE)
