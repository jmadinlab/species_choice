# Anyalsis
library(tripack)
source("R/functions.R")

# points selection

n <- 100 # total species
s <- 20  # species selected
dat <- data.frame(x=rnorm(n), y=rnorm(n))
plot(y ~ x, dat, col="grey")

# random
points(y ~ x, dat[sample(1:n, s, replace=FALSE),], col="red", pch=20, cex=0.6)

# most evenly spread
points(y ~ x, voronoiFilter(dat, s), col="blue", pch=20, cex=0.6)

