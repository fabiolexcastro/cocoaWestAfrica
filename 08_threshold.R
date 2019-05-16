
# Load libraries
library(pacman)
pacman::p_load(raster, rgdal, cclust, gtools, sp, rgeos, stringr, foreach, doMC, doSNOW)

# Initial setup
g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999)
run <- 'run1'

# Function to use
occ <- read.csv('../data/pnt/tbl/occ_rmOtl.csv')
prb <- raster('../rf/run1/results/raw/RF_5Prob_current.asc')

# Extract values
vls <- raster::extract(prb, occ[,1:2])
vls <- vls[!is.na(vls)]
qnt <- quantile(vls, seq(0, 1, 0.05)) %>% as.data.frame()
thr <- as.numeric(subset(qnt, rownames(qnt) == '5%'))
save(thr, file = paste0('../rds/run1/thr.rData'))


