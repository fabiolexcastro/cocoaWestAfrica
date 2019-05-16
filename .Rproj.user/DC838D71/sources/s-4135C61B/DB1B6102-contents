

# Load libraries
library(pacman)
pacman::p_load(raster, rgdal, cclust, gtools, sp, rgeos, stringr, foreach, doMC, doSNOW)

# Initial setup
g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999)
run <- 'run1'

# Function to use ---------------------------------------------------------
calcLimitations <- function(){
  # sp <- sps[1]
  prb <- paste0('../rf/run1/results/raw/') %>% 
    list.files(full.names = T, pattern = '.asc$') %>% 
    grep('Prob', ., value = TRUE) %>% 
    raster()
  cls <- paste0('../rf/run1/results/raw/') %>% 
    list.files(full.names = T, pattern = '.asc$') %>% 
    grep('Clust_c', ., value = TRUE) %>% 
    raster()
  ncl <- 5
  nbk <- 2
  
  print('The matrix')
  mprb <- matrix(c(0, thr, 0, thr, 1, 2), ncol = 3, byrow = T)
  mcls <- matrix(c(0, nbk + 0.5, 0, nbk + 0.5, nbk + ncl + 0.5, 1), ncol = 3, byrow = T)
  
  print('To reclassify')
  prb_rcl <- raster::reclassify(prb, mprb)
  cls_rcl <- raster::reclassify(cls, mcls)
  
  print('To make the difference')
  dff <- prb_rcl - cls_rcl
  rsl <- cls
  rsl[which(dff[] == -1)] <- nbk + ncl + 1
  rsl[which(dff[] == 2)]  <- nbk + ncl + 1
  
  print('To write the rasters')
  writeRaster(rsl, paste0('../rf/', run, '/results/raw/', 'RF_', ncl, '_Clust_lim_crn.asc'), overwrite = T)
  print('Done!')
}


# Load data ---------------------------------------------------------------
run <- 'run1'
load('../rds/run1/thr.rData')

