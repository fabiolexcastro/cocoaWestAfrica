

# Load data
require(pacman)
pacman::p_load(raster, rgdal, tidyverse, rgeos, gtools, stringr, velox, sf, pROC, foreach, doSNOW)

# Initial setup
g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999)
run <- 'run1'

# Function to use
calcMixed <- function(){
  # sp <- sps[1]
  print(paste0('To load the data'))
  prb <- paste0('../rf/', run, '/results/raw/') %>% 
    list.files(full.names = T, pattern = '.asc$') %>% 
    grep('Prob', ., value = TRUE) %>% 
    raster()
  cls <- paste0('../rf/', run, '/results/raw/') %>% 
    list.files(full.names = T, pattern = '.asc$') %>% 
    grep('Clust_lim_crn.asc', ., value = TRUE) %>% 
    raster()
  unc <- paste0('../rf/', run, '/results/raw/') %>% 
    list.files(full.names = T, pattern = '.asc$') %>% 
    grep('Unc', ., value = TRUE) %>% 
    raster()
  ncl <- 5
  occ <- readRDS('../rds/run1/occ_cluster.rds')
  load('../rds/run1/thr.rData')
  
  print('To make the extraction')
  thrUnc <- raster::extract(unc, occ[,c('X', 'Y')])
  thrUnc <- thrUnc[!is.na(thrUnc)]
  thrUnc <- quantile(thrUnc, 0.05) %>% as.numeric()
  save(thrUnc, file = paste0('../rds/', run, '/thrUnc_', '.rData'))
  load('../rds/run5/thrUnc_.rData')
  rsl <- cls
  rsl[which(unc[] < thrUnc & prb[] > thr)] <- max(unique(cls[]), na.rm = T) + 1 
  print('To write the raster')
  writeRaster(rsl, paste0('../rf/', run,'/results/raw/', 'RF_', ncl, 'Classes_unc_', gcm, '.asc'), overwrite = TRUE)
  print('Done!')
}

# Load data
ncl <- 5
gcm <- 'current'

# Apply the functions
calcMixed()








