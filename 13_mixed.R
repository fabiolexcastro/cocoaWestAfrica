
# Load data
require(pacman)
pacman::p_load(raster, rgdal, tidyverse, rgeos, gtools, stringr, velox, sf, pROC, foreach, doSNOW)

# Initial setup
g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999)
run <- 'run6'

# Function to use
calcMixed <- function(yr){
  # yr <- '_2030'
  print(paste0('To start ', yr))
  prb <- paste0('../rf/', run, '/results/raw/', yr) %>% 
    list.files(full.names = T, pattern = '.asc$') %>% 
    grep('Prob', ., value = TRUE)
  cls <- paste0('../rf/', run, '/results/raw/', yr) %>% 
    list.files(full.names = T, pattern = '.tif$') %>% 
    grep('Clust_lim_', ., value = TRUE)
  unc <- paste0('../rf/', run ,'/results/raw/', yr) %>% 
    list.files(full.names = T, pattern = '.asc$') %>% 
    grep('Unc', ., value = TRUE)
  ncl <- 5
  nbk <- 2
  load(paste0('../rds/', run, '/thr.rData'))
  load(paste0('../rds/', run, '/thrUnc_.rData'))
  gcm <- list.files('../data/asc/clm/ftr/2020_2049')
  
  print('To the GCMs')
  lapply(1:length(prb), function(k){
    
    print(prb[k])
    pr <- raster(prb[k])
    cl <- raster(cls[k])
    un <- raster(unc[k])
    gc <- gcm[k]
    
    print(paste0('To calculate ', gc))
    
    rsl <- cl
    rsl[which(un[] < thrUnc & pr[] > thr)] <- max(unique(cl[]), na.rm = T) + 1 
    print('To write the raster')
    writeRaster(rsl, paste0('../rf/', run, '/results/raw/', yr, '/RF_', ncl, 'Classes_unc_', gc, '.asc'), overwrite = TRUE)
    print('Done!')
  })
  print('Done!!!')
}

# Apply the function
calcMixed(yr = '2020_2049')
calcMixed(yr = '2040_2069')

list.files('../rf/run6/results/raw/2020_2049/')
list.files('../rf/run6/results/raw/2040_2069/')








