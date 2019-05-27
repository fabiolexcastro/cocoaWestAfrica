
# Load data
require(pacman)
pacman::p_load(raster, rgdal, tidyverse, rgeos, gtools, stringr, velox, sf, pROC, foreach, doSNOW)

# Initial setup
g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999)
run <- 'run6'

# Function to use
calcLimitations <- function(yr){
  # yr <- '_2030'
  print(paste0('To start ', ' ', yr))
  prb <- paste0('../rf/', run, '/results/raw/', yr) %>% 
    list.files(full.names = T, pattern = '.asc$') %>% 
    grep('Prob', ., value = TRUE)
  cls <- paste0('../rf/', run, '/results/raw/', yr) %>% 
    list.files(full.names = T, pattern = '.asc$') %>% 
    grep('Clust', ., value = TRUE)
  ncl <- 5
  nbk <- 2
  load(paste0('../rds/', run, '/thr.rData'))
  
  print('The matrix')
  mprb <- matrix(c(0, thr, 0, thr, 1, 2), ncol = 3, byrow = T)
  mcls <- matrix(c(0, nbk + 0.5, 0, nbk + 0.5, nbk + ncl + 0.5, 1), ncol = 3, byrow = T)
  gcms <- list.files('../data/asc/clm/ftr/2020_2049')
  
  print('To the GCMs')
  lapply(1:length(prb), function(k){
    
    print(prb[k])
    pr <- raster(prb[k])
    cl <- raster(cls[k])
    gc <- gcms[k]
    
    cl.rc <- raster::reclassify(cl, mcls)
    pr.rc <- raster::reclassify(pr, mprb)
    
    print(paste0('To make the difference ', gc))
    df <- pr.rc - cl.rc
    rs <- cl
    
    rs[which(df[] == -1)] <- nbk + ncl + 1
    rs[which(df[] == 2)]  <- nbk + ncl + 1
    
    print('To write the raster')
    writeRaster(rs, paste0('../rf/', run, '/results/raw/', yr, '/RF_', ncl, '_Clust_lim_', gc, '.tif'), overwrite = TRUE)
    
  })
  print('Done!!!')
}

# Apply the function
calcLimitations(yr = '2020_2049')
calcLimitations(yr = '2040_2069')

# Review the results
rslts <- list.files(paste0('../rf/', run, '/results/raw/2020_2049')) %>% grep('_lim_', ., value = TRUE)
rslts <- list.files(paste0('../rf/', run, '/results/raw/2040_2069')) %>% grep('_lim_', ., value = TRUE)
length(rslts)



