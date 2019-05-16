
# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(raster, rgdal, rgeos, gtools, tidyverse, stringr, velox, sf, foreach, doSNOW, parallel)

# Initial setup -----------------------------------------------------------
g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999)
vrs <- c(paste0('prec_', 1:12, '$'),
         paste0('tmax_', 1:12, '$'),
         paste0('tmin_', 1:12, '$'),
         paste0('tmean_', 1:12, '$'))
setwd('//mnt/workspace_cluster_9/Coffee_Cocoa2/_cocoaAfrica/r')

# Functions to use --------------------------------------------------------
extMask <- function(vr){
  # vr <- 'prec'
  fl <- grep(vr, fls, value = TRUE)
  st <- stack(fl)
  ct <- raster::crop(st, ext)
  ct <- unstack(ct)
  Map('writeRaster', x = ct, filename = paste0('../data/asc/clm/crn/', basename(fl), '.asc'), overwrite = TRUE)
  return(ct)
}
extMaskFtr <- function(gc, yr){
  # gc <- gcm[1]
  # yr <- yrs[1]
  print(paste0(yr, ' ', gc))
  gc <- paste0(gc, '/')
  fl <- grep(yr, fls, value = TRUE) %>%
    grep(gc, ., value = TRUE) %>%
    grep(paste0(vrs, collapse = '|'), ., value = TRUE)
  st <- stack(fl) 
  ct <- raster::crop(st, extent(msk))
  ct <- unstack(ct)
  dir.create(paste0('../data/asc/clm/ftr/', yr, '/', gc))
  pt <- paste0('../data/asc/clm/ftr/', yr, '/', gc)
  Map('writeRaster', x = ct, filename = paste0(pt, '/', basename(fl), '.asc'), overwrite = TRUE)
  print('Done!')
}

# Load extent --------------------------------------------------------------
ext <- readRDS(file = '../rds/extent.rds')

# Extract by mask - Current -----------------------------------------------
pth <- '//mnt/data_cluster_4/observed/gridded_products/worldclim/Global_2_5min'
vrs <- c(paste0('prec_', 1:12, '$'), paste0('tmin_', 1:12, '$'), paste0('tmean_', 1:12, '$'), paste0('tmax_', 1:12, '$'))
fls <- list.files(pth, full.names = TRUE) %>%
  grep(paste0(vrs, collapse = '|'), ., value = TRUE) %>%
  mixedsort()
var <- c('prec', 'tmin', 'tmean', 'tmax')
cl <- makeCluster(4)
registerDoSNOW(cl)
registerDoMC(4)
foreach(i = 1:length(var), .packages = c('raster', 'rgdal', 'gtools', 'foreach', 'sp', 'stringr', 'tidyverse'), .verbose = TRUE) %dopar% {
  print(var[i])
  extMask(vr = var[i])  
}
stopCluster(cl)

# Review ------------------------------------------------------------------
crn <- list.files('../data/asc/clm/crn', full.names = TRUE, pattern = '.asc$') %>% 
  stack()
msk <- crn[[1]] * 0 + 1
writeRaster(msk, '../data/asc/bse/msk.asc', overwrite = TRUE)

# Extract by mask future --------------------------------------------------
pth <- '//mnt/data_cluster_2/gcm/cmip5/downscaled/rcp60/global_2_5min/'
gcm <- list.files(pth)
yrs <- c('2020_2049', '2040_2069')
fls <- paste0(pth, '/',  gcm)
f30 <- paste0(fls, '/', 'r1i1p1', '/', yrs[1])
f50 <- paste0(fls, '/', 'r1i1p1', '/', yrs[2])
f30 <- list.files(f30, full.names = T) %>% 
  grep(paste0(vrs, collapse = '|'), ., value = T) %>% 
  mixedsort()
f50 <- list.files(f50, full.names = T) %>% 
  grep(paste0(vrs, collapse = '|'), ., value = T) %>% 
  mixedsort()
fls <- c(f30, f50)

registerDoMC(19)
foreach(i = 1:length(gcm), .packages = c('raster', 'rgdal', 'gtools', 'foreach', 'sp', 'stringr', 'tidyverse'), .verbose = TRUE) %dopar% {
  foreach(y = 1:length(yrs)) %do% {
    extMaskFtr(gc = gcm[i], yr = yrs[y])  
  }
}








