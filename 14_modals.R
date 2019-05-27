
# Load data
require(pacman)
pacman::p_load(raster, rgdal, tidyverse, rgeos, gtools, stringr, velox, sf, pROC, foreach, doSNOW)

# Initial setup
g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999)
run <- '_run5'

# Load data
yrs <- c('2020_2049', '2040_2069')

lapply(1:length(yrs), function(k){
  print(yrs[k])
  lyr <- paste0('../rf/', run, '/results/raw/', yrs[k]) %>% 
    list.files(full.names = T, pattern = '.asc$') %>% 
    grep('5Classes_unc_', ., value = T) %>%
    stack()
  mdl <- raster::modal(lyr)
  print('To write the final raster')
  writeRaster(mdl, paste0('../rf/', run, '/results/process/mixed/RF_5classes_unc_', yrs[k], '.asc'), overwrite = TRUE)
  print('Done!')
})
