

# Load libraries
require(pacman)
pacman::p_load(raster, rgdal, rgeos, stringr, tidyverse, velox, sf, rasterVis)

# Initial setup 
g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999)
run <- 'run6'

# Function to use
impGra <- function(crn, ftr){
  
  # crn <- crn
  # ftr <- f50
  
  msk <- crn * 0
  crd_df <- coordinates(crn)
  
  x <- raster::extract(crn, crd_df, cellnumbers = TRUE) %>% as_data_frame()
  ncell <- dplyr::select(x, cells)
  x <- select_(x, names(crn))
  colnames(x) <- 'current'
  
  x <- raster::extract(crn, crd_df, cellnumbers = TRUE) %>% as_data_frame()
  ncell <- dplyr::select(x, cells)
  x <- select_(x, names(crn))
  colnames(x) <- 'current'
  
  y <- raster::extract(ftr, crd_df[,c('x', 'y')], cellnumbers = TRUE) %>% as_data_frame()
  y <- select_(y, names(ftr))
  colnames(y) <- 'future'
  
  z <- data.frame(x, y, ncell) %>% as_tibble()
  
  print('To Results')
  rslts <- left_join(z, all_options, by = c('current', 'future'))
  labls <- as_tibble(labelss) %>% mutate(category = as.character(category))
  
  final <- left_join(rslts, labls, by = 'category') %>%
    dplyr::select(value) %>%
    pull(1)
  
  print(length(final))
  print(length(msk))
  
  rst <- raster::setValues(msk, final)
  print('Done!!!')
  return(rst)
  
}

# Load data
crn <- raster(paste0('../rf/', run, '/results/raw/RF_5Classes_unc_current.asc'))
f30 <- raster(paste0('../rf/', run, '/results/process/mixed/RF_5classes_unc_2020_2049.asc'))
f50 <- raster(paste0('../rf/', run, '/results/process/mixed/RF_5classes_unc_2040_2069.asc'))

all_options <- read_csv('../data/tbl/imp/classesImpGraLimMix.csv')
labelss <- data.frame(value = c(0, 1, 2, 3, 4, 5), category = c('Unsuit', 'cope', 'adjust', 'transform', 'opportunity', 'resilience'))

# Apply the function
i30 <- impGra(crn = crn, ftr = f30)
i50 <- impGra(crn = crn, ftr = f50)

i30 <- rst
i50 <- rst

dir.create('../rf/run6/results/process/impGra')
writeRaster(x = i30, filename = paste0('../rf/', run, '/results/process/impGra/impGra_30.tif'), overwrite = T)
writeRaster(x = i50, filename = paste0('../rf/', run, '/results/process/impGra/impGra_50.tif'), overwrite = T)



