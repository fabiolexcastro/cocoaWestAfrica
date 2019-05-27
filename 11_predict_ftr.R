

# Load libraries
library(pacman)
pacman::p_load(raster, rgdal, cclust, outliers, dismo, gtools, multcomp, sp, rgeos, outliers, FactoMineR, pROC, randomForest, stringr)

# Initial setup
g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 9999)
source('FunctionsRFclustering.R')
myproj <- CRS('+proj=longlat +datum=WGS84')
run <- 'run6'

# Load data
vrs <- readRDS(paste0('../rds/', run , '/vrs.rds'))
clusterdata - readRDS(paste0('../rds/', run, '/clusterdata.rds'))
load(paste0('../rds/', run, '/rflist_5.rdata'))

NumberOfClusters <- 5
ar5biofolder <- '../data/asc/clm/ftr'
yearlist <- list.files('../data/asc/clm/ftr')
gcmlist <- list.files('../data/asc/clm/ftr/2020_2049')
resultsfolder <- paste0('../rf/', run, '/results/raw')

rff <- do.call(randomForest::combine, rflist)

dir.create('../temp')
rasterOptions(tmpdir = '../temp')
cl <- makeCluster(19) 
registerDoSNOW(cl) 
y <- 2

foreach(i = 1:length(gcmlist), .packages = c('raster', 'rgdal', 'dplyr', 'gtools', 'foreach', 'randomForest', 'sp', 'stringr'), .verbose = TRUE) %dopar% {
  
  print(gcmlist[i]) 
  
  gcmfiles <- paste(ar5biofolder, yearlist[y], gcmlist[i], sep = '/') %>%
    list.files(., full.names = T, pattern = '.asc$') %>% 
    grep(paste0(vrs, collapse = '|'), ., value = T) %>%  
    mixedsort()
  
  climatelayers <- raster::stack(gcmfiles)
  climatevalues <- data.frame(getValues(climatelayers))
  
  print('Climate values')
  rasterProbs <- predict(rff, climatevalues, type = 'prob') # proximity = T
  rasterRF <- rowSums(rasterProbs[,c(3:(NumberOfClusters+2))])
  uncertainty <- apply(rasterProbs, 1, max)  
  
  rasterRFprob <- climatelayers[[1]]
  values(rasterRFprob) <- rasterRF 
  
  rasterRFuncertainty <- climatelayers[[1]]
  values(rasterRFuncertainty) <- uncertainty 
  
  rasterRF <- max.col(rasterProbs, 'first')
  rasterRFclass <- climatelayers[[1]]
  values(rasterRFclass) <- rasterRF
  
  # Write Raster
  print("Write Raster...")
  
  writeRaster(rasterRFclass, paste0('../rf/', run, '/results/raw/',  yearlist[y], '/RF_', NumberOfClusters, 'Clust_', gcmlist[i], '_', yearlist[y], '.asc', sep=''),  format = 'ascii', overwrite = T)
  writeRaster(rasterRFprob, paste0('../rf/', run, '/results/raw/', yearlist[y], '/RF_', NumberOfClusters, 'Prob_',  gcmlist[i], '_', yearlist[y], '.asc', sep=''),  format = 'ascii', overwrite = T)
  writeRaster(rasterRFuncertainty, paste('../rf/', run, '/results/raw/', yearlist[y], '/RF_', NumberOfClusters, 'Unc_', gcmlist[i], '_', yearlist[y], '.asc', sep=''),  format = 'ascii', overwrite = T)
  
  print('Done!')
  print(gcmlist[i])
  
}
stopCluster(cl)

dir.create('../rf/run6/results/raw/2020_2049')
dir.create('../rf/run6/results/raw/2040_2069')

