

# Load libraries
library(pacman)
pacman::p_load(raster, rgdal, cclust, outliers, dismo, gtools, multcomp, sp, rgeos, outliers, FactoMineR, pROC, randomForest, stringr)

# Initial setup
g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 9999)
run <- '_run5'
source('FunctionsRFclustering.R')
myproj <- CRS('+proj=longlat +datum=WGS84')
run <- 'run1'

# Load data
clusterdata <- readRDS(paste0('../rds/', run, '/clusterdata.rds'))
load(paste0('../rds/', run, '/rflist_5.rdata'))
rff <- do.call(randomForest::combine, rflist)
vrs <- paste0(vrs, '.asc')


gcmlist <- 'current'
ar5biofolder <- '../data/asc/clm/crn/'
resultsfolder <- paste0('../rf/', run, '/results/raw')
modelfolder<- paste0('../rf/', run, '/models')
gcm <- gcmlist
gcmfiles <- list.files(ar5biofolder, full.names = T, pattern = '.asc$') %>% 
  mixedsort() %>% 
  grep(paste0(vrs, collapse = '|'), ., value = T)
climatelayers <- stack(gcmfiles)
climatevalues <- data.frame(getValues(climatelayers))
NumberOfClusters <- 5

# To predict
rasterProbs <- predict(rff, climatevalues, type = 'prob')
rasterProbs_na <- na.omit(rasterProbs)
sum_rasterProbs_na <- apply(rasterProbs_na, 1, sum)

rasterRF <- rowSums(rasterProbs[,c(3:(NumberOfClusters+2))])
uncertainty <- apply(rasterProbs, 1, max)  

rasterRFprob <- climatelayers[[1]]
values(rasterRFprob) <- rasterRF 

rasterRFuncertainty <- climatelayers[[1]]
values(rasterRFuncertainty) <- uncertainty 

rasterRF <- max.col(rasterProbs, 'first')
rasterRFclass <- climatelayers[[1]]
values(rasterRFclass) <- rasterRF

# Write the final rasters
plot(rasterRFclass)
run
writeRaster(rasterRFclass, paste0('../rf/', run, '/results/raw/RF_5Clust_current.asc'), format = 'ascii', overwrite = T)
writeRaster(rasterRFprob, paste0('../rf/', run, '/results/raw/RF_5Prob_current.asc'), format = 'ascii', overwrite = T)
writeRaster(rasterRFuncertainty, paste0('../rf/', run, '/results/raw/RF_5Unc_current.asc'), format = 'ascii', overwrite = T)

 