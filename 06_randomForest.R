
# Load libraries
require(pacman)
pacman::p_load(tidyverse, raster, rgdal, cclust, dismo, gtools, sp, rgeos, FactoMineR, pROC, randomForest, Hmisc, velox, sf, tmaptools)

# Initial setup
g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 9999)
run <- 'run2'
source('FunctionsRFclustering.R')
myproj <- CRS('+proj=longlat +datum=WGS84')

# Functions to use
rf.clust <- function(occ, nforest, ntrees, nVars, nclasses){
  # occ = back_swd; nforest = 50; ntrees = 500; nVars = 8; nclasses = 2
  datRF_presences <- occ[,3:ncol(occ)] %>% as.data.frame()
  print(nrow(datRF_presences))
  attach(datRF_presences)
  no.forests <- nforest
  no.trees <- ntrees
  distRF_presences <- RFdist(datRF_presences, mtry1 = nVars, no.trees, no.forests, addcl1 = T, addcl2 = F, imp = T, oob.prox1 = T)
  no.presencesclasses <- nclasses
  labelRF <- pamNew(distRF_presences$cl1, no.presencesclasses)
  print(table(labelRF))
  clusterdata <- hclust(as.dist(distRF_presences$cl1), method = 'ward.D2')
  return(list(labelRF, clusterdata))
}


# Load data ---------------------------------------------------------------
run <- 'run2'
vrs <- paste0('bio', 1:19)
vrs <- c(vrs, 'bio20_100mm', paste0('bio', 21:33))
vrs <- paste0(vrs, '.asc')
lyr <- list.files('../data/asc/clm/crn', pattern = '.asc$', full.names = TRUE) %>% 
  grep(paste0(vrs, collapse = '|'), ., value = TRUE) %>% 
  mixedsort() %>% 
  stack()
occ <- read.csv('../data/pnt/tbl/run2/occ_rmOtl.csv')
clusteredpresdata <- readRDS('../rds/run2/occ_cluster.rds')
msk <- raster('../data/asc/bse/msk.asc')

# Background process ------------------------------------------------------

# Background sampling probability
prd <- read_csv('../data/tbl/prd/FAOSTAT_data_5-15-2019.csv') %>% 
  as_tibble() %>% 
  dplyr::select(Area, Element, Item, Year, Unit, Value) %>% 
  group_by(Area, Element, Unit) %>% 
  dplyr::summarize(avg = mean(Value)) %>%
  ungroup() %>% 
  dplyr::filter(Element == 'Production') %>% 
  mutate(Area = as.character(Area))
adm <- st_read('../data/shp/all_countries.shp')
msk_shp <- st_read('../data/shp/mask_polygon.shp') 
adm <- crop_shape(adm, msk_shp)
occ_sft <- st_as_sf(occ, coords = c('X', 'Y'))
st_crs(occ_sft) <- st_crs(adm)
occ_sft <- st_intersection(occ_sft, adm) %>% 
  dplyr::select(ENGLISH)
count <- occ_sft %>% pull(ENGLISH) %>% as.character() %>% table() %>% as_tibble() %>% setNames(c('Country', 'n')) %>% mutate(Country = factor(Country))
adm <- full_join(adm, count, by = c('ENGLISH' = 'Country')) %>% 
  dplyr::select(NAME, COUNTRY, ENGLISH, n)

# Bias
bias <- inner_join(adm, prd, by = c('ENGLISH' = 'Area')) %>% 
  mutate(prob = n/avg)
bias.shape <- as(bias, 'Spatial')

# Velox
r <- lyr[[1]]
r[] <- 1
r <- velox(r)
r$rasterize(bias.shape, field = 'prob', band = 1)
rl <- r$as.RasterLayer(band = 1)
bias.raster <- rl
writeRaster(bias.raster, '../data/asc/bse/bias_raster_run2.tif', overwrite = TRUE)

# Bias raster - process
SPspecies <- SpatialPoints(occ[,1:2]) 
crs(SPspecies) <- myproj
back_raster <- bias.raster
speciescell <- raster::extract(bias.raster, SPspecies, cellnumber = TRUE)
back_raster[speciescell[,1]]  <- NA #remove the cell with presences
samplesize <- round(min(summary(as.factor(clusteredpresdata$cluster))) / 2, 0) 
NumberOfClusters <- max(clusteredpresdata$cluster) 
ratio <- NumberOfClusters/1
numberofpresences <- nrow(clusteredpresdata) #numberofpresences <- NumberOfClusters * samplesize

back_raster <- resample(back_raster, lyr[[21]]) %>%
  raster::crop(., lyr[[21]]) %>%
  raster::mask(., lyr[[21]])
crs(back_raster) <- myproj
back <- randomPoints(back_raster, 1*numberofpresences, prob = T) %>%
  as_data_frame()
coordinates(back) <- ~ x + y
back_swd  <- raster::extract(lyr, back) %>% cbind(coordinates(back), .)

nrow(back_swd) == nrow(back_swd[complete.cases(back_swd),])
write.csv(back_swd, paste0('../data/pnt/tbl/run2/back_swd.csv'), row.names = FALSE)
back_swd <- as.data.frame(back_swd)

coordinates(back_swd) <- ~ x + y
writeOGR(back_swd, dsn = '../data/pnt/shp/run2', layer = 'back_swd', driver = 'ESRI Shapefile')

back_swd <- as.data.frame(back_swd)

# Cluster analysis to pseudoabsences
bckclust <- rf.clust(occ = back_swd, nforest = 50, ntrees = 500, nVars = 8, nclasses = 2)
datRF <- as.data.frame(back_swd[,3:ncol(back_swd)])
attach(datRF)
no.forests <- 50#raw = 25
no.trees <- 500
distRF <- RFdist(datRF, mtry1 = 8, no.trees, no.forests, addcl1 = T, addcl2 = F, imp =T, oob.prox1 = T)# mtry1 = 4 raw  # es la cantidad de variables a utilizar en cada no
no.absenceclasses <- 2
labelRF <- pamNew(distRF$cl1, no.absenceclasses)
detach(datRF)
classdata <- cbind(pb = as.factor(labelRF), back_swd[,3:ncol(back_swd)])

cl <- pull(clusteredpresdata, 3)
clusteredpresdata <- raster::extract(lyr, clusteredpresdata[,1:2]) %>%
  cbind(clusteredpresdata[,1:2], ., cluster = cl)
head(clusteredpresdata)

# Join presences and background
presvalue_swd  <- clusteredpresdata[,3:ncol(clusteredpresdata)] %>%
  cbind(pb = (clusteredpresdata$cluster + no.absenceclasses), .) %>%
  na.omit() %>%
  as.data.frame() %>%
  mutate(cluster = cluster + no.absenceclasses)
unique(presvalue_swd$cluster)
presvalue_swd <- dplyr::select(presvalue_swd, pb, bio1:bio33)
presvalue_swd <- mutate(presvalue_swd, pb = as.factor(pb))

classdata_2 <- cbind(pb = as.data.frame(classdata)$pb, classdata[,2:ncol(classdata)]) # Background
dim(classdata_2); dim(presvalue_swd)
allclasses_swd <- rbind(classdata_2, presvalue_swd)
unique(allclasses_swd$pb)
head(allclasses_swd)

write.csv(allclasses_swd, paste0('../data/pnt/tbl/all_classes_swd.csv'), row.names = FALSE)

# Model Random Forest
vrs <- gsub('.asc', '', vrs)
model1 <- as.formula(paste('factor(pb) ~', paste(vrs, collapse = '+', sep =' ')))
rflist <- vector('list', 50) 
auc <- vector('list', 50)

for(repe in 1:50){ # 50 bosques
  
  print(repe)
  pressample <- list()
  
  for (i in 1:(NumberOfClusters+no.absenceclasses)){
    
    if(any(i==c(1:no.absenceclasses))) { 
      
      rows <- sample(rownames(allclasses_swd[allclasses_swd$pb==i,]), 
                     size = samplesize*NumberOfClusters/2/no.absenceclasses)
    } else {
      rows <- sample(rownames(allclasses_swd[allclasses_swd$pb==i,]), size=samplesize)
    }
    pressample[[i]] <- allclasses_swd[rows,] 
  }
  
  species <- na.omit(do.call(rbind, pressample)) 
  head(species)
  Samplesplit <- sample(rownames(species)) 
  
  envtrain <- species[Samplesplit[1:(0.8*nrow(species))],] 
  envtest <- species[Samplesplit[(0.8*nrow(species)):nrow(species)],] 
  
  rfmodel <- randomForest(model1, data = envtrain, ntree = 500, na.action = na.omit, nodesize = 2) 
  
  save(rfmodel, file = paste('../rf/', run, '/models/', NumberOfClusters, 'Prob_' , 'rep_' ,repe, '.rdata' ,sep=''))
  rflist[[repe]] <- rfmodel
  
  # AUC 
  predicted <- as.numeric(predict(rfmodel, envtest))
  observed <- as.vector(envtest[,'pb'])
  auc[[repe]] <- auc(observed, predicted) 
  rm(rfmodel)
  
  cat(auc[[repe]] ,'\n')
  
}

auc <- unlist(auc)
mean(auc)
rff <- do.call(randomForest::combine, rflist)
importance <- as.data.frame(rff$importance)

save(rflist, file = paste('../rds/', run, '/rflist_', NumberOfClusters, '.rdata', sep = ''))
save(importance, file = paste0('../rds/', run, '/importanceRF.rData'))
save(auc, file = paste0('../rds/', run, '/aucRF_dist.rData'))
save(rff, file = paste0('../rds/', run, '/rff_dist.rData'))

load('../_rds/_run1/aucRF_dist.rData')
auc
mean(auc)
