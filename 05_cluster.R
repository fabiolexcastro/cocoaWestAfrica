
# Load libraries
require(pacman)
pacman::p_load(raster, rgdal, rgeos, stringr, velox, sf, randomForest, outliers, gtools, tidyverse)

# Initial setup -----------------------------------------------------------
g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999)
source('FunctionsRFclustering.R')

# Functions to use --------------------------------------------------------
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
dup_cell <- function(msk, occ){
  
  src <- as_data_frame(occ) %>% pull(1)
  occ <- st_coordinates(x = occ) %>% as.data.frame()
  occ$source <- src
  
  cellNum <- raster::extract(msk, occ[,1:2], cellnumbers = T) 
  cells <- xyFromCell(msk, cellNum[,'cells'])
  dupvec <- duplicated(cells[,c('x', 'y')])
  occ_rmDupCell <- tbl_df(occ[!dupvec,])
  occ_DupCell <- tbl_df(occ[dupvec,])
  return(list(occ_rmDupCell, occ_DupCell))
}

# Load data ---------------------------------------------------------------
pnt <- st_read('../data/shp/pnt/run1/points_raw.shp')
msk <- raster('../data/asc/bse/msk.asc')
msk_shp <- shapefile('../data/shp/mask_polygon.shp')
vrs <- paste0('bio', 1:19)
vrs <- c(vrs, 'bio20_100mm', paste0('bio', 21:33))
vrs <- paste0(vrs, '.asc')
lyr <- list.files('../data/asc/clm/crn', pattern = '.asc$', full.names = TRUE) %>% 
  grep(paste0(vrs, collapse = '|'), ., value = TRUE) %>% 
  mixedsort() %>% 
  stack()

# Mask 5 minutos
msk5 <- raster('//dapadfs/data_cluster_2/gcm/cmip5/downscaled/rcp60/global_5min/bcc_csm1_1/r1i1p1/2040_2069/bio_1') * 0 + 1
msk5 <- raster::crop(msk5, msk_shp) %>% raster::mask(., msk_shp)

# Removing the duplicated by cells ----------------------------------------
pnt_dup <- dup_cell(msk = msk, occ = pnt)[[1]]
pnt_dup5 <- dup_cell(msk = msk5, occ = pnt)[[1]]

# Extracting the values from th estack ------------------------------------
pnt_swd <- raster::extract(lyr, pnt_dup[,1:2]) %>% 
  cbind(pnt_dup, .) %>%
  as.tibble() %>% 
  drop_na()
pnt_swd %>% 
  pull(3) %>% 
  table
norm <- scores(pnt_swd[,4:ncol(pnt_swd)], 'z') ; norm_na <- norm
norm_na[abs(norm_na)>3.5]  <- NA 
normpoints <- cbind(pnt_swd[,1:3], norm_na) %>%
  tbl_df() %>% 
  drop_na()
nrow(pnt_swd)
nrow(normpoints)
write.csv(normpoints, '../data/pnt/tbl/run2/normpoints.csv', row.names = FALSE)

normpoints %>% 
  pull(3) %>% 
  table()

pnt_cln <- raster::extract(lyr, normpoints[,1:2]) %>% 
  cbind(normpoints[,1:3], .) %>%
  as.tibble() 
write.csv(pnt_cln, '../data/pnt/tbl/run2/occ_rmOtl.csv', row.names = FALSE)

pnt_cln %>% 
  pull(3) %>% 
  table()

# Clustering --------------------------------------------------------------
occ <- pnt_cln
env_values <- as.matrix(occ[,4:36]); nrow(env_values)
datRF <- as.data.frame(occ[,4:ncol(occ)]); nrow(datRF)
d <- dist(datRF, method = "euclidean")
rfClust <- rf.clust(occ = occ, nforest = 25, ntrees = 100, nVars = 8, nclasses = 5)
labelRF <- rfClust[[1]]
table(labelRF)
clusterdata <- rfClust[[2]]
classdata <- cbind(pb = as.factor(labelRF), occ[,3:ncol(occ)])
clusteredpresdata <- cbind(occ, cluster = labelRF) %>% na.omit() %>% tbl_df()
no.clusters <- 5

dir.create('../rds/run2', recursive = TRUE)
dir.create('../data/pnt/tbl/run2', recursive = TRUE)
run <- 'run2'

write.csv(clusteredpresdata, '../data/pnt/tbl/run2/clusteredpresdata.csv', row.names = FALSE)

saveRDS(object = datRF, file = paste0('../rds/', run, '/datRF.rds'))
saveRDS(object = clusterdata, file = paste0('../rds/', run, '/clusterdata.rds'))
saveRDS(object = env_values, file = paste0('../rds/', run, '/env_values.rds'))

# Cluster occurences
grouped <- as.data.frame(cbind(env_values, cluster = as.factor(labelRF))) %>% tbl_df()
occ_cluster <- cbind(occ[,1:2], grouped[,'cluster'])
saveRDS(object = occ_cluster, file = paste0('../rds/', run, '/occ_cluster.rds'))
saveRDS(object = grouped, file = paste0('../rds/', run, '/grouped.rds'))

# Cluster tables to shapefile cluster
nrow(occ_cluster)
head(occ_cluster)
coordinates(occ_cluster) <- ~ X + Y
dir.create('../data/pnt/shp/run2')
writeOGR(obj = occ_cluster, dsn = '../data/pnt/shp/run2', layer = 'occ_cluster', driver = 'ESRI Shapefile')







