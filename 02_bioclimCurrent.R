
# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(raster, rgdal, rgeos, stringr, velox, sf, tidyverse, gtools, dismo)

# Initial setup -----------------------------------------------------------
g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999)
setwd('//mnt/workspace_cluster_9/Coffee_Cocoa2/_cocoaAfrica/r')
source('functionsR.R')

# Load data ---------------------------------------------------------------
msk <- raster('../data/asc/bse/msk.asc')

# Climate data ------------------------------------------------------------
fls <- list.files('../data/asc/clm/crn', full.names = TRUE, patter = '.asc') %>%
  mixedsort()
var <- c('prec', 'tmin', 'tmean', 'tmax')
for(i in 1:length(var)){
  eval(parse(text = paste0(var[i], ' <- grep(var[', i, '], fls, value = T) %>% mixedsort() %>% stack()')))
}

# Solar radiation ---------------------------------------------------------
xtr <- list.files('//mnt/workspace_cluster_9/Coffee_Cocoa2/_colombiaETP/_data/_tif/ET_SolRad', full.names = TRUE) %>%
  grep('/et_', ., value = TRUE) %>%
  mixedsort() %>% 
  stack() 
xtr <- raster::crop(xtr, extent(msk)) 
xtr <- resample(x = xtr, y = msk, method = 'ngb')
xtr <- xtr*c(31,29,31,30,31,30,31,31,30,31,30,31)
nms <- names(xtr)
xtr <- unstack(xtr)
Map('writeRaster', x = xtr, paste0('../data/asc/clm/crn/', nms, '.asc'), overwrite = TRUE)
xtr <- stack(xtr)

extent(xtr) == extent(msk)

# Calculating the ETP variables -------------------------------------------
etp <- 0.0013 * 0.408 * xtr * (tmean + 17) * (tmax - tmean - 0.0123 * prec) ^ 0.76
etp <- unstack(etp)
Map('writeRaster', x = etp, filename = paste0('../data/asc/clm/crn/etp_', 1:12, '.asc'), overwrite = TRUE)
etp <- stack(etp)

# To calculate bioclimatic variables --------------------------------------
biostack <- biovars(prec,tmin,tmax)
nms <- names(biostack)
biostack <- unstack(biostack)
Map('writeRaster', x = biostack, filename = paste0('../data/asc/clm/crn/', nms, '.asc'), overwrite = TRUE)

prc_mtx <- as.matrix(prec)
gcmbiofolder <- '../data/asc/clm/crn/'
bio20ies <-  t(apply(prc_mtx, 1, cumDry))
mm <- c(40,50,100)
for(i in 1:3){ 
  bio20 <- msk
  values(bio20) <- bio20ies[,i]
  writeRaster(bio20, paste(gcmbiofolder,"/","bio20_",mm[i],"mm",sep=""),
              format = "ascii",
              overwrite = T, 
              NAflag=-9999)
} 

# Deficit stack
dfc_stk <- prec - etp
DefAndTemp <- cbind(as.matrix(dfc_stk),as.matrix(tmean),as.matrix(tmax))
biovalues <-  t(apply(DefAndTemp,1,cumTemp))
biovalues <- round(biovalues, 0)
nms <- paste0('bio', 21:24)

lapply(1:ncol(biovalues), function(i){
  lyr <- msk
  values(lyr) <- biovalues[,i]
  writeRaster(lyr, paste0('../data/asc/clm/crn/', nms[i], '.asc'), overwrite = FALSE)
  print('Done!')
})

# ETP variables
ETPAndPrec <- cbind(as.matrix(etp),as.matrix(prec),as.matrix(tmean))
etpbios <-  t(apply(ETPAndPrec,1,etpvars))
etpbios <- round(etpbios, 0)
nms <- paste0('bio', 25:33)

lapply(1:ncol(etpbios), function(i){
  lyr <- msk  
  values(lyr) <- etpbios[,i]
  writeRaster(lyr, paste0('../data/asc/clm/crn/', nms[i], '.asc'), overwrite = TRUE)
  print('Done!')
})

# Review files
fls <- list.files('../data/asc/clm/crn', full.names = TRUE) %>% 
  grep('bio', ., value = T)
length(fls)

