
# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(raster, rgdal, rgeos, stringr, velox, sf, tidyverse, gtools, dismo)

# Initial setup -----------------------------------------------------------
g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999)

setwd('//mnt/workspace_cluster_9/Coffee_Cocoa2/_cocoaAfrica/r')
source('functionsR.R')

# Functions to use --------------------------------------------------------
makeBioclim <- function(gc, yr){
  # gc <- gcm[1]; yr <- yrs[1]
  gc <- paste0(gc, '/')
  fle <- paste0('../data/asc/clm/ftr/', yr, '/', gc) %>%
    list.files(., full.names = TRUE, pattern = '.asc$') %>%
    mixedsort()
  var <- c('prec', 'tmin', 'tmean', 'tmax')
  for(i in 1:length(var)){
    print(var[i])
    eval(parse(text = paste0(var[i], ' <- grep(var[', i, '], fle, value = T) %>% mixedsort() %>% stack()')))
  }
  biostack <- biovars(prec, tmin, tmax)
  nms <- names(biostack)
  biostack <- unstack(biostack)
  Map('writeRaster', x = biostack, filename = paste0('../data/asc/clm/ftr/', yr, '/', gc, '/bio', 1:19, '.asc'))
  
}
makeBioclim_2 <- function(gc, yr){
  # gc <- gcm[1]; yr <- yrs[1]
  gc <- paste0(gc, '/')
  fls <- paste0('../data/asc/clm/ftr/', yr, '/', gc) %>%
    list.files(., full.names = TRUE, pattern = '.asc$') %>% 
    mixedsort()
  var <- c('prec', 'tmin', 'tmean', 'tmax', 'etp')
  for(i in 1:length(var)){
    eval(parse(text = paste0(var[i], ' <- grep(var[', i, '], fls, value = T) %>% mixedsort() %>% stack()')))
  }
  msk <- prec[[1]] * 0 + 1
  prc_mtx <- as.matrix(prec)
  gcmbiofolder <- paste0('../data/asc/clm/ftr/', yr, '/', gc)
  bio20ies <- t(apply(prc_mtx, 1, cumDry))
  mm <- c(40, 50, 100)
  for(i in 1:3){
    bio20 <- msk
    values(bio20) <- bio20ies[,i]
    writeRaster(bio20, 
                filename = paste(gcmbiofolder, '/', 'bio20_', mm[i], "mm", sep = ''),
                format = 'ascii',
                overwrite = T,
                NAflag = -9999)
  }
  
  print('Deficit stacks')
  dfc_stk <- prec - etp
  DefAndTemp <- cbind(as.matrix(dfc_stk),as.matrix(tmean),as.matrix(tmax))
  biovalues <- t(apply(DefAndTemp,1,cumTemp))
  biovalues <- round(biovalues, 0)
  nms <- paste0('bio', 21:24)
  
  lapply(1:ncol(biovalues), function(i){
    lyr <- msk
    values(lyr) <- biovalues[,i]
    writeRaster(lyr, paste('../data/asc/clm/ftr/', yr, '/', gc, '/', nms[i], sep = ''), format = 'ascii', overwrite = T, NAflag = -9999)
  })
  
  ETPAndPrec <- cbind(as.matrix(etp),as.matrix(prec),as.matrix(tmean))
  etpbios <- t(apply(ETPAndPrec,1,etpvars))
  nms <- paste0('bio', 25:33)
  
  lapply(1:ncol(etpbios), function(i){
    print(nms[i])
    lyr <- msk
    values(lyr) <- etpbios[,i]
    writeRaster(lyr, paste('../data/asc/clm/ftr/', yr, '/', gc, '/', nms[i], sep = ''), format = 'ascii', overwrite = TRUE, NAflag = -9999)
  })

  print(paste0('Done ', gc)) 
  
}
calcETP <- function(gc, yr){
  # gc <- gcm[1]
  # yr <- yrs[1]
  fls <- paste0('../data/asc/clm/ftr/', yr, '/', gc) %>%
    list.files(., full.names = TRUE, pattern = '.asc$') %>%
    mixedsort()
  prec <- grep('prec', fls, value = TRUE) %>% stack()
  tmax <- grep('tmax', fls, value = TRUE) %>% stack()
  tmean <- grep('tmean', fls, value = TRUE) %>% stack()
  tmin <- grep('tmin', fls, value = TRUE) %>% stack()
  etp <- 0.00013 * 0.408 * xtr * (tmean + 17 ) * (tmax - tmean - 0.0123 * prec) ^ 0.76
  etp <- round(etp, 0)
  etp <- unstack(etp)
  Map('writeRaster', x = etp, filename = paste0('../data/asc/clm/ftr/', yr, '/', gc, '/etp_', 1:12, '.asc'), overwrite = TRUE)
  print('Done!')
}

# Load data ---------------------------------------------------------------
msk <- raster('../data/asc/bse/msk.asc')
gcm <- list.files('../data/asc/clm/ftr/2020_2049')
yrs <- list.files('../data/asc/clm/ftr')
xtr <- list.files('../data/asc/clm/crn', full.names = TRUE, pattern = 'et_sol') %>% 
  mixedsort() %>% 
  stack()

# Apply the functions -----------------------------------------------------

# Bioclimatic 1:19
registerDoMC(19)
foreach(i = 1:length(gcm), .packages = c('raster', 'rgdal', 'tidyverse', 'dismo', 'rgeos', 'gtools'), .verbose = TRUE) %dopar% {
  makeBioclim(gc = gcm[i], yr = yrs[2])
}

# ETP variables
foreach(i = 1:length(gcm), .packages = c('raster', 'rgdal', 'dplyr', 'gtools', 'foreach', 'sp', 'stringr', 'dismo'), .verbose = TRUE) %dopar% {
  calcETP(gc = gcm[i], yr = yrs[2])
}

# Biocliamatic 20:33
foreach(i = 1:length(gcm), .packages = c('raster', 'rgdal', 'dplyr', 'gtools', 'foreach', 'sp', 'stringr', 'dismo'), .verbose = TRUE) %dopar% {
  makeBioclim_2(gc = gcm[i], yr = yrs[2])
}






