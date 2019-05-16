
# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(raster, rgdal, rgeos, gtools, tidyverse, stringr, velox, sf, foreach, doSNOW, parallel, rgbif)

# Initial setup -----------------------------------------------------------
g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999)

# Load data ---------------------------------------------------------------
gbl <- read_csv('../data/pnt/tbl/cocoa_global.csv')
all <- shapefile('../data/shp/all_countries.shp')
pnt <- shapefile('../data/pnt/shp/IturiCocoaFields_geo.shp')
msk <- raster('../data/asc/bse/msk.asc')

# Download cocoa data from GBIF -------------------------------------------
gbf <- occ_data(scientificName = 'Theobroma cacao',
                limit = 20000,
                hasCoordinate = TRUE,
                hasGeospatialIssue = FALSE)[[2]]
coordinates(gbf) <- ~ decimalLongitude + decimalLatitude
crs(gbf) <- crs(pnt)
crs(msk) <- crs(pnt)

msk_pol <- rasterToPolygons(msk, dissolve = T)
writeOGR(obj = msk_pol, dsn = '../data/shp', layer = 'mask_polygon', driver = 'ESRI Shapefile')

countries <- crop(all, msk_pol)
writeOGR(obj = countries, dsn = '../data/shp', layer = 'countries_africa', driver = 'ESI Shapefile')
countries_vec <- countries@data$ENGLISH

# GBIF --------------------------------------------------------------------
gbf_zne <- gIntersection(msk_pol, gbf)
crd_gbf <- coordinates(gbf_zne) %>%
  as_tibble() %>% 
  mutate(source = 'gbif database')
plot(msk_pol)
points(crd_gbf$x, crd_gbf$y, pch = 16, col = 'red')

# Point fields ------------------------------------------------------------
pnt_fld <- coordinates(pnt) %>% 
  as_tibble() %>% 
  mutate(source = 'itury cocoa fields') %>% 
  setNames(colnames(crd_gbf))

# Data from cocoa global --------------------------------------------------
countries_global <- gbl %>% 
  distinct(Country) %>% 
  pull()
gbl_afr <- gbl %>% filter(Country %in% countries_vec) 

points(gbl_afr$Longitude, gbl_afr$Latitude, pch = 16, col = 'red')

pnt_afr <- gbl_afr %>% 
  dplyr::select(Longitude, Latitude) %>% 
  mutate(source = 'Global database') %>% 
  setNames(colnames(pnt_fld))

# Rbind -------------------------------------------------------------------
nrow(crd_gbf)
nrow(pnt_afr)
nrow(pnt_fld)
pnts <- rbind(crd_gbf, pnt_fld, pnt_afr)
write.csv(pnts, '../data/pnt/tbl/run1/all_points.csv', row.names = TRUE)

# Table to shapefile
dir.create('../data/shp/pnt/run1')
coordinates(pnts) <- ~ x + y
writeOGR(obj = pnts, dsn = '../data/shp/pnt/run1', layer = 'points_raw', driver = 'ESRI Shapefile')













