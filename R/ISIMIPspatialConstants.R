# ISIMIP project spatial constants
require(terra)
require(data.table)
source("R/ISIMIPconstants.R")
extent_noAntarctica <- ext(-180, 180, -60, 90) #-60 gets rid of Antarctica for global
northernHemExtent <- ext( -180, 180, 0, 90)
southernHemExtent <-ext( -180, 180, -60, 0)

woptList <- list(gdal=c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL = 6", "NUM_THREADS=ALL_CPUS"))

# projection choices -----
RobinsonProj <-  "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
# RobinsonProj <- "epsg:54030"

GoodeHomolosineProj <- "+proj=goode +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" # see https://proj.org/operations/projections/goode.html
crsRob <- RobinsonProj
crsGoode <- GoodeHomolosineProj
crslatlong <- "+proj=longlat +datum=WGS84 +no_defs"

#coastline data
coastline <- vect("data-raw/regionInformation/ne_50m_coastline/ne_50m_coastline.shp")
coastline_cropped <- crop(coastline, extent_noAntarctica )
coastline_cropped_Rob <- project(coastline_cropped, RobinsonProj)
coastline_cropped_igh <- project(coastline_cropped, GoodeHomolosineProj)

coastline_cropped_Rob_sf <- sf::st_as_sf(coastline_cropped_Rob)

#landOnlyMask <- rast(paste0(locOfRawDataFiles, "landseamask.nc")) # 0 is water; 1 is land
landOnlyMaskNoAntarctica <- rast(paste0(locOfRawDataFiles, "landseamask_no_antarctica.nc"))
landOnlyMaskNoAntarctica <- crop(landOnlyMaskNoAntarctica, extent_noAntarctica)
# Note: the landOnlyMaskNoAntarctica uses 1 for land and 0 for ocean. To use it as a mask do something like the following

# degree symbol - "°C"