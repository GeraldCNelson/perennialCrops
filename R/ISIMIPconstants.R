# constants to be used in the various ISIMIP data crunching scripts The idea is to save some space in each of the scripts and to be a central place to look for what constants are available

#compression code for use with writeRaster
require(terra)
require(data.table)

woptList <- list(gdal=c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL = 6", "NUM_THREADS=ALL_CPUS"))

# choice variables, often used in for loops
# sspChoices <- c("ssp126", "ssp585") 
# modelChoices <- c( "GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL") 
# modelChoices_lower <- tolower(modelChoices)
# startYearChoices <-  c(2041, 2081) 
# hemispheres <- c("NH", "SH")
# extent_NH <- ext( -180, 180, 0, 90)
# extent_SH <-ext( -180, 180, -60, 0) #-60 gets rid of Antarctica for SH
# 
# yearRange <- 19

# file locations -----
locOfClimFiles <- "climdata/"
path_data <- "data/perennials/"
locOfCPFiles <- "data/chillPortions/chill_portions/"
locOfDataFiles_THI <- "data/THI/"
locOfResultsFiles <- "results/"
lofOfGraphicsFiles <- "graphics/"
locOfRawDataFiles <- "data-raw/"
locOfHarvestDataFiles <- "data-raw/crops/HarvestedAreaYield175Crops_Geotiff/GeoTiff/"
locOfgddsFiles <- "data/growingDegreeDays/"
path_data <- "data/runs/"

# general test values
modelChoice <-  "IPSL-CM6A-LR"
modelChoice_lower <- tolower(modelChoice)
k <- "ssp585"
l <- 2041
yearNumber <- l + 2
yearSpan <- paste0(l, "_", l + yearRange)
hem <- "SH"
