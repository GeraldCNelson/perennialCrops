source("R/perennialFunctions.R")

speciesNames <- gsub(varietiesChoiceSuffix, "", speciesChoices) 
# # to figure out what is wrong with olive
# speciesNames = "olive"
# source("R/perennialsPrep.R") # creates the data tables majorCropValues_main, majorCropValues_lo and majorCropValues_hi. Included in perennialFunctions.R
woptList <- list(gdal=c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL = 6", "NUM_THREADS=ALL_CPUS"))
terraOptions(memfrac = 0.9, copies = 1, progress = 10, tempdir="data/ISIMIP", verbose = FALSE) # memfrac = 2, 
#coastline data
coastline <- vect("data-raw/regionInformation/ne_50m_coastline/ne_50m_coastline.shp")
coastline_cropped <- crop(coastline, extent_noAntarctica )
coastline_cropped_Rob <- project(coastline_cropped, crsRob)
coastline_cropped_Rob_sf <- sf::st_as_sf(coastline_cropped_Rob)
colList <- c("#e66101","#fdae61", "#abd9e9", "#2c7bb6") # old
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") # from https://stackoverflow.com/questions/65013406/how-to-generate-30-distinct-colors-that-are-color-blind-friendly
colList <- c("#D55E00", "#56B4E9", "#009E73", "#F0E442")
yearRange <- 19
locOfClimFiles <- "climdata"

climVal <- -2 # changed from 0 because subject experts say plants can tolerate this
test_logic <- paste0("x > ", climVal)
logicDirection <- ">"
if (logicDirection == ">") ldtext <-"gt"
if (logicDirection == "<") ldtext <-"lt"
runlengthChoices <- c(100) # at least 100 days of tmin > 0
climateVariable <- "tasmin"
runsParms <- c(climVal, test_logic, logicDirection, ldtext, runlengthChoices, climateVariable)

  f_convert_graphics <- function(rawIn, desc) {
    # project spatRaster to Robinson and convert to data frame
    rawIn_rob <- project(rawIn, crsRob)
    rawIn_rob_df <- as.data.frame(rawIn_rob, xy = TRUE)
    names(rawIn_rob_df) <- c("x", "y", "value")
    rawIn_rob_df$type = desc
    rawIn_rob_df$value <- round(rawIn_rob_df$value, 0)
    rawIn_rob_df$value[rawIn_rob_df$value == 0] <- NA
    return(rawIn_rob_df)
  }

yearSpan_early <- "1991_2010"
yearSpan_mid <- "2041_2060"
yearSpan_end <- "2081_2100"
speciesChoices <- unique(cropVals$cropName)

