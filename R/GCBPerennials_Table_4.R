# code to produce the GCB perennials paper Table 4.
library(terra)
library(data.table)
path_data <- "data/perennials/nonlimiting_all_"
path_harvestFiles <- "data-raw/crops/HarvestedAreaYield175Crops_Geotiff/GeoTiff/"

#source("R/perennialFunctions.R") # loads some functions, packages and some needed constants
speciesChoices <- c("almond_main", "apple_main",  "cherry_main", "olive_main",  "grape_main")
extent_noAntarctica <- ext(-180, 180, -60, 90)

f_getArea <- function(r, layer) {
  r_sub <- subset(r, layer)
  r_sub[r_sub < 1 | r_sub > 1] <- NA # keep only locations with value of 1
  r_sub_area <- (expanse(r_sub, unit = "km"))/1000 # convert to 1000 sq km
  r_sub_area <- round(r_sub_area, 0)
}

f_harvestArea <- function(speciesName, minArea) {
  rInArea <- rast(paste0(path_harvestFiles, speciesName,"/",  speciesName, "_HarvestedAreaHectares.tif"))
  rInArea <- crop(rInArea, extent_noAntarctica)
  harvestArea_earlyCent <- aggregate(rInArea, fact = 6, fun = "sum") # convert 5 arc minutes to 1/2 degrees
  maskMin <- switch(
    speciesName,
    "almond" = minArea,
    "apple" = minArea,
    "cherry" = minArea,
    "grape" = minArea,
    "olive" = minArea
  )
  harvestArea_earlyCent[harvestArea_earlyCent < maskMin] <- 0 # set minimum area to be greater than minarea hectares per grid cell
  return(harvestArea_earlyCent)
}

f_areacalcs <- function(speciesChoice) {
  speciesName <- gsub("_main", "", speciesChoice) 
  r_combined_hist <- rast(paste0(path, speciesChoice, "_", "historical", "_", "good", "_", "1991_2010", ".tif"))
  r_combined_ssp585_mid <- rast(paste0(path, speciesChoice, "_", "ssp585", "_", "good", "_", "2041_2060", ".tif"), lyrs = 1)
  r_combined_ssp585_end <- rast(paste0(path, speciesChoice, "_", "ssp585", "_", "good", "_", "2081_2100", ".tif"), lyrs = 1)
  suitableArea_historical <- r_combined_hist[[1]]
  
   # now get harvested area map
  harvestArea_earlyCent <- f_harvestArea(speciesName, minArea = 1)
  harvestArea_earlyCent[harvestArea_earlyCent > 0] <- 1
  harvestArea_earlyCent[harvestArea_earlyCent <= 0] <- NA 
  
  # area calculations -----
  commonArea <- harvestArea_earlyCent * r_combined_hist # - locations where all suitability requirements are met in recent period and the crop was harvested early century
  commonArea_area <- f_getArea(commonArea, 1) # area in common
  harvestArea_earlyCent_area <- f_getArea(harvestArea_earlyCent, 1) # harvested area early century (units are 1000 sq km from f_getArea)
  r_combined_hist_area <- f_getArea(r_combined_hist, 1)
  areas_combined <- list(speciesName, harvestArea_earlyCent_area, r_combined_hist_area, commonArea_area)
  return(areas_combined)
}

areaDat <- data.table(species = character(), earlyCenturyHarvestArea = numeric(), recentPeriodSuitability = numeric(), earlyCentCommon = numeric())
for (speciesChoice in speciesChoices) {
  print(speciesChoice)
  temp <- f_areacalcs(speciesChoice)
  areaDat <- rbind(temp, areaDat)
}
setorder(areaDat, cols = "species")

write.csv(areaDat, "results/earlyCenturyAreas.csv")

