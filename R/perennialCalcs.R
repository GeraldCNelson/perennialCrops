# code to do various perennial crop calculations
# varietiesChoiceSuffix <- "_main" # whether to use lo, main, or hi cp values. Needs to be before the sourcing code next. Moved to perennialFunctions.R
source("R/perennialFunctions.R")

# extreme cold calcs, scenarios -----
for (k in sspChoices) {
  for (l in startYearChoices) {
    for (hem in hemispheres) {
      for (speciesChoice in speciesChoices) {
        f_extremeCold(k, l, speciesChoice, hem, varietiesChoiceSuffix)
      }
    }
  }
}

# extreme cold calcs, historical -----
k <- "historical"
l <- 1991
for (hem in hemispheres) {
  for (speciesChoice in speciesChoices){
    f_extremeCold(k, l, speciesChoice, hem, varietiesChoiceSuffix)
  }
}

# heat damage, historical -----
k <- "historical"
l <- 1991
for (hem in hemispheres) {
  for (speciesChoice in speciesChoices) {
    f_heatDamage(k, l, speciesChoice, hem, varietiesChoice, suitabilityLevel)
  }
}

# heat damage, scenarios -----
for (k in sspChoices) {
  for (l in startYearChoices) {
    for (hem in hemispheres) {
      for (speciesChoice in speciesChoices) {
        f_heatDamage(k, l, speciesChoice, hem, varietiesChoice, suitabilityLevel)
      }
    }
  }
}

# frost damage, scenarios -----
for (k in sspChoices) {
  for (l in startYearChoices) {
    for (hem in hemispheres) {
      for (speciesChoice in speciesChoices) {
        f_frostDamage(k, l, speciesChoice, hem, varietiesChoice, suitabilityLevel)
      }
    }
  }
}

# frost damage, historical -----
k <- "historical"
l <- 1991
for (hem in hemispheres) {
  for (speciesChoice in speciesChoices) {
    f_frostDamage(k, l, speciesChoice, hem, varietiesChoice, suitabilityLevel)
  }
}
# 
# # gdd calcs ----
# gddFilesCompleted <- list.files(locOfgddsFiles,  full.names = TRUE)
# gddFilesCompleted <- gddFilesCompleted[!grepl("aux.xml", gddFilesCompleted, fixed = TRUE)]
# gddFilesCompleted <- gsub("//", "/", gddFilesCompleted)
# gddFilesCompleted <- NULL # temporary to make sure everything is redone
# gdds, historical -----
k <- "historical"
l <- 1991
for (modelChoice in modelChoices) {
  for (hem in hemispheres) {
    print(system.time(f_computeGDDs(k, l, modelChoice, cropVals)))
  }
}

# gdds, scenarios -----
for (k in sspChoices) {
  for (l in startYearChoices) {
    for (modelChoice in modelChoices) {
      for (hem in hemispheres) {
        f_computeGDDs(k, l, modelChoice, cropVals) 
      }
    }
  }
}
# just almond and winegrapes here because the gdds are the same for them and the other crops. this is true for both main and lo chill portions varieties
# gdd sums, historical ------
k <- "historical"
l <- 1991
for (speciesChoice in c(paste0("almond", varietiesChoiceSuffix), paste0("winegrape", varietiesChoiceSuffix))) {
  for (hem in hemispheres) {
    f_gddSums(k, l, speciesChoice, hem, runsParms)
  }
}

# gdd sums, scenarios ------
for (speciesChoice in c(paste0("almond", varietiesChoiceSuffix), paste0("winegrape", varietiesChoiceSuffix))) { #
  for (k in sspChoices) {
    for (l in startYearChoices) {
      for (hem in hemispheres) {
        f_gddSums(k, l, speciesChoice, hem, runsParms)
      }
    }
  }
}

# added because almond, apple, and cherry have the same gdd requirements and winegrape and olive also have the same
# first get all the gdd file names and pull out just the ones for almond and winegrape
gddFilesCompleted <- list.files(locOfgddsFiles,  full.names = TRUE)
gddFilesCompleted <- gddFilesCompleted[!grepl("aux.xml", gddFilesCompleted, fixed = TRUE)]
gddFilesCompleted <- gsub("//", "/", gddFilesCompleted)
gddSumFilesCompleted <- gddFilesCompleted[!grepl("gddSum_mean_", gddFilesCompleted, fixed = TRUE)]
gddSumFilesCompleted <- gddSumFilesCompleted[!grepl("ensemble_", gddSumFilesCompleted, fixed = TRUE)]
gddSumFilesCompleted <- gddSumFilesCompleted[grepl("gddSum_", gddSumFilesCompleted, fixed = TRUE)]

gddFiles_almond <- gddSumFilesCompleted[grepl("almond", gddSumFilesCompleted, fixed = TRUE)]
gddFiles_winegrape <- gddSumFilesCompleted[grepl("winegrape", gddSumFilesCompleted, fixed = TRUE)]

# because the gdds are the same for main and lo, the next bit of code just creates a new  filename and copies the main to the new

for (i in gddFiles_almond) {
  gdd_cherry <- gsub("almond", "cherry", i)
  
  file.copy(from = i, to = gdd_cherry, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
  gdd_apple <- gsub("almond", "apple", i)
  file.copy(from = i, to = gdd_apple, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
  print(gdd_cherry)
  print(i)
  al <- rast(i)
  plot(al, 1)
  al <- rast(gdd_cherry)
  plot(al, 1)
  
}
for (i in gddFiles_winegrape) {
  gdd_olive <- gsub("winegrape", "olive", i)
  file.copy(from = i, to = gdd_olive, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
}

for ( i in gddSumFilesCompleted) {
  i_lo <- gsub("_main", "_lo", i) 
  file.copy(from = i, to = i_lo, overwrite = TRUE, copy.mode = TRUE)
}

# GDD sum, means by model, historical -----
k <- "historical"
l <- 1991
f_gddSum_mean(k, l, speciesChoices, modelChoices, hemispheres)

# GDD sum, means by model, scenarios -----
for (k in sspChoices) {
  for (l in startYearChoices) {
    f_gddSum_mean(k, l, speciesChoices, modelChoices, hemispheres)
  }
}

# ensemble GDD sum, historical -----
k <- "historical"
l <- 1991
f_ensemble_GDD_sum_mean(k, l, yearRange, speciesChoices, cropVals)

# ensemble GDD sum, scenarios -----
for (k in sspChoices) {
  for (l in startYearChoices) {
    yearSpan <- paste0(l, "_", l + yearRange)
    f_ensemble_GDD_sum_mean(k, l, yearRange, speciesChoices, cropVals)
  }
}

# gdds suitable, scenarios -----
for (k in sspChoices) {
  for (l in startYearChoices) {
    for (hem in hemispheres) {
      for (speciesChoice in speciesChoices) {
        print(paste0("speciesChoice: ", speciesChoice, ", ssp choice: ", k, ", start year: ", l))
        f_gddsSuitability(k, l, speciesChoice, hem) 
      }
    }
  }
}

# gdds suitable, historical -----
k <- "historical"
l <- 1991
for (hem in hemispheres) {
  for (speciesChoice in speciesChoices) {
    print(paste0("speciesChoice: ", speciesChoice, ", ssp choice: ", k, ", start year: ", l))
    f_gddsSuitability(k, l, speciesChoice, hem) 
  }
}

# combined damage, scenarios -----
# code to read in 1/0 metrics files and produce 1/0 tifs where the crop is potentially growable. The chill portions files are created in the chillPortions.R script
# combined damage, historical -----
k <- "historical"
l <- 1991
for (speciesChoice in speciesChoices) {
  print(paste0("speciesChoice: ", speciesChoice, ", ssp choice: ", k, ", start year: ", l))
  f_combinedDamage(k, l, speciesChoice, suitabilityLevel) 
} 

for (k in sspChoices) {
  for (l in startYearChoices) {
    for (speciesChoice in speciesChoices) {
      print(paste0("speciesChoice: ", speciesChoice, ", ssp choice: ", k, ", start year: ", l))
      f_combinedDamage(k, l, speciesChoice, suitabilityLevel) 
    }
  }
}

# area calculations -----
{
  dt_area <- data.table(species = character(), cultivar = character(), chillPortions= numeric(), hemisphere = character(), ssp = character(), yearSpan = character(), quality = character(), area_suitable = numeric(), rasterName = character())
  suitabilityLevels <- "good"
  for (k in sspChoices) {
    for (l in startYearChoices) {
      yearSpan <- paste0(l, "_", l + yearRange)
      for (speciesChoice in speciesChoices) {
        CPfruit <- cropVals[cropName == speciesChoice, chill_portions]
        cultivar <-  cropVals[cropName == speciesChoice, cultivar]
        
        for (hem in hemispheres) {
          for (suitabilityLevel in suitabilityLevels) {
            fileName_in <- paste0(locOfDataFiles_perennials, "nonlimiting_all_", speciesChoice, "_", k, "_", suitabilityLevel, "_", hem, "_", yearSpan, ".tif")
            rastName <- paste0("suitable_", speciesChoice, "_", hem, "_", k, "_", suitabilityLevel, "_", yearSpan)
            r <- rast(fileName_in)
            # r has all 5 suitability layers, get just the combined one
            r_combSuit <- r$combinedSuit
            r_combSuit[r_combSuit == 0] <- NA # only want to get area of cells that have a value of 1; ie, suitable
            r_area <- expanse(r_combSuit, unit = "km")
            print(r_area)
            dt_area <- rbind(dt_area, list(speciesChoice, cultivar, CPfruit, hem, k, yearSpan, suitabilityLevel, r_area, rastName))
          }
        }
      }
    }
  }
  
  k <- "historical"
  l <- 1991
  yearSpan <- paste0(l, "_", l + yearRange)
  for (speciesChoice in speciesChoices) {
    CPfruit <- cropVals[cropName == speciesChoice, chill_portions]
    cultivar <-  cropVals[cropName == speciesChoice, cultivar]
    for (hem in hemispheres) {
      for (suitabilityLevel in suitabilityLevels) {
        fileName_in <- paste0(locOfDataFiles_perennials, "nonlimiting_all_", speciesChoice, "_", k, "_", suitabilityLevel, "_", hem, "_", yearSpan, ".tif")
        r <- rast(fileName_in) 
        # r has all 5 suitability layers, get just the combined one
        r_combinedSuit <- r$combinedSuit
        r_combinedSuit[r_combinedSuit == 0] <- NA # only want to get area of cells that have a value of 1; ie, suitable
        r_area <- expanse(r_combinedSuit, unit = "km")
        dt_area <- rbind(dt_area, list(speciesChoice, cultivar, CPfruit, hem, k, yearSpan, suitabilityLevel, r_area, rastName))
      }
    }
  }
  fileName_out <- paste0(locOfDataFiles_perennials, "areaCalcs", varietiesChoiceSuffix, ".csv")
  write.csv(dt_area, file = fileName_out, row.names = FALSE)
  print(paste0("fileName out: ", fileName_out))
}

# area deltas, need to run code above first -----
dt_area_delta <- data.table(species = character(), cultivar = character(), chillPortions= numeric(), hemisphere = character(), ssp_base = character(), ssp = character(), yearSpan = character(), quality = character(), area_base = numeric(), area_delta = numeric(), delta_share = numeric())
for (speciesChoice in speciesChoices) {
  CPfruit <- cropVals[cropName == speciesChoice, chill_portions]
  cultivar <-  cropVals[cropName == speciesChoice, cultivar]
  for (hem in hemispheres) {
    for (suitabilityLevel in suitabilityLevels) {
      r_historical <- rast(paste0(locOfDataFiles_perennials, "nonlimiting_all_", speciesChoice, "_", "historical", "_", suitabilityLevel, "_", hem, "_", "1991_2010.tif"))
      r_historical[r_historical == 0] <- NA
      area_base <- expanse(r_historical, unit = "km")
      for (l in startYearChoices) {
        yearSpan <- paste0(l, "_", l + yearRange)
        for (k in sspChoices) {
          rastName <- paste0(locOfDataFiles_perennials, "nonlimiting_all_", speciesChoice, "_",  k, "_", suitabilityLevel, "_", hem, "_", yearSpan, ".tif")
          r_delta <- rast(rastName) - r_historical
          r_delta[r_delta == 0] <- NA
          r_delta_area <- expanse(r_delta, unit = "km")
          delta_ratio <- 100 * r_delta_area/area_base # added 100 * here to convert to percent November 9, 2021
          dt_area_delta <- rbind(dt_area_delta, list(speciesChoice, cultivar, CPfruit, hem, "historical", k, yearSpan, suitabilityLevel, area_base, r_delta_area, delta_ratio))
        }
      }
    }
  }
}

fileName_out <- paste0(locOfDataFiles_perennials, "areaCalcs_delta", varietiesChoiceSuffix, ".csv")
write.csv(dt_area_delta, file = fileName_out, row.names = FALSE)
print(paste0("fileName out: ", fileName_out))

# harvest and suitability area calcs ------
# area common to mid and end century -----
dt_area_common_mid <- data.table(species = character(), cultivar = character(), chillPortions= numeric(), ssp = character(), hemisphere = character(), yearSpan = character(), area_common = numeric())
dt_area_common_end <- data.table(species = character(), cultivar = character(), chillPortions= numeric(), ssp = character(), hemisphere = character(), yearSpan = character(), area_common = numeric())

for (speciesChoice in speciesChoices) {
  CPfruit <- cropVals[cropName == speciesChoice, chill_portions]
  cultivar <-  cropVals[cropName == speciesChoice, cultivar]
  speciesName <- gsub(varietiesChoiceSuffix, "", speciesChoice) # needed for the harvested area data
  # the harvested area data is just for grapes so need to get rid of wine in the names
  if (speciesChoice == paste0("winegrape", varietiesChoiceSuffix)) speciesName <- "grape"
  fileName_in <- paste0(locOfHarvestDataFiles, speciesName,"/",  speciesName, "_HarvestedAreaHectares.tif")
  rInArea <- rast(fileName_in)
  harvestArea_earlyCent <- aggregate(rInArea, fact = 6, fun = "sum") # convert 5 arc minutes to 1/2 degrees
  # crop to no Antarctica
  harvestArea_earlyCent <- crop(harvestArea_earlyCent, extent_noAntarctica)
  # mask to land only
  harvestArea_earlyCent <- mask(harvestArea_earlyCent, landOnlyMaskNoAntarctica, maskvalues = 0)
  maskMin <- switch(
    speciesName,
    "almond" = 10,
    "apple" = 10,
    "cherry" = 10,
    "grape" = 10,
    "olive" = 10
  )
  #  harvestArea_earlyCent[harvestArea_earlyCent < maskMin] <- NA # set minimum area to be greater than maxmin hectares per grid cell
  harvestArea_earlyCent[harvestArea_earlyCent < maskMin] <- 0 # set minimum area to be greater than maxmin hectares per grid cell
  harvestArea_earlyCent[harvestArea_earlyCent > 0] <- 1
  harvestArea_earlyCent_NH <- crop(harvestArea_earlyCent, extent_NH)
  harvestArea_earlyCent_SH <- crop(harvestArea_earlyCent, extent_SH)
  
  #suitable areas -----
  # suitable area end century
  fileName_in_ssp126_SH_end <- paste0(locOfDataFiles_perennials, "nonlimiting_all_", speciesChoice, "_", "ssp126", "_", suitabilityLevel, "_", "SH", "_", "2081_2100", ".tif")
  fileName_in_ssp126_NH_end <- paste0(locOfDataFiles_perennials, "nonlimiting_all_", speciesChoice, "_", "ssp126", "_", suitabilityLevel, "_", "NH", "_", "2081_2100", ".tif")
  fileName_in_ssp585_SH_end <- paste0(locOfDataFiles_perennials, "nonlimiting_all_", speciesChoice, "_", "ssp585", "_", suitabilityLevel, "_", "SH", "_", "2081_2100", ".tif")
  fileName_in_ssp585_NH_end <- paste0(locOfDataFiles_perennials, "nonlimiting_all_", speciesChoice, "_", "ssp585", "_", suitabilityLevel, "_", "NH", "_", "2081_2100", ".tif")
  
  suitableArea_ssp126_SH_end <- rast(fileName_in_ssp126_SH_end)
  suitableArea_ssp126_NH_end <- rast(fileName_in_ssp126_NH_end)
  suitableArea_ssp585_SH_end <- rast(fileName_in_ssp585_SH_end)
  suitableArea_ssp585_NH_end <- rast(fileName_in_ssp585_NH_end)
  
  # suitable area mid century
  fileName_in_ssp126_SH_mid <- paste0(locOfDataFiles_perennials, "nonlimiting_all_", speciesChoice, "_", "ssp126", "_", suitabilityLevel, "_", "SH", "_", "2041_2060", ".tif")
  fileName_in_ssp126_NH_mid <- paste0(locOfDataFiles_perennials, "nonlimiting_all_", speciesChoice, "_", "ssp126", "_", suitabilityLevel, "_", "NH", "_", "2041_2060", ".tif")
  fileName_in_ssp585_SH_mid <- paste0(locOfDataFiles_perennials, "nonlimiting_all_", speciesChoice, "_", "ssp585", "_", suitabilityLevel, "_", "SH", "_", "2041_2060", ".tif")
  fileName_in_ssp585_NH_mid <- paste0(locOfDataFiles_perennials, "nonlimiting_all_", speciesChoice, "_", "ssp585", "_", suitabilityLevel, "_", "NH", "_", "2041_2060", ".tif")
  
  suitableArea_ssp126_SH_mid <- rast(fileName_in_ssp126_SH_mid)
  suitableArea_ssp126_NH_mid <- rast(fileName_in_ssp126_NH_mid)
  suitableArea_ssp585_SH_mid <- rast(fileName_in_ssp585_SH_mid)
  suitableArea_ssp585_NH_mid <- rast(fileName_in_ssp585_NH_mid)
  
  # get just the overall suitability
  suitableArea_ssp126_SH_mid <- suitableArea_ssp126_SH_mid$combinedSuit
  suitableArea_ssp126_NH_mid <- suitableArea_ssp126_NH_mid$combinedSuit
  suitableArea_ssp585_SH_mid <- suitableArea_ssp585_SH_mid$combinedSuit
  suitableArea_ssp585_NH_mid <- suitableArea_ssp585_NH_mid$combinedSuit
  
  suitableArea_ssp126_SH_end <- suitableArea_ssp126_SH_end$combinedSuit
  suitableArea_ssp126_NH_end <- suitableArea_ssp126_NH_end$combinedSuit
  suitableArea_ssp585_SH_end <- suitableArea_ssp585_SH_end$combinedSuit
  suitableArea_ssp585_NH_end <- suitableArea_ssp585_NH_end$combinedSuit
  
  #combined
  suitableArea_ssp126_mid <- merge(suitableArea_ssp126_SH_mid, suitableArea_ssp126_NH_mid)
  suitableArea_ssp585_mid <- merge(suitableArea_ssp585_SH_mid, suitableArea_ssp585_NH_mid)
  
  suitableArea_ssp126_end <- merge(suitableArea_ssp126_SH_end, suitableArea_ssp126_NH_end)
  suitableArea_ssp585_end <- merge(suitableArea_ssp585_SH_end, suitableArea_ssp585_NH_end)
  
  r_combined_ssp126_mid <- c(harvestArea_earlyCent, suitableArea_ssp126_mid)
  r_combined_ssp585_mid <- c(harvestArea_earlyCent, suitableArea_ssp585_mid)
  
  commonArea_ssp126_NH_mid <- suitableArea_ssp126_NH_mid * harvestArea_earlyCent_NH # changed to * June 9, 2021
  commonArea_ssp126_SH_mid <- suitableArea_ssp126_SH_mid * harvestArea_earlyCent_SH
  commonArea_ssp585_NH_mid <- suitableArea_ssp585_NH_mid * harvestArea_earlyCent_NH # changed to * June 9, 2021
  commonArea_ssp585_SH_mid <- suitableArea_ssp585_SH_mid * harvestArea_earlyCent_SH
  commonArea_ssp126_NH_mid[commonArea_ssp126_NH_mid == 0] <- NA
  commonArea_ssp126_SH_mid[commonArea_ssp126_SH_mid == 0] <- NA
  commonArea_ssp585_NH_mid[commonArea_ssp585_NH_mid == 0] <- NA
  commonArea_ssp585_SH_mid[commonArea_ssp585_SH_mid == 0] <- NA
  
  commonArea_ssp126_NH_end <- suitableArea_ssp126_NH_end * harvestArea_earlyCent_NH # changed to * June 9, 2021
  commonArea_ssp126_SH_end <- suitableArea_ssp126_SH_end * harvestArea_earlyCent_SH
  commonArea_ssp585_NH_end <- suitableArea_ssp585_NH_end * harvestArea_earlyCent_NH # changed to * June 9, 2021
  commonArea_ssp585_SH_end <- suitableArea_ssp585_SH_end * harvestArea_earlyCent_SH
  commonArea_ssp126_NH_end[commonArea_ssp126_NH_end == 0] <- NA
  commonArea_ssp126_SH_end[commonArea_ssp126_SH_end == 0] <- NA
  commonArea_ssp585_NH_end[commonArea_ssp585_NH_end == 0] <- NA
  commonArea_ssp585_SH_end[commonArea_ssp585_SH_end == 0] <- NA
  
  
  cellSize(commonArea_ssp126_NH_mid, mask=TRUE, unit="km")
  cellSize(commonArea_ssp126_SH_mid, mask=TRUE, unit="km")
  cellSize(commonArea_ssp585_NH_mid, mask=TRUE, unit="km")
  cellSize(commonArea_ssp585_SH_mid, mask=TRUE, unit="km")
  
  cellSize(commonArea_ssp126_NH_end, mask=TRUE, unit="km")
  cellSize(commonArea_ssp126_SH_end, mask=TRUE, unit="km")
  cellSize(commonArea_ssp585_NH_end, mask=TRUE, unit="km")
  cellSize(commonArea_ssp585_SH_end, mask=TRUE, unit="km")
  
  dt_area_common_mid <- rbind(dt_area_common_mid, list(speciesChoice, cultivar, CPfruit, "ssp126", "NH", "2041_2060", expanse(commonArea_ssp126_NH_mid, unit = "km")))
  dt_area_common_mid <- rbind(dt_area_common_mid, list(speciesChoice, cultivar, CPfruit, "ssp126", "SH", "2041_2060", expanse(commonArea_ssp126_SH_mid, unit = "km")))
  dt_area_common_mid <- rbind(dt_area_common_mid, list(speciesChoice, cultivar, CPfruit, "ssp585", "NH", "2041_2060", expanse(commonArea_ssp585_NH_mid, unit = "km")))
  dt_area_common_mid <- rbind(dt_area_common_mid, list(speciesChoice, cultivar, CPfruit, "ssp585", "SH", "2041_2060", expanse(commonArea_ssp585_SH_mid, unit = "km")))
  
  dt_area_common_end <- rbind(dt_area_common_end, list(speciesChoice, cultivar, CPfruit, "ssp126", "NH", "2081_2100", expanse(commonArea_ssp126_NH_end, unit = "km")))
  dt_area_common_end <- rbind(dt_area_common_end, list(speciesChoice, cultivar, CPfruit, "ssp126", "SH", "2081_2100", expanse(commonArea_ssp126_SH_end, unit = "km")))
  dt_area_common_end <- rbind(dt_area_common_end, list(speciesChoice, cultivar, CPfruit, "ssp585", "NH", "2081_2100", expanse(commonArea_ssp585_NH_end, unit = "km")))
  dt_area_common_end <- rbind(dt_area_common_end, list(speciesChoice, cultivar, CPfruit, "ssp585", "SH", "2081_2100", expanse(commonArea_ssp585_SH_end, unit = "km")))
} 

dt_area_common <- rbind(dt_area_common_mid, dt_area_common_end)
dt_area_common[, ssp_year := paste0(ssp, "_", yearSpan)]
dt_area_common[, yearSpan := NULL]
# fileName_out <- paste0(locOfDataFiles_perennials, "areaCalcs_common", "_", k, varietiesChoiceSuffix, ".csv")
# write.csv(dt_area_common, file = fileName_out, row.names = FALSE)
# print(paste0("fileName out: ", fileName_out))

# create table of area changes ------
dt_area <- as.data.table(read.csv(file = paste0(locOfDataFiles_perennials, "areaCalcs", varietiesChoiceSuffix, ".csv")))

# delete rows that are for suitability bad or acceptable and mid century - keeping just early and mid century. And mid century, Oct 28
#dt_area <- dt_area[!quality %in% c("bad", "acceptable") & !yearSpan %in% "2041_2060",]
dt_area <- dt_area[!quality %in% c("bad", "acceptable") ,]
dt_area[,ssp_year := paste0(ssp, "_", yearSpan)]
dt_area[, c("quality", "rasterName", "ssp", "yearSpan") := NULL]
dt_area_wide <-        dcast(dt_area,        species + cultivar + chillPortions + hemisphere ~ ssp_year, value.var = "area_suitable")
dt_area_common_wide <- dcast(dt_area_common, species + cultivar + chillPortions + hemisphere ~ ssp_year, value.var = "area_common")
sspNames <- c("ssp126_2041_2060", "ssp126_2081_2100", "ssp585_2041_2060", "ssp585_2081_2100")
sspNames_new = paste0(sspNames, "_common")
setnames(dt_area_common_wide, old = sspNames, new =sspNames_new)
# dt_area_common[, ssp := NULL]
combined <- merge(dt_area_wide, dt_area_common_wide, all = TRUE)

# convert units
colsToConvert <- c("historical_1991_2010", sspNames, sspNames_new)
combined[, (colsToConvert) := lapply(.SD, '/', 1000), .SDcols = (colsToConvert)] # convert to 1000 sq km
# added 100 * below to convert to percent Nov 9, 2021
combined[, ratiomid2Early_126 := 100 * (-1 + ssp126_2041_2060/historical_1991_2010)]
combined[, ratiomid2Early_585 := 100 * (-1 + ssp585_2041_2060/historical_1991_2010)]
combined[, ratioCommon2Early_126 := 100 * (ssp126_2041_2060_common/historical_1991_2010)]
combined[, ratioCommon2Early_585 := 100 * (ssp585_2041_2060_common/historical_1991_2010)]
combined[, ratioLossOfEarly_126 := 100 * ((historical_1991_2010 - ssp126_2041_2060)/historical_1991_2010)]
combined[, ratioLossOfEarly_585 := 100 * ((historical_1991_2010 - ssp585_2041_2060)/historical_1991_2010)]
ratioColumns <- c("ratiomid2Early_126", "ratiomid2Early_585", "ratioCommon2Early_126", "ratioCommon2Early_585", "ratioLossOfEarly_126", "ratioLossOfEarly_585")
sumColumns <- c("historical_1991_2010", "ssp126_2041_2060", "ssp585_2041_2060", "ssp126_2041_2060_common", "ssp585_2041_2060_common", "ssp126_2081_2100", "ssp585_2081_2100", "ssp126_2081_2100_common", "ssp585_2081_2100_common")
ssp126Columns <- c("ssp126_2041_2060", "ssp126_2041_2060_common", "ratiomid2Early_126", "ratioCommon2Early_126", "ratioLossOfEarly_126") 
combined[, (ratioColumns):= round(.SD, 2), .SDcols = ratioColumns]
combined[, (sumColumns):= round(.SD, 0), .SDcols = sumColumns]
combined[, chillPortions := NULL] # column not needed for presentation table

write.csv(combined, file = paste0(locOfDataFiles_perennials, "sumTable", varietiesChoiceSuffix, ".csv"), row.names = FALSE)

# prepare summary table ------
# main varieties -----
library(flextable)
library(officer)
sumTable_main <- as.data.table(read.csv(file = paste0(locOfDataFiles_perennials, "sumTable", "_main", ".csv")))
sumTable_main[, species := gsub("_main", "", species)]
sumTable_main[, species := stringr::str_to_title(species)][species == "Winegrape", species := "Wine grape"]

# sumTable_main[, c(ssp126Columns, "ratioCommon2Early_585", "ratioLossOfEarly_585") := NULL]
namesOrderOrig <- names(sumTable_main)
namesOrderNew <- c("species", "cultivar", "hemisphere", "historical_1991_2010", "ssp126_2041_2060", "ssp126_2081_2100", "ssp585_2041_2060", "ssp585_2081_2100", 
                   "ssp126_2041_2060_common", "ssp585_2041_2060_common", 
                   "ssp126_2081_2100_common", "ssp585_2081_2100_common", 
                   "ratiomid2Early_126", "ratiomid2Early_585", 
                   "ratioCommon2Early_126", "ratioCommon2Early_585", 
                   "ratioLossOfEarly_126", "ratioLossOfEarly_585")
setcolorder(sumTable_main, namesOrderNew)
sumTable_main_flex <- flextable(sumTable_main)
# sumTable_main_flex <- set_header_labels(sumTable_main_flex, values = list(
#   species = "species", cultivar = "cultivar", hemisphere = "hemisphere", 
#   historical = "early century area", ssp585 = "mid century area, SSP585", ssp585_common = "common area, SSP mid century and early century", ratiomid2Early_585 = "area ratio, mid spp585 to beginning"))
typology_all <- data.frame(
  col_keys_all = names(sumTable_main),
  what = c("Species", "Cultivar", "Hemi-\nsphere", #3
           "Area (sq. km)", "Area (sq. km)","Area (sq. km)","Area (sq. km)","Area (sq. km)","Area (sq. km)","Area (sq. km)","Area (sq. km)","Area (sq. km)", #9
           "Area ratio", "Area ratio", "Area ratio", "Area ratio", "Area ratio", "Area ratio"), #15
  measure = c("Species", "Cultivar", "Hemi-\nsphere", #3
              "Early century", "mid century, SSP1-2.6", "end century, SSP1-2.6", "mid century, SSP5-8.5", "end century, SSP5-8.5", #5
              "mid century, SSP1-2.6, common", "end century, SSP1-2.6, common", "mid century, SSP5-8.5, common", "end century, SSP5-8.5, common", #3
              "Mid to early century, SSP1-2.6", "Mid to early century, SSP5-8.5", "Common to early century, SSP1-2.6", "Common to early century, SSP5-8.5", "Loss of early to 126", "Loss of early to 585"), # 6
  stringsAsFactors = FALSE )

sumTable_main_flex <- set_header_df(sumTable_main_flex, mapping = typology_all, key = "col_keys_all")
sumTable_main_flex <- merge_h(sumTable_main_flex, part = "header")
sumTable_main_flex <- merge_v(sumTable_main_flex, j = c("species", "cultivar", "hemisphere", "ratiomid2Early_585"), part = "header")
sumTable_main_flex <- theme_vanilla(sumTable_main_flex)
sumTable_main_flex <- fix_border_issues(sumTable_main_flex)
sumTable_main_flex <- autofit(sumTable_main_flex)
sumTable_main_flex <- width(sumTable_main_flex, j = 7, width = 1.5)
sumTable_main_flex <- width(sumTable_main_flex, j = 6, width = 1.0)
sumTable_main_flex <- width(sumTable_main_flex, j = 3, width = 0.75)
#sumTable_main_flex <- fit_to_width(sumTable_main_flex, 8, inc = 1L, max_iter = 20)

# sumTable_main_flex <- footnote(sumTable_main_flex, i = NULL, j = c(5,7),
#                           value = as_paragraph(
#                             c("SSP585")
#                           ),
#                           ref_symbols = "a", part = "header", inline = FALSE, sep = "; "
# )
sumTable_main_flex <- align(sumTable_main_flex, align = "center", i = 1, j = 4, part = "header")
sumTable_main_flex <- height(sumTable_main_flex, height = .25, part = "body")
#sumTable_main_flex <- add_footer_row(sumTable_main_flex, values = "This is a note in footer" )
#sumTable_main_flex <- merge_at(sumTable_main_flex, j = 1:5, part = "footer")

#area only flextable -----
col_keys_area = names(sumTable_main)[!grepl("ratio", names(sumTable_main), fixed = TRUE)]
sumTable_main_flex_area <- flextable(sumTable_main, col_keys = col_keys_area)
typology_area <- data.frame(
  col_keys = col_keys_area,
  # what = c("Area (sq. km)", "Area (sq. km)","Area (sq. km)", #3
  #          "Area (sq. km)", "Area (sq. km)","Area (sq. km)","Area (sq. km)","Area (sq. km)","Area (sq. km)","Area (sq. km)","Area (sq. km)","Area (sq. km)"), #9
  measure = c("Species", "Cultivar", "Hemi-\nsphere", #3
              "Early century", "mid century, SSP1-2.6", "end century, SSP1-2.6", "mid century, SSP5-8.5", "end century, SSP5-8.5", 
              "mid century, SSP1-2.6, common", "end century, SSP1-2.6, common", "mid century, SSP5-8.5, common", "end century, SSP5-8.5, common"), # 12
  stringsAsFactors = FALSE )

sumTable_main_flex_area <- set_header_df(sumTable_main_flex_area, mapping = typology_area)
sumTable_main_flex_area <- merge_h(sumTable_main_flex_area, part = "header")
# sumTable_main_flex_area <- merge_v(sumTable_main_flex_area, j = c("species", "cultivar", "hemisphere", "ratiomid2Early_585"), part = "header")
sumTable_main_flex_area <- theme_vanilla(sumTable_main_flex_area)
sumTable_main_flex_area <- fix_border_issues(sumTable_main_flex_area)
sumTable_main_flex_area <- autofit(sumTable_main_flex_area)
sumTable_main_flex_area <- width(sumTable_main_flex_area, j = 7, width = 1.0)
sumTable_main_flex_area <- width(sumTable_main_flex_area, j = 6, width = 1.0)
sumTable_main_flex_area <- width(sumTable_main_flex_area, j = 3, width = 0.75)
sumTable_main_flex_area <- align(sumTable_main_flex_area, align = "center", i = 1, j = 4, part = "header")
sumTable_main_flex_area <- height(sumTable_main_flex_area, height = .25, part = "body")
sumTable_main_flex_area <- add_footer(sumTable_main_flex_area, values = "This is a note in footer" )

prsect_area <- prop_section(
  page_size = page_size(width = 8, height = 11, orient = "landscape"),
  page_margins = page_mar(
    bottom = 1, top = 1,
    right = 1, left = 1,
    header = 0.5, footer = 0.5, gutter = 0.5
  ),
  type = NULL,
  section_columns = NULL
)
save_as_docx(sumTable_main_flex_area, values = NULL, path = "results/flextable_area_main.docx", pr_section = prsect_area)

flextable_dim(sumTable_main_flex_area)

dim_pretty(sumTable_main_flex_area)

#ratios only flextable -----
col_keys_ratios = c("species", "cultivar", "hemisphere", names(sumTable_main)[grepl("ratio", names(sumTable_main), fixed = TRUE)])
sumTable_main_flex_ratios <- flextable(sumTable_main, col_keys = col_keys_ratios)
typology_ratios <- data.frame(
  col_keys = col_keys_ratios,
  measure = c("Species", "Cultivar", "Hemi-\nsphere", #3
              "Mid to early century, SSP1-2.6", "Mid to early century, SSP5-8.5", "Common to early century, SSP1-2.6", "Common to early century, SSP5-8.5", 
              "Loss of early to SSP1-2.6", "Loss of early to SSP5-8.5"), #6
  stringsAsFactors = FALSE )

sumTable_main_flex_ratios <- set_header_df(sumTable_main_flex_ratios, mapping = typology_ratios)
sumTable_main_flex_ratios <- merge_h(sumTable_main_flex_ratios, part = "header")
# sumTable_main_flex_ratios <- merge_v(sumTable_main_flex_ratios, j = c("species", "cultivar", "hemisphere", "ratiomid2Early_585"), part = "header")
sumTable_main_flex_ratios <- theme_vanilla(sumTable_main_flex_ratios)
sumTable_main_flex_ratios <- fix_border_issues(sumTable_main_flex_ratios)
sumTable_main_flex_ratios <- autofit(sumTable_main_flex_ratios)
sumTable_main_flex_ratios <- width(sumTable_main_flex_ratios, j = 7, width = 1.0)
sumTable_main_flex_ratios <- width(sumTable_main_flex_ratios, j = 6, width = 1.0)
sumTable_main_flex_ratios <- width(sumTable_main_flex_ratios, j = 3, width = 0.75)
sumTable_main_flex_ratios <- align(sumTable_main_flex_ratios, align = "center", i = 1, j = 4, part = "header")
sumTable_main_flex_ratios <- height(sumTable_main_flex_ratios, height = .25, part = "body")
prsect_ratios <- prop_section(
  page_size = page_size(width = 8, height = 11, orient = "landscape"),
  page_margins = page_mar(
    bottom = 1, top = 1,
    right = 1, left = 1,
    header = 0.5, footer = 0.5, gutter = 0.5
  ),
  type = NULL,
  section_columns = NULL
)
save_as_docx(sumTable_main_flex_ratios, values = NULL, path = "results/flextable_ratios_main.docx", pr_section = prsect_ratios)

flextable_dim(sumTable_main_flex_ratios)

dim_pretty(sumTable_main_flex_ratios)


# results for lo chill portions cultivars -----
sumTable_lo <- as.data.table(read.csv(file = paste0(locOfDataFiles_perennials, "sumTable", "_lo", ".csv")))
sumTable_lo[, species := gsub("_lo", "", species)]
sumTable_lo[, species := stringr::str_to_title(species)][species == "Winegrape", species := "Wine grape"]

namesOrderOrig <- names(sumTable_lo)
namesOrderNew <- c("species", "cultivar", "hemisphere", "historical_1991_2010", "ssp126_2041_2060", "ssp126_2081_2100", "ssp585_2041_2060", "ssp585_2081_2100", 
                   "ssp126_2041_2060_common", "ssp585_2041_2060_common", 
                   "ssp126_2081_2100_common", "ssp585_2081_2100_common", 
                   "ratiomid2Early_126", "ratiomid2Early_585", 
                   "ratioCommon2Early_126", "ratioCommon2Early_585", 
                   "ratioLossOfEarly_126", "ratioLossOfEarly_585")
setcolorder(sumTable_lo, namesOrderNew)
sumTable_lo_flex <- flextable(sumTable_lo)
# sumTable_lo_flex <- set_header_labels(sumTable_lo_flex, values = list(
#   species = "species", cultivar = "cultivar", hemisphere = "hemisphere", 
#   historical = "early century area", ssp585 = "mid century area, SSP585", ssp585_common = "common area, SSP mid century and early century", ratiomid2Early_585 = "area ratio, mid spp585 to beginning"))
typology_all <- data.frame(
  col_keys_all = names(sumTable_lo),
  what = c("Species", "Cultivar", "Hemi-\nsphere", #3
           "Area (sq. km)", "Area (sq. km)","Area (sq. km)","Area (sq. km)","Area (sq. km)","Area (sq. km)","Area (sq. km)","Area (sq. km)","Area (sq. km)", #9
           "Area ratio", "Area ratio", "Area ratio", "Area ratio", "Area ratio", "Area ratio"), #15
  measure = c("Species", "Cultivar", "Hemi-\nsphere", #3
              "Early century", "mid century, SSP1-2.6", "end century, SSP1-2.6", "mid century, SSP5-8.5", "end century, SSP5-8.5", 
              "mid century, SSP1-2.6, common", "end century, SSP1-2.6, common", "mid century, SSP5-8.5, common", "end century, SSP5-8.5, common", #3
              "Mid to early century, SSP1-2.6", "Mid to early century, SSP5-8.5", "Common to early century, SSP1-2.6", "Common to early century, SSP5-8.5", "Loss of early to 126", "Loss of early to 585"), # 12
  stringsAsFactors = FALSE )

sumTable_lo_flex <- set_header_df(sumTable_lo_flex, mapping = typology_all, key = "col_keys_all")
sumTable_lo_flex <- merge_h(sumTable_lo_flex, part = "header")
sumTable_lo_flex <- merge_v(sumTable_lo_flex, j = c("species", "cultivar", "hemisphere", "ratiomid2Early_585"), part = "header")
sumTable_lo_flex <- theme_vanilla(sumTable_lo_flex)
sumTable_lo_flex <- fix_border_issues(sumTable_lo_flex)
sumTable_lo_flex <- autofit(sumTable_lo_flex)
sumTable_lo_flex <- fit_to_width(sumTable_lo_flex, max_width = 10, inc = 1L, max_iter = 20, unit = "in")

sumTable_lo_flex <- width(sumTable_lo_flex, j = 7, width = 1.5)
sumTable_lo_flex <- width(sumTable_lo_flex, j = 6, width = 1.0)
sumTable_lo_flex <- width(sumTable_lo_flex, j = 3, width = 0.75)
#sumTable_lo_flex <- fit_to_width(sumTable_lo_flex, 8, inc = 1L, max_iter = 20)

# sumTable_lo_flex <- footnote(sumTable_lo_flex, i = NULL, j = c(5,7),
#                           value = as_paragraph(
#                             c("SSP585")
#                           ),
#                           ref_symbols = "a", part = "header", inline = FALSE, sep = "; "
# )
sumTable_lo_flex <- align(sumTable_lo_flex, align = "center", i = 1, j = 4, part = "header")
sumTable_lo_flex <- height(sumTable_lo_flex, height = .25, part = "body")

#area only flextable -----
col_keys_area = names(sumTable_lo)[!grepl("ratio", names(sumTable_lo), fixed = TRUE)]
sumTable_lo_flex_area <- flextable(sumTable_lo, col_keys = col_keys_area)
typology_area <- data.frame(
  col_keys = col_keys_area,
  # what = c("Area (sq. km)", "Area (sq. km)","Area (sq. km)", #3
  #          "Area (sq. km)", "Area (sq. km)","Area (sq. km)","Area (sq. km)","Area (sq. km)","Area (sq. km)","Area (sq. km)","Area (sq. km)","Area (sq. km)"), #9
  measure = c("Species", "Cultivar", "Hemi-\nsphere", #3
              "Early century", "mid century, SSP1-2.6", "end century, SSP1-2.6", "mid century, SSP5-8.5", "end century, SSP5-8.5", 
              "mid century, SSP1-2.6, common", "end century, SSP1-2.6, common", "mid century, SSP5-8.5, common", "end century, SSP5-8.5, common"), # 12
  stringsAsFactors = FALSE )

sumTable_lo_flex_area <- set_header_df(sumTable_lo_flex_area, mapping = typology_area)
sumTable_lo_flex_area <- merge_h(sumTable_lo_flex_area, part = "header")
# sumTable_lo_flex_area <- merge_v(sumTable_lo_flex_area, j = c("species", "cultivar", "hemisphere", "ratiomid2Early_585"), part = "header")
sumTable_lo_flex_area <- theme_vanilla(sumTable_lo_flex_area)
sumTable_lo_flex_area <- fix_border_issues(sumTable_lo_flex_area)
sumTable_lo_flex_area <- autofit(sumTable_lo_flex_area)
sumTable_lo_flex_area <- width(sumTable_lo_flex_area, j = 7, width = 1.0)
sumTable_lo_flex_area <- width(sumTable_lo_flex_area, j = 6, width = 1.0)
sumTable_lo_flex_area <- width(sumTable_lo_flex_area, j = 3, width = 0.75)
sumTable_lo_flex_area <- align(sumTable_lo_flex_area, align = "center", i = 1, j = 4, part = "header")
sumTable_lo_flex_area <- height(sumTable_lo_flex_area, height = .25, part = "body")
prsect_area <- prop_section(
  page_size = page_size(width = 8, height = 11, orient = "landscape"),
  page_margins = page_mar(
    bottom = 1, top = 1,
    right = 1, left = 1,
    header = 0.5, footer = 0.5, gutter = 0.5
  ),
  type = NULL,
  section_columns = NULL
)
save_as_docx(sumTable_lo_flex_area, values = NULL, path = "results/flextable_area_lo.docx", pr_section = prsect_area)

flextable_dim(sumTable_lo_flex_area)

dim_pretty(sumTable_lo_flex_area)

#ratios only flextable -----
col_keys_ratios = c("species", "cultivar", "hemisphere", names(sumTable_lo)[grepl("ratio", names(sumTable_lo), fixed = TRUE)])
sumTable_lo_flex_ratios <- flextable(sumTable_lo, col_keys = col_keys_ratios)
typology_ratios <- data.frame(
  col_keys = col_keys_ratios,
  measure = c("Species", "Cultivar", "Hemi-\nsphere", #3
              "Mid to early century, SSP1-2.6", "Mid to early century, SSP5-8.5", "Common to early century, SSP1-2.6", "Common to early century, SSP5-8.5", 
              "Loss of early to SSP1-2.6", "Loss of early to SSP5-8.5"), #6
  stringsAsFactors = FALSE )

sumTable_lo_flex_ratios <- set_header_df(sumTable_lo_flex_ratios, mapping = typology_ratios)
sumTable_lo_flex_ratios <- merge_h(sumTable_lo_flex_ratios, part = "header")
# sumTable_lo_flex_ratios <- merge_v(sumTable_lo_flex_ratios, j = c("species", "cultivar", "hemisphere", "ratiomid2Early_585"), part = "header")
sumTable_lo_flex_ratios <- theme_vanilla(sumTable_lo_flex_ratios)
sumTable_lo_flex_ratios <- fix_border_issues(sumTable_lo_flex_ratios)
sumTable_lo_flex_ratios <- autofit(sumTable_lo_flex_ratios)
sumTable_lo_flex_ratios <- width(sumTable_lo_flex_ratios, j = 7, width = 1.0)
sumTable_lo_flex_ratios <- width(sumTable_lo_flex_ratios, j = 6, width = 1.0)
sumTable_lo_flex_ratios <- width(sumTable_lo_flex_ratios, j = 3, width = 0.75)
sumTable_lo_flex_ratios <- align(sumTable_lo_flex_ratios, align = "center", i = 1, j = 4, part = "header")
sumTable_lo_flex_ratios <- height(sumTable_lo_flex_ratios, height = .25, part = "body")
prsect_ratios <- prop_section(
  page_size = page_size(width = 8, height = 11, orient = "landscape"),
  page_margins = page_mar(
    bottom = 1, top = 1,
    right = 1, left = 1,
    header = 0.5, footer = 0.5, gutter = 0.5
  ),
  type = NULL,
  section_columns = NULL
)
save_as_docx(sumTable_lo_flex_ratios, values = NULL, path = "results/flextable_ratios_lo.docx", pr_section = prsect_ratios)


#suitable locations graphics, historical -----
k <- "historical"
l <- 1991
for (speciesChoice in speciesChoices) {
  print(paste0("Suitability level:  ", suitabilityLevel, ", speciesChoice: ", speciesChoice, ", ssp choice: ", k, ", start year: ", l))
  f_suitableLocsGraphics(k, l, speciesChoice, suitabilityLevel)
}

#suitable locations graphics, scenarios -----
for (k in sspChoices) {
  for (l in startYearChoices) {
    for (speciesChoice in speciesChoices) {
      print(paste0("Suitability level:  ", suitabilityLevel, ", speciesChoice: ", speciesChoice, ", ssp choice: ", k, ", start year: ", l))
      f_suitableLocsGraphics(k, l, speciesChoice, suitabilityLevel)
    }
  }
}

#create facet maps for good suitability -----
for (speciesChoice in speciesChoices) {
  f_suitableLocsGraphics_clean(speciesChoice)
}

# demo slides ------
k = "historical"
l = 1991
yearSpan <- paste0(l, "_", l + yearRange)
speciesChoice <- "cherry_main"
legendTitle <- "Suitable"
defaultColor <- c("green", "red")
CPfruit <- cropVals[cropName == speciesChoice, chill_portions]
gddsFruit <- cropVals[cropName == speciesChoice, gdd]
excessiveHeat <- cropVals[cropName == speciesChoice, summer_heat_threshold]

caption <- paste0("Note: ", suitabilityLevel, " growing conditions for ", strsplit(speciesChoice, "_")[[1]][1], ", ", cropVals[cropName == speciesChoice, cultivar], " variety", " include chill portions of at least ", CPfruit, ", at least ", gddsFruit, " growing degree days, \nno more than ", frostRiskDays[2], " days of frost risk during spring and no more than ", heatRiskDays[2], " days of excessive heat (above ", excessiveHeat,  "°C) in summer.")

fileName_nonlimiting_all_NH_in <- paste0(locOfDataFiles_perennials, "nonlimiting_all_", speciesChoice, "_",  k, "_", suitabilityLevel, "_", "NH", "_", yearSpan, ".tif")
suitable_combined_NH <- rast(fileName_nonlimiting_all_NH_in)
suitable_NH <- subset(suitable_combined_NH, "combinedSuit") 
extremeCold_NH <- subset(suitable_combined_NH, "extremeColdSuit") 
unsuitableFrost_NH <- subset(suitable_combined_NH, "springFrostSuit")
unsuitableHeat_NH <- subset(suitable_combined_NH, "heatSuit")
chillPortionsCutoff_NH <- subset(suitable_combined_NH, "chillPortionsSuit")
gddsSuitable_NH <- subset(suitable_combined_NH, "gddsSuit") 

fileName_nonlimiting_all_SH_in <- paste0(locOfDataFiles_perennials, "nonlimiting_all_", speciesChoice, "_",  k, "_", suitabilityLevel, "_", "SH", "_", yearSpan, ".tif")
suitable_combined_SH <- rast(fileName_nonlimiting_all_SH_in)
suitable_SH <- subset(suitable_combined_SH, "combinedSuit") 
extremeCold_SH <- subset(suitable_combined_SH, "extremeColdSuit") 
unsuitableFrost_SH <- subset(suitable_combined_SH, "springFrostSuit")
unsuitableHeat_SH <- subset(suitable_combined_SH, "heatSuit")
chillPortionsCutoff_SH <- subset(suitable_combined_SH, "chillPortionsSuit")
gddsSuitable_SH <- subset(suitable_combined_SH, "gddsSuit") 

extremeCold <- merge(extremeCold_NH,  extremeCold_SH)
unsuitableFrost <- merge(unsuitableFrost_NH, unsuitableFrost_SH)
unsuitableHeat <- merge(unsuitableHeat_NH, unsuitableHeat_SH)
chillPortionsCutoff <- merge(chillPortionsCutoff_NH, chillPortionsCutoff_SH)
gddsSuitable <- merge(gddsSuitable_NH, gddsSuitable_SH)


# demo slide info
k = "historical"
l <- 1991
yearSpan <- paste0(l, "_", l + yearRange)
r <- extremeCold
col <- rev(defaultColor)
tminExtremeVal <- cropVals[cropName == speciesChoice, low_temp_threshold]
titleText <- paste0("Temperatures below ", tminExtremeVal, "°C for at least one day can kill ", strsplit(speciesChoice, "_")[[1]][1], ". \nResults are for the ", k, " scenario, ", "period ", gsub("_", "-", yearSpan), ". \nGreen areas are not limited by this metric.")
fileName_out <- paste0(lofOfGraphicsFiles, "perennials/demoSlide_extremeColdMask_", speciesChoice, "_", suitabilityLevel, "_", k, "_", yearSpan, ".png")
f_graphics_demoSlides(r, titleText, caption, fileName_out, col)

r <- chillPortionsCutoff # 1 is where chill portions is adequate 
col <- rev(defaultColor)
titleText <- paste0("Locations where chill portions are adequate for ", strsplit(speciesChoice, "_")[[1]][1], ", ", cropVals[cropName == speciesChoice, cultivar], " variety", ". \nResults are for the ", k, " scenario, ", "period ", gsub("_", "-", yearSpan))
fileName_out <- paste0(lofOfGraphicsFiles, "perennials/demoSlide_chillPortionsCutoff_", speciesChoice, "_", suitabilityLevel, "_", k, "_", yearSpan, ".png")
f_graphics_demoSlides(r, titleText, caption, fileName_out, col)

{
  # these next lines of code must stay in this order.
  # titleText <- "Locations where the chill portions are not affected by extreme cold"
  # fileName_out <- paste0(lofOfGraphicsFiles, "perennials/demoSlide_chillPortionsGood", speciesChoice, "_", suitabilityLevel, "_", k, "_", yearSpan, ".png")
  # f_graphics_demoSlides(r, titleText, caption, fileName_out, col)
  
  r <- unsuitableFrost
  col <- rev(defaultColor)
  titleText = paste0("Locations where the frost days during spring are not limiting for ", strsplit(speciesChoice, "_")[[1]][1], ". \nResults are for the ", k, " scenario, ", "period ", gsub("_", "-", yearSpan))
  fileName_out <- paste0(lofOfGraphicsFiles, "perennials/demoSlide_springFrostRisk_", speciesChoice, "_", suitabilityLevel, "_", k, "_", yearSpan, ".png")
  f_graphics_demoSlides(r, titleText, caption, fileName_out, col)
  
  r <- unsuitableHeat # 1 is where heat Ct is above the unsuitable level
  col <- rev(defaultColor)
  titleText = paste0("Locations where summer hot days for ",  strsplit(speciesChoice, "_")[[1]][1], " are not limiting.", "\nResults are for the ", k, " scenario, ", "period ", gsub("_", "-", yearSpan))
  fileName_out <- paste0(lofOfGraphicsFiles, "perennials/demoSlide_summerHeat_", speciesChoice, "_", suitabilityLevel, "_", k, "_", yearSpan, ".png")
  f_graphics_demoSlides(r, titleText, caption, fileName_out, col)
  
  r <- gddsSuitable # 1 is where heat Ct is above the unsuitable level
  col <- rev(defaultColor)
  titleText = paste0("Locations where the growing degree days requirement for ",  strsplit(speciesChoice, "_")[[1]][1], " is met.", "\nResults are for the ", k, " scenario, ", "period ", gsub("_", "-", yearSpan))
  fileName_out <- paste0(lofOfGraphicsFiles, "perennials/demoSlide_gdds_", speciesChoice, "_", suitabilityLevel, "_", k, "_", yearSpan, ".png")
  f_graphics_demoSlides(r, titleText, caption, fileName_out, col)
}

# r <- heatCt
# col <- defaultColor
# titleText = "Locations in red are where the number of summer hot days are not suitable."
# fileName_out <- paste0(lofOfGraphicsFiles, "perennials/demoSlide_summerHeat_", speciesChoice, "_", suitabilityLevel, "_", k, "_", yearSpan, ".png")
# f_graphics_demoSlides(r, titleText, caption, fileName_out, col)

# powerpoint, suitable regions -----
{library(officer)
  library(magrittr)
  
  # for two stacked graphs such and NH and SH
  # defaultWidth <- 10
  # defaultHeight <- 4
  # defaultLeft <- 0
  # defaultTopNH <- 0.5
  # defaultTopSH <- 4
  
  # for one global graph
  defaultWidth <- 9
  defaultHeight <- 7
  defaultLeft <- .5
  defaultTop <- 1
  
  
  # f_perennialStressPpt <- function(season) {
  #   if (season == "spring") fileNameStart <- paste0("springFrost_hem_")
  #   if (season == "summer") fileNameStart <- paste("summerHeat_hem_")
  #   fileName_NH <- paste0(lofOfGraphicsFiles, "perennials/", fileNameStart, "NH", "_", k, "_", yearSpan, ".png")
  #   fileName_SH <- paste0(lofOfGraphicsFiles, "perennials/", fileNameStart, "SH", "_", k, "_", yearSpan, ".png")
  #   
  #   extImg_NH <- external_img(src = fileName_NH, width = defaultWidth, height = defaultHeight)
  #   extImg_SH <- external_img(src = fileName_SH, width = defaultWidth, height = defaultHeight)
  #   
  #   my_pres <- add_slide(x = my_pres, layout = 'Title Only', master = 'Office Theme')
  #   my_pres <- ph_with(x = my_pres, value = extImg_NH, location = ph_location(left = defaultLeft, top = defaultTopNH, width = defaultWidth, height = defaultHeight - 0.5) )
  #   my_pres <- ph_with(x = my_pres, value = extImg_SH, location = ph_location(left = defaultLeft, top = defaultTopSH, width = defaultWidth, height = defaultHeight - 0.5) )
  #   return(my_pres)
  # }
  
  # presentation text and figures 
  titleString <- paste0("Climate change and temperature limits for perennial crops")
  contentString <- paste0("Powerpoint produced on ", Sys.Date())
  
  dataText1 <- "The climate data set used in these graphics was prepared initially by the ISIMIP project (www.isimip.org) using CMIP6 data. " 
  dataText2 <- "This analysis uses the ISIMIP3b output data sets (https://www.isimip.org/news/isimip3ab-protocol-released/, as corrected in Spring 2021). "
  dataText3 <- "It includes data from 5 earth system models (GFDL-ESM4, UKESM1-0-LL, MPI-ESM1-2-HR, MRI-ESM2-0, and IPSL-CM6A-LR) and three scenarios (ssp126, ssp370 and ssp585). In this powerpoint, only results using ssp 126 and ssp585 are presented. " 
  dataText4 <- "The data from a 20 year period for the individual models are averaged. "
  #dataText5 <- "Crop area is based on the SAGE cropping calendar data set, based in the early 2000s. Areas that are not cropped have an NA value and are displayed in white. Areas in gray have chilling less than the minimum requirements. Areas in yellow have chilling hours between the lower and upper range of the requirements."
  dataText <- c(dataText1, dataText2, dataText3, dataText4) #, dataText5)
  
  introText1 <- "Temperate perennial fruits face five temperature-based challenges - temperature extremes of cold and heat can severely damage or kill the plant, a minimum period of time during the winter season where temperatures are below freezing (a chilling period) is required, a minimum number of growing degree days is needed and periods during the growing season where high temperatures can cause damage in critical growth periods. These challenges, and specific values, are described in more detail in the table and text on the next page."
  introText2 <- "The initial slides in this ppt show how these weather-based challenges combine to identify locations across the globe today where these temperatures metrics are not limiting.  The remaining slides are specific to individual perennial fruits and show how non-limited locations change in the future under different climate scenarios ."
  introText3 <- 'The next slide shows the values for each of the metrics and their time frame.' 
  introText4 <- 'The definitions of the metrics are: "low temp threshold" - temperature that can kill a plant after a day\'s exposure; "frost threshold" - minimum temperature at which frost damage occurs; "summer heat threshold" - temperature at which heat damages occurs; "GDD base" and "GDD upper optimum" - range of temperatures where growing degree days accumulate; "Required GDDs" - minimum number of growing degree days to reach maturity.'
  introText5 <- paste0("The frost (daily minimum temperature of less than -2°C) damage window extmids for ", springFrostLength, " days, beginning on day ", springStart_NH, " in the northern hemisphere and day ", springStart_SH, " in the southern hemisphere. Frost is not limiting if the number of frost days is less than ", frostRiskDays[2], 
                       ". The heat damage window extmids for ", heatDamageLength, " days, beginning on day ", heatDamageStart_NH, " in the northern hemisphere and day ", heatDamageStart_SH, " in the southern hemisphere. ", 
                       "The chill portion window extmids for ", chillPortionWindow, " days, beginning on day ", chillPortionStart_NH, " in the northern hemisphere and day ", chillPortionStart_SH, " in the southern hemisphere. High heat is not limiting if the number of hot days is less than ", heatRiskDays[2], 
                       ". The growing degree window begins in the first day of a window of at least 100 days with temperatures above the frost threshold.")
  introText6 <- "The table below illustrates the hemisphere-specific impacts of climate change on temperature-based suitability metrics on four perennial species. Early century non-limited area is for the period 1991-2010. mid century suitable area is for the period 2081-2100 for the SSP585 scenario. Climate change results both in shifts in non-limited areas towards the poles, a net increase in non-limited areas in the northern hemisphere and a net loss in non-limited areas in the southern hemisphere."
  cropValues <- copy(cropVals)
  setcolorder(cropValues, neworder =  c("cropName","cultivar", "low_temp_threshold", "frost_threshold", "summer_heat_threshold", "gdd", "gddtb", "GDD_opt", "chill_portions"))
  cropValues[, cropName := gsub("_main", "", cropName)]
  setnames(cropValues, old = names(cropValues), new = gsub("_", " ", names(cropValues)))
  setnames(cropValues, old = c("cropName", "gdd", "gddtb", "GDD opt", "chill portions"), new = c("crop name", "Required GDDs", "GDD base", "GDD upper optimum", "Required chill portions"))
  
  fp_1 <- fp_text(bold = TRUE, color = "pink", font.size = 0)
  fp_2 <- fp_text(bold = FALSE, font.size = 12)
  fp_3 <- fp_text(italic = TRUE, color = "black", font.size = 14)
  
  blIntro <- block_list(
    fpar(
      ftext(introText1, fp_2)),
    fpar(),
    fpar(
      ftext(introText2, fp_2))
  )
  
  blIntro2 <- block_list(
    fpar(
      ftext(introText4, fp_2)),
    fpar(),
    fpar(
      ftext(introText5, fp_2))
  )
  
  blIntro3 <- block_list(
    fpar(
      ftext(introText6, fp_2))
  )
  
  blData <- block_list(
    #  fpar(ftext("hello world", fp_1)),
    fpar(
      ftext(dataText1, fp_2),
      ftext(dataText2, fp_2),
      ftext(dataText3, fp_2),
      ftext(dataText4, fp_2) #,
      #    ftext(dataText5, fp_2)
    ))
}

my_pres <- read_pptx()
my_pres <- add_slide(x = my_pres, layout = 'Title Slide', master = 'Office Theme')
my_pres <- ph_with(x = my_pres, value = titleString, location = ph_location_type(type = "ctrTitle"))
my_pres <- ph_with(x = my_pres, value = contentString, location = ph_location_type(type = "subTitle"))

# fileName_RTable <- paste0(lofOfGraphicsFiles, "perennials/", "RebeccaTable.png")
# extImg_Rtable <- external_img(src = fileName_RTable, width = defaultWidth, height = defaultHeight)
my_pres <- add_slide(my_pres, layout = "Title and Content", master = "Office Theme")
my_pres <-  ph_with(x = my_pres, value = "Introduction", location = ph_location_type(type = "title"))
#my_pres <- ph_with(x = my_pres, value = extImg_Rtable, location = ph_location(left = defaultLeft, top = defaultHeight/2, width = 8, height = 3) )
my_pres <- ph_with(x = my_pres, value = blIntro, location = ph_location_type(type = "body") )

my_pres <- add_slide(my_pres, layout = "Title and Content", master = "Office Theme")
my_pres <-  ph_with(x = my_pres, value = "Metrics", location = ph_location_type(type = "title"))
my_pres <- ph_with(x = my_pres, value = blIntro2, location = ph_location_type(type = "body") )
my_pres <- ph_with(x = my_pres, value = flextable(cropVals, cwidth = 1.0), location = ph_location(left = defaultLeft, top = defaultHeight/1.45, width = 8, height = 3) )
my_pres <- add_slide(my_pres, layout = "Title and Content", master = "Office Theme")
my_pres <-  ph_with(x = my_pres, value = "Area impacts of climate change", location = ph_location_type(type = "title"))
my_pres <- ph_with(x = my_pres, value = blIntro3, location = ph_location_type(type = "body") )
my_pres <- ph_with(x = my_pres, value = sumTable_main_flex, location = ph_location(left = defaultLeft, top = defaultHeight/2.5, width = 6, height = 4) )


# add graphics for the demo slides
k = "historical"
l <- 1991
yearSpan <- paste0(l, "_", l + yearRange)
extremeColdMask <- external_img(src =  paste0(lofOfGraphicsFiles, "perennials/demoSlide_extremeColdMask_cherry_main_good_", k, "_", yearSpan, ".png"), width = defaultWidth, height = defaultHeight)
chillPortionsCutoff <- external_img(src =  paste0(lofOfGraphicsFiles, "perennials/demoSlide_chillPortionsCutoff_cherry_main_good_", k, "_", yearSpan, ".png"), width = defaultWidth, height = defaultHeight)
chillPortions <- external_img(src =  paste0(lofOfGraphicsFiles, "perennials/demoSlide_chillPortionsCutoff_cherry_main_good_", k, "_", yearSpan, ".png"), width = defaultWidth, height = defaultHeight)
springFrost <- external_img(src =  paste0(lofOfGraphicsFiles, "perennials/demoSlide_springFrostRisk_cherry_main_good_", k, "_", yearSpan, ".png"), width = defaultWidth, height = defaultHeight)
#summerHeatBad <- external_img(src =  paste0(lofOfGraphicsFiles, "perennials/demoSlide_summerHeat_cherry_bad_ssp126_2041_2060.png", width = defaultWidth, height = defaultHeight)
summerHeatGood <- external_img(src =  paste0(lofOfGraphicsFiles, "perennials/demoSlide_summerHeat_cherry_main_good_", k, "_", yearSpan, ".png"), width = defaultWidth, height = defaultHeight)
gddsGood <- external_img(src =  paste0(lofOfGraphicsFiles, "perennials/demoSlide_gdds_cherry_main_good_", k, "_", yearSpan, ".png"), width = defaultWidth, height = defaultHeight)

# demo slides in ppt -----
my_pres <- add_slide(x = my_pres, layout = 'Section Header', master = 'Office Theme')
my_pres <- ph_with(x = my_pres, value = "Demonstration Slides", location = ph_location_type(type = "title"))

my_pres <- add_slide(x = my_pres, layout = 'Title Only', master = 'Office Theme')
my_pres <- ph_with(x = my_pres, value = extremeColdMask, location = ph_location(left = defaultLeft, top = defaultTop, width = defaultWidth, height = defaultHeight - 0.5) )

my_pres <- add_slide(x = my_pres, layout = 'Title Only', master = 'Office Theme')
my_pres <- ph_with(x = my_pres, value = chillPortionsCutoff, location = ph_location(left = defaultLeft, top = defaultTop, width = defaultWidth, height = defaultHeight - 0.5) )

# my_pres <- add_slide(x = my_pres, layout = 'Title Only', master = 'Office Theme')
# my_pres <- ph_with(x = my_pres, value = chillPortions, location = ph_location(left = defaultLeft, top = defaultTop, width = defaultWidth, height = defaultHeight - 0.5) )

my_pres <- add_slide(x = my_pres, layout = 'Title Only', master = 'Office Theme')
my_pres <- ph_with(x = my_pres, value = springFrost, location = ph_location(left = defaultLeft, top = defaultTop, width = defaultWidth, height = defaultHeight - 0.5) )

my_pres <- add_slide(x = my_pres, layout = 'Title Only', master = 'Office Theme')
my_pres <- ph_with(x = my_pres, value = summerHeatGood, location = ph_location(left = defaultLeft, top = defaultTop, width = defaultWidth, height = defaultHeight - 0.5) )

my_pres <- add_slide(x = my_pres, layout = 'Title Only', master = 'Office Theme')
my_pres <- ph_with(x = my_pres, value = gddsGood, location = ph_location(left = defaultLeft, top = defaultTop, width = defaultWidth, height = defaultHeight - 0.5) )

ensembleTitle <- paste("Locations where temperature metrics are non-limiting for perennial fruits")
#ensembleBody <- "Results for locations with good suitability." #each quality level are presented sequentially; ie 'good locations for all species, followed by 'ok' locations for all species, and then 'bad' locations for all species."
my_pres <- add_slide(x = my_pres, layout = 'Section Header', master = 'Office Theme')
my_pres <- ph_with(x = my_pres, value = ensembleTitle, location = ph_location_type(type = "title"))

# do historical first, then ssps and future periods
for (suitabilityLevel in "good") { #
  for (speciesChoice in speciesChoices) {
    ensembleTitle <- paste("Non-temperature-limited locations for ", strsplit(speciesChoice, "_")[[1]][1]) #, ", Locations quality: ", suitabilityLevel
    my_pres <- add_slide(x = my_pres, layout = 'Section Header', master = 'Office Theme')
    my_pres <- ph_with(x = my_pres, value = ensembleTitle, location = ph_location_type(type = "title"))
    
    k <- "historical"
    l <- 1991
    yearSpan <- paste0(l, "_", l + yearRange)
    my_pres <- f_suitableLocsPpt(k, l, yearSpan, suitabilityLevel, speciesChoice)
    
    for (k in sspChoices) {
      for (l in startYearChoices) {
        my_pres <- f_suitableLocsPpt(k, l, yearSpan, suitabilityLevel, speciesChoice)
      }
    }
  }
}

my_pres <- add_slide(my_pres, layout = "Title and Content", master = "Office Theme")
my_pres <-  ph_with(x = my_pres, value = "Data Sources", location = ph_location_type(type = "title"))
my_pres <- ph_with(x = my_pres, value = blData, location = ph_location_type(type = "body") )

print(my_pres, target = "presentations/cmip6/perennials/summerHeatandspringFrost.pptx") %>% browseURL

# plot points on a raster -----
k <- "historical"
l <- 1991
yearSpan <- paste0(l, "_", l + yearRange)
speciesChoice <- "winegrape_main"
suitabilityLevel = "good"
cropLocations <- as.data.table(read_excel("data-raw/crops/perennials/cropLocations.xlsx"))
#locsVector <- cropLocations[, location := NULL]
speciesName <- gsub("_main", "", speciesChoice)
if (speciesChoice == "winegrape_main") speciesName <- "wine grape"
locsVector <- cropLocations[species == speciesName,]
setnames(locsVector, old = c("lon", "lat"), new = c("x", "y"))
fileName_nonlimiting_all_in <- paste0(locOfDataFiles_perennials, "nonlimiting_all_", speciesChoice, "_",  k, "_", suitabilityLevel, "_", yearSpan, ".tif")
suit_all <- rast(fileName_nonlimiting_all_in)
suit_all <- suit_all$combinedSuit
main = paste0("Locations not temperature limited for ", speciesName, ", period: ", k)
plot(suit_all, main = main)
points(x = locsVector$x, y = locsVector$y)

suit_all_NH <- crop(suit_all, extent_NH)
plot(suit_all_NH, main = main)
points(x = locsVector$x, y = locsVector$y)

extent_Eur <- c(0, 20, 40, 55)
suit_all_eur <- crop(suit_all, extent_Eur)
plot(suit_all_eur, main = main)
points(x = locsVector$x, y = locsVector$y)

# test of graticule package
library(graticule)
xx <- seq(from = extent_Eur[1], to = extent_Eur[2], by = .5)
yy = seq(from = extent_Eur[3], to = extent_Eur[4], by = .5)
#prj <- "+proj=lcc +lon_0=150 +lat_0=-80 +lat_1=-85 +lat_2=-75 +ellps=WGS84"
prj <- "+proj=longlat +datum=WGS84 +no_defs"
plot(suit_all_eur)
grat <- graticule(lons = xx, lats = yy,  xlim = c(extent_Eur[1], extent_Eur[2]), c(ylim = extent_Eur[3], extent_Eur[4]), proj = prj)
plot(grat, add = TRUE)
labs <- graticule_labels(lons = xx, lats = yy, xline = extent_Eur[1], yline = extent_Eur[3],  proj = prj)
op <- par(xpd = NA)
#text(labs, lab = parse(text = labs$lab), pos = c(2, 1)[labs$islon + 1], adj = 1.2)
#text(subset(labs, labs$islon), lab = parse(text = labs$lab[labs$islon]), pos = 3)
#text(subset(labs, !labs$islon), lab = parse(text = labs$lab[!labs$islon]), pos = 2)
text(subset(labs, labs$islon), lab = parse(text = labs$lab[labs$islon]), pos = 1, cex=0.75, srt = 90)
text(subset(labs, !labs$islon), lab = parse(text = labs$lab[!labs$islon]), pos = 2, cex=0.75)
par(op)
points(x = locsVector$x, y = locsVector$y)

