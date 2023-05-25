# functions for the perennials calculations

require(terra)
require(sf)
require(ggplot2)
# source("R/ISIMIPconstants.R")
# source("R/ISIMIPspatialConstants.R")
library(Rcpp)
sourceCpp("R/cpp/gdd.cpp")

theme_custom <- function() {
  theme_bw() +
    theme(
      legend.text.align = 1,
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(size = 12, hjust = 0.5),
      plot.caption = element_text(hjust = 0, vjust = 7.0, size = 7)
    )
}

options(warn = 1)
# file locations -----
path_clim <- "climdata/"
path_gdds <- "data/growingDegreeDays/"

# choose whether to do the base cps, or the lo  cp requirements varieties -----
suitabilityLevels <- c("good")
yearRangeSH <- 18 # one less year because of 6 month offset

# runsParms list created here -----
climVal <- -2 # changed from 0 because subject experts say plants can tolerate this
test_logic <- paste0("x > ", climVal)
logicDirection <- ">"
if (logicDirection == ">") ldtext <- "gt"
if (logicDirection == "<") ldtext <- "lt"
runlengthChoices <- c(100) # at least 100 days of tmin > 0
climateVariable <- "tasmin"
runsParms <- c(climVal, test_logic, logicDirection, ldtext, runlengthChoices, climateVariable)

# constants, perennials -----
minimumGrwSeasonLength <- 100
# All day numbers are Julian calendar days
# The spring frost window is NH = 60:120 (1Mar:30Apr); 227:288 (SH = 15Aug:15Oct)
springFrostLength <- 60
heatDamageLength <- 60
chillPortionWindow <- 214
springStart_NH <- 60 # March 1 in 2019
springStart_SH <- 227 # Aug 15 in 2019
heatDamageStart_NH <- 182 # July 1
heatDamageStart_SH <- 1 # Jan 1
chillPortionStart_NH <- 272
chillPortionStart_SH <- 92
extremeColdCutoff <- -30
# suitability day counts - first two numbers are good window, second two numbers are acceptable window; third two numbers are bad window; fourth 2 numbers are unsuitable range. These are not used in the GCB paper version of the calculations
frostRiskDays <- c(0, 6, 7, 20, 21, 45, 46, springFrostLength)
heatRiskDays <- c(0, 12, 13, 30, 31, 45, 46, heatDamageLength)

crsRob <- "+proj=robin"
coastline <- vect("data-raw/regionInformation/ne_50m_coastline/ne_50m_coastline.shp")
coastline_cropped <- crop(coastline, extent_noAntarctica )
coastline_cropped_Rob <- project(coastline_cropped, crsRob)
coastline_cropped_Rob_sf <- sf::st_as_sf(coastline_cropped_Rob)


# functions -----
gddFilesList <- function() {
  gddSumFilesCompleted <- list.files(path_gdds, full.names = TRUE)
  gddSumFilesCompleted <- gddSumFilesCompleted[grepl("gddSum", gddSumFilesCompleted, fixed = TRUE)]
  gddSumFilesCompleted <- gddSumFilesCompleted[!grepl("aux.xml", gddSumFilesCompleted, fixed = TRUE)]
  gddSumFilesCompleted <- gsub("//", "/", gddSumFilesCompleted)
}

f_harvestArea <- function(speciesName, minArea) {
  rInArea <- rast(paste0(path_harvest_data, speciesName, "/", speciesName, "_HarvestedAreaHectares.tif"))
  rInArea <- crop(rInArea, extent_noAntarctica)
  harvestArea_earlyCent <- aggregate(rInArea, fact = 6, fun = "sum") # convert 5 arc minutes to 1/2 degrees
  maskMin <- switch(speciesName,
                    "almond" = minArea,
                    "apple" = minArea,
                    "cherry" = minArea,
                    "grape" = minArea,
                    "olive" = minArea
  )
  harvestArea_earlyCent[harvestArea_earlyCent < maskMin] <- 0 # set minimum area to be greater than minArea hectares per grid cell
  return(harvestArea_earlyCent)
}

f_h <- function(regionExtent, defaultWidth) {
  x <- regionExtent[2] - regionExtent[1]
  y <- regionExtent[4] - regionExtent[3]
  h <- defaultWidth / (x / y)
}

f_getArea <- function(r, layer) {
  r_sub <- subset(r, layer)
  r_sub[r_sub < 1 | r_sub > 1] <- NA # keep only locations with value of 1
  r_sub_area <- (expanse(r_sub, unit = "km")) / 1000 # convert to 1000 sq km
  r_sub_area <- round(r_sub_area, 0)
}

f_range <- function(x, rangeVals) {
  # set locations with values outside the range to 999 as an indicator where land is for later adjustment
  x[x < rangeVals[1]] <- 999
  x[x > rangeVals[length(rangeVals)]] <- 999 # length needed because rangeVals can sometimes have more than 2 entries.
  return(x)
}

f_convertToSHdays <- function(yearNum, calIn) {
  if (lubridate::leap_year(yearNum)) {
    dayCt <- 366
  } else {
    dayCt <- 365
  }
  NHYear <- rep(1:dayCt, 1)
  SHyrStart <- rep(182:dayCt, 1) # july 1
  SHyrmid <- rep(1:181, 1)
  SHYear <- c(SHyrStart, SHyrmid)
  daysLookup <- data.frame(cbind(NHYear, SHYear))
  uniqueVals <- unique(calIn)
  for (val in 1:length(uniqueVals)) {
    calIn[calIn == uniqueVals[val]] <- daysLookup[uniqueVals[val], "SHYear"]
  }
  return(calIn)
}

f_readRast_count <- function(modelChoice, k, l, yearSpan, hem, climateVarChoice, threshold, layersToKeep, probVal) {
  print(paste0("climateVarChoice: ", climateVarChoice))
  print(paste0("threshold: ", threshold))
  fileName_in <- paste0(path_clim, modelChoice, "_", climateVarChoice, "_", k, "_", yearSpan, ".tif")
  print(paste0("climate fileName_in: ", fileName_in))
  r <- rast(fileName_in)
  # convert day numbers to calendar days to subset
  # indices values are from 1 to the number of years in the subsetted file, usually 20.
  # layersToKeep are the layer numbers to keep in each year.
  datesToKeep <- indices <- c()
  for (yearNum in l:(l + yearRange)) {
    temp <- as.Date(layersToKeep, origin = paste0(yearNum, "-01-01"))
    indices_yearNum <- rep(yearNum - l + 1, length(layersToKeep))
    indices <- c(indices, indices_yearNum)
    temp <- paste0("X", temp)
    datesToKeep <- c(datesToKeep, temp)
  }
  r <- subset(r, datesToKeep)
  r <- crop(r, get(paste0("extent_", hem)))
  # for spring frost damage
  if (climateVarChoice == "tasmin") f_ct <- function(x) (sum(x < threshold))
  if (climateVarChoice == "tasmax") f_ct <- function(x) (sum(x > threshold))
  # print(paste0("length r: ", length(names(r)), ", length indices: ", length(indices)))
  
  print(system.time(tempCt <- tapp(r, indices, f_ct)))
  print(tempCt)
  return(tempCt)
}

f_readRast_extreme <- function(modelChoice_lower, k, l, yearSpan, hem, climateVarChoice, threshold, layersToKeep, probVal) {
  yearSpan <- paste0(l, "_", l + yearRange)
  fileName_in <- paste0(path_clim, modelChoice_lower, "_", climateVarChoice, "_", k, "_", yearSpan, ".tif")
  print(paste0("climate fileName_in: ", fileName_in))
  r <- rast(fileName_in) |> crop(get(paste0("extent_", hem))) # this is slow because r has 7000+ layers
  print(r)
  indices <- format(as.Date(names(r), format = "X%Y-%m-%d"), format = "%Y") # %Y is year as 4 digit number
  indices <- as.numeric(indices)
  indices <- indices - l + 1
  extremeTemp <- tapp(r, indices, funDir) # indices are from 1 to the number of years in the input file. 365 1s, then 365 2s, etc. The funDir is the minimum or maximum value of the temp var. extremeTemp is the highest or lowest temp value in each year in each cell. #also slow
}

f_extremeCold <- function(k, l, speciesName, hem, modelChoices_lower, cropVals) {
  gddFilesCompleted <- list.files(path_data, full.names = TRUE)
  gddFilesCompleted <- gddFilesCompleted[!grepl("aux.xml", gddFilesCompleted, fixed = TRUE)]
  gddFilesCompleted <- gsub("//", "/", gddFilesCompleted)
  yearSpan <- paste0(l, "_", l + yearRange)
  climateVarChoice <- "tasmin"
  funDir <- "min"
  probVal <- 0.80
  speciesChoice <- paste0(speciesName, "_main")
  
  fileName_out <- paste0(path_data, "extremeCold_cutoff_", speciesName, "_", k, "_", hem, "_", yearSpan, ".tif")
  if (fileName_out %in% gddFilesCompleted) {
    print(paste0("Already done: ", fileName_out))
  } else {
    system.time(x <- lapply(modelChoices_lower, f_readRast_extreme, k, l, yearSpan, hem, climateVarChoice, threshold, layersToKeep, probVal)) # read in tasmin for the relevant period and all ESMs
    r <- rast(x)
    print(r)
    extremeColdCutoff <- cropVals[cropName == speciesChoice, low_temp_threshold]
    
    # now do ensemble mean and cutoff
    print(paste0("Working on extreme cold for speciesName: ", speciesName, ", working on ssp: ", k, ", start year ", l, ", hemisphere ", hem))
    print(system.time(extremeCold_quant <- quantile(r, probs = probVal, na.rm = TRUE)))
    extremeCold_quant[extremeCold_quant < extremeColdCutoff] <- 0 #  extreme cold limited
    extremeCold_quant[extremeCold_quant >= extremeColdCutoff] <- 1 # not extreme cold limited
    
    print(system.time(writeRaster(extremeCold_quant, filename = fileName_out, overwrite = TRUE, wopt = woptList)))
    print(paste0("outF extreme cold: ", fileName_out))
    plot(extremeCold_quant, main = speciesName, col = "green")
    print(paste0("fileName_out: ", fileName_out))
    return(extremeCold_quant)
  }
}

f_gdd <- function(temp, topt_min, topt_max) {
  return(clamp(temp, topt_min, topt_max, values = TRUE) - topt_min)
}

f_computeGDDs <- function(k, l, speciesName, modelChoice, topt_min, topt_max) {
  modelChoice_lower <- tolower(modelChoice)
  gddFilesCompleted <- list.files(path_gdds, full.names = TRUE)
  gddFilesCompleted <- gddFilesCompleted[!grepl("aux.xml", gddFilesCompleted, fixed = TRUE)]
  gddFilesCompleted <- gsub("//", "/", gddFilesCompleted)
  yearSpan <- paste0(l, "_", l + yearRange)
  outF <- paste0(path_gdds, modelChoice_lower, "_", "gdd", "_", speciesName, "_", k, "_", yearSpan, ".tif")
  if (outF %in% gddFilesCompleted) {
    print(paste0("Already done: ", outF))
  } else {
    yearSpan <- paste0(l, "_", l + yearRange)
    modelChoice_lower <- tolower(modelChoice)
    print(paste0("start year: ", l, ", ssp: ", k, " model: ", modelChoice, ", start year: ", l, ", speciesName: ", speciesName))
    fileName_tas_in <- paste0(path_clim, modelChoice_lower, "_", "tas_cropped", "_", k, "_", yearSpan, ".tif")
    tas <- rast(fileName_tas_in)
    print(paste0("Working on: ", outF))
    print(paste0("crop: ", speciesName, " topt_min: ", topt_min, " topt_max: ", topt_max, " outF: ", outF))
    # f_gdd is a cpp function loaded near the beginning of this code
    print(system.time(gdd <- f_gdd(tas, topt_min, topt_max)))
    print(system.time(writeRaster(gdd, filename = outF, names = names(tas), overwrite = TRUE)))
    return(gdd)
  }
}

f_gddSums <- function(k, l, speciesName, hem, runsParms, modelChoice) {
  gddSumFilesCompleted <- gddFilesList()
  yearSpan <- paste0(l, "_", l + yearRange)
  modelChoice_lower <- tolower(modelChoice)
  # growing season directions for the perennials paper
  # -- Frost free season. Last frost day calculated from day-of-year 1 to first subsequent frost
  # these values needed to generate the file names from runs calculations done in f_runs
  # logicDirection <- runsParms[3] # > for perennials paper; not needed here
  ldtext <- runsParms[4] # combines with logic direction; gt for the perennials
  climVal <- runsParms[1] # -2 # threshold for counting the day
  climateVariable <- runsParms[6] # tasmin
  runlength <- runsParms[5] # runlength must be at least this long #100
  # daily gdd values
  fileName_gdd_in <- paste0(path_gdds, modelChoice_lower, "_", "gdd", "_", speciesName, "_", k, "_", yearSpan, ".tif")
  gdds <- rast(fileName_gdd_in)
  
  # get hemisphere-specific
  gdds_hem <- crop(gdds, get(paste0("extent_", hem)))
  # gdds are daily for the 20 year period
  if (hem == "SH") {
    startDate <- paste0(l, "-07-01")
    endDate <- paste0(l + yearRange - 1, "-06-30") # in southern hemisphere search July 1 to June 30 of the next year.
    fileName_gddSums_out <- paste0(path_gdds, "gddSum", "_", modelChoice_lower, "_", hem, "_", speciesName, "_", k, "_", yearSpan, ".tif")
  }
  if (hem == "NH") {
    startDate <- paste0(l, "-01-01")
    endDate <- paste0(l + yearRange, "-12-31")
    fileName_gddSums_out <- paste0(path_gdds, "gddSum", "_", modelChoice_lower, "_", hem, "_", speciesName, "_", k, "_", yearSpan, ".tif")
  }
  if (fileName_gddSums_out %in% gddSumFilesCompleted) {
    print(paste0("Already done with ", fileName_gddSums_out))
  } else {
    indices <- seq(as.Date(startDate), as.Date(endDate), by = "days")
    indicesChar <- paste0("X", indices)
    indicesYr <- unique(as.numeric(format(indices, "%Y")))
    indicesYr <- indicesYr[1:yearRange]
    fileName_startDay1_in <- paste0(path_data, "startday_1_", climateVariable, "_", modelChoice_lower, "_run_", runlength, "_lim_", ldtext, climVal, "_", hem, "_", k, "_", yearSpan, ".tif")
    fileName_endDay1_in <- paste0(path_data, "endday_1_", climateVariable, "_", modelChoice_lower, "_run_", runlength, "_lim_", ldtext, climVal, "_", hem, "_", k, "_", yearSpan, ".tif")
    startDay <- rast(fileName_startDay1_in)
    endDay <- rast(fileName_endDay1_in)
    
    # now do calc by year
    for (yearNumber in 1:nlyr(startDay)) {
      print(paste0("Working on species choice: ", speciesName, ", sssp: ", k, ", startYear: ", l, ", hem: ", hem, ", model: ", modelChoice, ", yearNumber: ", yearNumber))
      startDay_yr <- subset(startDay, yearNumber)
      endDay_yr <- subset(endDay, yearNumber)
      startYear <- l + yearNumber - 1
      
      if (hem == "SH") {
        startDate <- paste0(startYear, "-07-01")
        endDate <- paste0(startYear + 1, "-06-30")
      } # in southern hemisphere search July 1 to June 30 of the next year. NH is just the calendar year
      if (hem == "NH") {
        startDate <- paste0(startYear, "-01-01")
        endDate <- paste0(startYear, "-12-31")
      }
      fileName_gddSums_out <- paste0(path_gdds, "gddSum", "_", modelChoice_lower, "_", hem, "_", speciesName, "_", k, "_", yearSpan, ".tif")
      
      indices <- seq(as.Date(startDate), as.Date(endDate), by = "days")
      indicesChar <- paste0("X", indices)
      if ((hem == "SH" & yearNumber < 20) | (hem == "NH")) {
        gdds_yr <- subset(gdds_hem, indicesChar)
      }
      # rapp needs to have a start day that is greater than 0. Next two lines sets all 0 values to 1
      startDay_yr[startDay_yr == 0] <- 1
      endDay_yr[endDay_yr == 0] <- 1
      endDay_yr[endDay_yr > nlyr(gdds_yr)] <- nlyr(gdds_yr) # needed because the end day can be # 367 if we're in a tropical region; end day is the frost day that ends the run
      print(system.time(sum_gdds <- rapp(gdds_yr, startDay_yr, endDay_yr, "sum")))
      plot(sum_gdds, main = paste0("Sum of gdds in a run of at least 100 days, scenario: ", k, ", period: ", yearSpan, ", year number: ", yearNumber, ", model: ", modelChoice_lower, ", crop: ", speciesName))
      if (yearNumber == 1) {
        period_sums <- sum_gdds
      } else {
        period_sums <- c(period_sums, sum_gdds)
      }
    }
    gc()
    period_sums
    print(system.time(writeRaster(period_sums, filename = fileName_gddSums_out, overwrite = TRUE, wopt = woptList)))
    flush.console()
    print(paste0("fileName_gddSums_out: ", fileName_gddSums_out))
  }
}

f_gddSum_mean <- function(k, l, speciesNames, modelChoice, hemispheres) {
  yearSpan <- paste0(l, "_", l + yearRange)
  for (speciesName in speciesNames) {
    modelChoice_lower <- tolower(modelChoice)
    for (hem in hemispheres) {
      fileName_gddSums_in <- paste0(path_gdds, "gddSum", "_", modelChoice_lower, "_", hem, "_", speciesName, "_", k, "_", yearSpan, ".tif")
      r_in <- rast(fileName_gddSums_in)
      outF <- paste0(path_gdds, "gddSum_mean", "_", modelChoice_lower, "_", hem, "_", speciesName, "_", k, "_", yearSpan, ".tif")
      test <- app(r_in, mean, filename = outF, overwrite = TRUE, wopt = woptList)
      print(paste0("outF: ", outF))
    }
  }
}

f_readRast_gddSum <- function(modelChoice, speciesName, k, l, hem) {
  yearSpan <- paste0(l, "_", l + yearRange)
  modelChoice_lower <- tolower(modelChoice)
  fileName_in <- paste0(path_gdds, "gddSum_mean", "_", modelChoice_lower, "_", hem, "_", speciesName, "_", k, "_", yearSpan, ".tif")
  print(paste0("speciesName: ", speciesName, ", k: ", k, ", modelChoice: ", modelChoice, ", fileName in: ", fileName_in))
  r <- rast(fileName_in)
}

f_ensemble_GDD_sum_mean <- function(k, l, yearRange, speciesNames, cropVals) {
  yearSpan <- paste0(l, "_", l + yearRange)
  for (speciesName in speciesNames) {
    speciesChoice <- paste0(speciesName, "_main")
    gddsRequired <- cropVals[cropName == speciesChoice, gdd]
    for (hem in hemispheres) {
      x <- lapply(modelChoices, f_readRast_gddSum, speciesName, k, l, hem)
      r <- rast(x)
      indices_day <- rep(seq(1, nlyr(x[[1]]), 1), 5) # 5 is number of models; if omitted should get the same result
      outF <- paste0(path_gdds, "ensemble_gddSum_mean", "_", hem, "_", speciesName, "_", k, "_", yearSpan, ".tif")
      print(paste0("speciesName: ", speciesName, ", ensemble ssp: ", k, ", start year: ", l, ", fileName out: ", outF))
      print(system.time(r_mean <- app(r, indices_day, fun = "mean", na.rm = TRUE, filename = outF, overwrite = TRUE, wopt = woptList)))
      main <- paste0("Perennial: ", speciesName, ", GDDs required: ", gddsRequired, ", hemisphere: ", hem, ", ssp: ", k, ", period: ", yearSpan)
      plot(r_mean, main = main, axes = FALSE)
    }
  }
}

f_gddsSuitability <- function(k, l, speciesName, hem, cropVals) {
  yearSpan <- paste0(l, "_", l + yearRange)
  speciesChoice <- paste0(speciesName, "_main")
  # get gdds ensemble mean
  fileName_in <- paste0(path_gdds, "ensemble_gddSum_mean", "_", hem, "_", speciesName, "_", k, "_", yearSpan, ".tif")
  gdds <- rast(fileName_in)
  # get gdd requirements
  gddsRequired <- cropVals[cropName == speciesChoice, gdd]
  gddsSuitable <- gdds
  gddsSuitable[gddsSuitable < gddsRequired] <- 0
  gddsSuitable[gddsSuitable >= gddsRequired] <- 1
  outF <- paste0(path_gdds, "gdds_not_limiting", "_", hem, "_", speciesName, "_", k, "_", yearSpan, ".tif")
  print(system.time(writeRaster(gddsSuitable, filename = outF, overwrite = TRUE, wopt = woptList)))
  print(paste0("outF suitable gdds: ", outF))
  
  plot(gddsSuitable, main = paste0("Adequate GDDs for ", speciesName, ", minimum required ", gddsRequired, ", hemisphere: ", hem, ", scenario: ", k, ", period: ", yearSpan))
  return(gddsSuitable)
}
# functions for runs -----
f_runs <- function(x, runlength, logicString) {
  # browser()
  
  # number of layers to be returned
  # element 1 - number of runs that meet the logic criterion
  # element 2 - the total length of all runs that meet the run criterion
  # element 3 - first start day number
  # element 4 - first end day number
  runResult <- c(NA, NA, NA, NA)
  if (is.nan(as.numeric(x[1]))) {
    return(runResult)
  }
  #   browser()
  seqLengthCode <- paste0("1{", runlength, ",}")
  # A regular expression  to get the first item of gregexpr. It says look for  run_length times See http://xenon.stanford.edu/~xusch/regexp/
  # seqLengthCode - what's the minimum length of a run that meets the logic entry
  # logicString - what condition is evaluated ; e.g., x > 0
  # g[[1]] - if positive, start position of run; if there are two entries, there are 2 runs that satisfy the condition and the values are the starting position of each run
  # attributes(g)$match.length - if positive, length of elements in run
  g <- gregexpr(seqLengthCode, paste(+eval(parse(text = logicString)), collapse = ""))[[1]] # The + converts TRUE and FALSE to 1 and 0
  if ((g[1] == -1)) { # no need to write to growing season if g returns -1, return 0s
    runResult <- c(0, 0, 0, 0)
    return(runResult)
  }
  startDays <- unlist(g) # start day(s) of the run(s)
  runlengths <- as.numeric(attributes(g)$match.length)
  startDay_r1 <- startDays[1]
  endDay_r1 <- startDay_r1 + runlengths[1]
  runLengths <- sum(runlengths)
  runResult <- c(length(startDays), runLengths, startDay_r1, endDay_r1)
  #  print(paste0("runResult: ", runResult))
  return(runResult)
}

f_runsSetup <- function(k, l, modelChoice, runsParms) {
  logicDirection <- runsParms[3]
  climVal <- runsParms[1]
  climateVariable <- runsParms[6]
  ldtext <- runsParms[4]
  runlength <- runsParms[5]
  yearSpan <- paste0(l, "_", l + yearRange)
  modelChoice_lower <- tolower(modelChoice)
  logicString <- paste0("x ", logicDirection, " ", climVal)
  fileName_in <- paste0(path_clim, modelChoice_lower, "_", climateVariable, "_", k, "_", yearSpan, ".tif")
  r <- rast(fileName_in)
  for (hem in hemispheres) {
    endYear <- l + yearRange # for NH
    if (hem == "SH") endYear <- l + yearRange - 1
    for (yearNumber in l:endYear) {
      print(paste0("working on ssp: ", k, ", start year: ", l, ", model choice: ", modelChoice, ", hemisphere: ", hem, ", year: ", yearNumber))
      if (hem == "SH") {
        startDate <- paste0(yearNumber, "-07-01")
        endDate <- paste0(yearNumber + 1, "-06-30")
      } # in southern hemisphere search July 1 to June 30 of the next year.
      if (hem == "NH") {
        startDate <- paste0(yearNumber, "-01-01")
        endDate <- paste0(yearNumber, "-12-31")
      }
      indices <- seq(as.Date(startDate), as.Date(endDate), by = "days")
      indicesChar <- paste0("X", indices)
      r_yr <- subset(r, indicesChar)
      r_yr <- crop(r_yr, get(paste0("extent_", hem)))
      #            browser()
      print(system.time(r_runs <- app(r_yr, f_runs, runlength, logicString)))
      if (yearNumber == l) {
        print(paste0("yearNumber: ", yearNumber))
        runs_ct <- subset(r_runs, 1)
        runs_length <- subset(r_runs, 2)
        startday_1 <- subset(r_runs, 3)
        endday_1 <- subset(r_runs, 4)
      } else {
        runs_ct <- c(runs_ct, subset(r_runs, 1))
        runs_length <- c(runs_length, subset(r_runs, 2))
        startday_1 <- c(startday_1, subset(r_runs, 3))
        endday_1 <- c(endday_1, subset(r_runs, 4))
      }
    }
    names(runs_ct) <- l:endYear
    names(runs_length) <- l:endYear
    names(startday_1) <- l:endYear
    names(endday_1) <- l:endYear
    
    fileName_ct_out <- paste0(path_data, "runs_ct_", climateVariable, "_", modelChoice_lower, "_run_", runlength, "_lim_", ldtext, climVal, "_", hem, "_", k, "_", yearSpan, ".tif")
    fileName_length_out <- paste0(path_data, "runs_length_", climateVariable, "_", modelChoice_lower, "_run_", runlength, "_lim_", ldtext, climVal, "_", hem, "_", k, "_", yearSpan, ".tif")
    fileName_startDay1_out <- paste0(path_data, "startday_1_", climateVariable, "_", modelChoice_lower, "_run_", runlength, "_lim_", ldtext, climVal, "_", hem, "_", k, "_", yearSpan, ".tif")
    fileName_endDay1_out <- paste0(path_data, "endday_1_", climateVariable, "_", modelChoice_lower, "_run_", runlength, "_lim_", ldtext, climVal, "_", hem, "_", k, "_", yearSpan, ".tif")
    
    writeRaster(runs_ct, filename = fileName_ct_out, overwrite = TRUE, wopt = woptList)
    writeRaster(runs_length, filename = fileName_length_out, overwrite = TRUE, wopt = woptList)
    writeRaster(startday_1, filename = fileName_startDay1_out, overwrite = TRUE, wopt = woptList)
    writeRaster(endday_1, filename = fileName_endDay1_out, overwrite = TRUE, wopt = woptList)
    print(paste0("fileName_ct_out: ", fileName_ct_out))
    print(paste0("fileName_length_out: ", fileName_length_out))
    print(paste0("fileName_startDay1_out: ", fileName_startDay1_out))
    print(paste0("fileName_endDay1_out: ", fileName_endDay1_out))
  }
  mainstart <- paste0("Start day of at least 100 consecutive days\n where temp is greater than ", climVal, "°C, \nssp: ", k, ", yearspan: ", yearSpan, ", model: ", modelChoice)
  mainend <- paste0("End day of at least 100 consecutive days\n where temp is greater than ", climVal, "°C, \nssp: ", k, ", yearspan: ", yearSpan, ", model: ", modelChoice)
  plot(startday_1[[1]], main = mainstart)
  plot(endday_1[[1]], main = mainend)
}

f_readRast_runs <- function(modelChoice, k, l, yearSpan, hem, climateVarChoice, threshold, layersToKeep, probVal) {
  logicDirection <- runsParms[3]
  climVal <- runsParms[1]
  climateVariable <- runsParms[6]
  ldtext <- runsParms[4]
  #  yearSpan <- paste0(l, "_", l + yearRange)
  modelChoice_lower <- tolower(modelChoice)
  fileName_in <- paste0(path_data, runType, "_", climateVariable, "_", modelChoice_lower, "_run_", runlength, "_lim_", ldtext, climVal, "_", hem, "_", k, "_", yearSpan, ".tif")
  print(paste0("runType: ", runType, ", k: ", k, ", modelChoice: ", modelChoice, ", fileName in: ", fileName_in))
  r <- rast(fileName_in)
}

# chillPortions calcs done in chillPortions.R -----

f_frostDamage <- function(k, l, speciesName, hem, modelChoices_lower, suitabilityLevel, cropVals) {
  yearSpan <- paste0(l, "_", l + yearRange)
  speciesChoice <- paste0(speciesName, "_main")
  # frost damage day windows ---
  spLyrStart <- get(paste0("springStart_", hem))
  spLyrend <- spLyrStart + springFrostLength
  spLyrs <- spLyrStart:spLyrend
  
  # this really doesn't need to have a range. If the second value says what the cutoff is; lower values are fine
  frDays <- switch(suitabilityLevel,
                   "good" = frostRiskDays[1:2],
                   "acceptable" = frostRiskDays[3:4],
                   "bad" = frostRiskDays[5:6],
                   "unsuitable" = frostRiskDays[7:8]
  )
  
  climateVarChoice <- "tasmin"
  threshold <- cropVals[cropName == speciesChoice, frost_threshold]
  layersToKeep <- spLyrs
  probVal <- 0.90
  system.time(x <- lapply(modelChoices_lower, f_readRast_count, k, l, yearSpan, hem, climateVarChoice, threshold, layersToKeep, probVal))
  r <- rast(x)
  r[r <= frDays[2]] <- 1
  r[r > frDays[2]] <- 0
  system.time(fr <- quantile(r, probs = probVal, na.rm = TRUE)) # note: if all layers have the same value quantile returns that value
  fr[fr > 0] <- 1
  
  titleText_fr <- paste0("Green indicates locations where frost days are not limiting for ", strsplit(speciesName, "_")[[1]][1], " during the ", k, " scenario", ", ", gsub("_", "-", yearSpan), ". \nUnsuitable frost risk is more than ", frDays[2], " frost days (-2°C) during the spring frost window.")
  plot(fr, main = titleText_fr)
  fileName_fr_out <- paste0(path_data, "frostDamage_", speciesName, "_", k, "_", suitabilityLevel, "_", hem, "_", yearSpan, ".tif")
  writeRaster(fr, fileName_fr_out, overwrite = TRUE, wopt = woptList)
  print(paste0(" frost damage fileName out: ", fileName_fr_out))
  print("--------------------------------------------")
}

f_heatDamage <- function(k, l, modelChoices_lower, speciesName, hem, suitabilityLevel, cropVals) {
  yearSpan <- paste0(l, "_", l + yearRange)
  speciesChoice <- paste0(speciesName, "_main")
  hdLyrStart <- get(paste0("heatDamageStart_", hem))
  hdLyrend <- hdLyrStart + heatDamageLength
  hdLyrs <- hdLyrStart:hdLyrend
  hdDays <- switch(suitabilityLevel,
                   "good" = heatRiskDays[1:2],
                   "acceptable" = heatRiskDays[3:4],
                   "bad" = heatRiskDays[5:6],
                   "unsuitable" = heatRiskDays[7:8]
  )
  climateVarChoice <- "tasmax"
  threshold <- cropVals[cropName == speciesChoice, summer_heat_threshold]
  layersToKeep <- hdLyrs
  probVal <- 0.90
  system.time(x <- lapply(modelChoices_lower, f_readRast_count, k, l, yearSpan, hem, climateVarChoice, threshold, layersToKeep, probVal))
  r <- rast(x)
  r[r <= hdDays[2]] <- 1
  r[r > hdDays[2]] <- 0
  system.time(hd <- quantile(r, probs = probVal, na.rm = TRUE)) # note: if all layers have the same value quantile returns that value
  hd[hd > 0] <- 1
  titleText_hd <- paste0("1 indicates locations where extreme summer heat is not limiting for ", strsplit(speciesName, "_")[[1]][1], " during the ", k, " scenario", ", ", gsub("_", "-", yearSpan), ". \nUnsuitable heat is more than ", hdDays[2], " days above ", threshold, "°C during the summer window.")
  
  plot(hd, main = titleText_hd)
  fileName_hd_out <- paste0(path_data, "heatDamage_", speciesName, "_", k, "_", suitabilityLevel, "_", hem, "_", yearSpan, ".tif")
  writeRaster(hd, fileName_hd_out, overwrite = TRUE, wopt = woptList)
  print(paste0(" heat damage fileName out: ", fileName_hd_out))
  print("--------------------------------------------")
}

f_combinedDamage <- function(k, l, speciesName, suitabilityLevel, chillLevel) {
  # combine suitability metrics from chill portions, extreme cold, spring frost, and summer heat; locations with value 1 is suitable
  yearSpan <- paste0(l, "_", l + yearRange)
  speciesChoice <- paste0(speciesName, chillLevel)
  # read in all the rasters needed
  for (hem in hemispheres) {
    fileTail <- paste0(speciesChoice, "_", k, "_", hem, "_", yearSpan, ".tif")
    fileTailSuit <- paste0(speciesChoice, "_", k, "_", suitabilityLevel, "_", hem, "_", yearSpan, ".tif")
    extremeColdCt <- rast(paste0(path_data, "extremeCold_cutoff_", fileTail <- paste0(speciesName, "_", k, "_", hem, "_", yearSpan, ".tif")))
    frostCt <- rast(paste0(path_data, "frostDamage_", fileTailSuit))
    heatCt <- rast(paste0(path_data, "heatDamage_", fileTailSuit))
    chillPortionsCutoff <- rast(paste0("data/chillPortions/chill_portions/", "ensemble_chill_cutoff_", paste0(speciesChoice, "_", k, "_", hem, "_", yearSpan, ".tif")))
    gdds_suitable <- rast(paste0(path_gdds, "gdds_not_limiting", "_", hem, "_", speciesName, "_", k, "_", yearSpan, ".tif"))
    
    print(paste0("working on combined damage ", speciesName, " in hemisphere ", hem, ", year ", l, ", scenario ", k))
    
    r_combined <- c(extremeColdCt, frostCt, heatCt, chillPortionsCutoff, gdds_suitable)
    names(r_combined) <- c("extremeColdSuit", "springFrostSuit", "heatSuit", "chillPortionsSuit", "gddsSuit")
    r_suitable <- app(r_combined, prod)
    names(r_suitable) <- "combinedSuit"
    r_all <- c(r_suitable, r_combined)
    r_suit_hem <- paste0("r_nonlimiting_", hem)
    r_suit_all_hem <- paste0("r_nonlimiting_all_", hem)
    assign(r_suit_hem, r_suitable)
    assign(r_suit_all_hem, r_all)
    # write out hemisphere-specific suitable all files
    fileName_nonlimiting_all_hem_out <- paste0(path_data, "nonlimiting_all_", speciesName, "_", k, "_", suitabilityLevel, "_", hem, "_", yearSpan, ".tif")
    fileName_nonlimiting_all_out <- paste0(path_data, "nonlimiting_all_", speciesName, "_", k, "_", suitabilityLevel, "_", yearSpan, ".tif")
    print(system.time(writeRaster(r_all, filename = fileName_nonlimiting_all_hem_out, overwrite = TRUE, wopt = woptList)))
    flush.console()
    print(paste0("fileName_nonlimiting_all_hem_out: ", fileName_nonlimiting_all_hem_out))
  }
  r_nonlimiting_all_globe <- merge(r_nonlimiting_all_NH, r_nonlimiting_all_SH) # all the metrics
  names(r_nonlimiting_all_globe) <- names(r_nonlimiting_all_NH)
  r_nonlimiting_globe <- f_getArea(r_nonlimiting_all_globe$combinedSuit, layer = 1) # just the combined nonlimiting value
  writeRaster(r_nonlimiting_all_globe, filename = fileName_nonlimiting_all_out, overwrite = TRUE, wopt = woptList)
  flush.console()
  print(paste0("fileName_nonlimiting_all_out: ", fileName_nonlimiting_all_out))
  fileName_nonlimiting_all_df_out <- paste0(path_data, "nonlimiting_all_", speciesName, "_", k, "_", suitabilityLevel, "_", yearSpan, ".csv")
  r_nonlimiting_all_globe_df <- as.data.frame(r_nonlimiting_all_globe, xy = TRUE, na.rm = FALSE)
  fileName_nonlimiting_all_out <- paste0(path_data, "nonlimiting_all_", speciesName, "_", k, "_", suitabilityLevel, "_", yearSpan, ".tif")
  write.csv(r_nonlimiting_all_globe_df, fileName_nonlimiting_all_df_out)
  titleText <- paste0("Locations where suitability is ", suitabilityLevel, " for ", strsplit(speciesName, "_")[[1]][1], ",\nduring the ", k, " scenario", ", period ", gsub("_", "-", yearSpan))
  pal <- colorRampPalette(c("red", "green"))
  
  plot(r_nonlimiting_all_globe$combinedSuit, main = titleText, xlab = FALSE, axes = FALSE, col = pal(2))
  plot(coastline_cropped, add = TRUE)
}

f_suitableLocsPpt <- function(k, l, speciesName, suitabilityLevel) {
  yearSpan <- paste0(l, "_", l + yearRange)
  fileName_in <- paste0(path_graphics, "perennials/", speciesName, "_", suitabilityLevel, "_", k, "_", yearSpan, ".pdf")
  extImg <- external_img(src = fileName_in, width = defaultWidth, height = defaultHeight)
  my_pres <- add_slide(x = my_pres, layout = "Title Only", master = "Office Theme")
  my_pres <- ph_with(x = my_pres, value = extImg, location = ph_location(left = defaultLeft, top = defaultTop, width = defaultWidth, height = defaultHeight - 0.5))
  return(my_pres)
}

f_suitableLocsGraphics <- function(k, l, speciesName, suitabilityLevel) {
  crsRob <- "+proj=robin"
  yearSpan <- paste0(l, "_", l + yearRange)
  legendTitle <- suitabilityLevel
  speciesChoice <- paste0(speciesName, "_main")
  
  # crop out areas where summer heat or spring frost are greater than the unsuitable values, either unsuitable_springFreezeDays or unsuitable_summerHotDays
  CPfruit <- cropVals[cropName == speciesChoice, chill_portions]
  summerHeat <- cropVals[cropName == speciesChoice, summer_heat_threshold]
  cultivar <- cropVals[cropName == speciesChoice, cultivar]
  gddsFruit <- cropVals[cropName == speciesChoice, gdd]
  
  titleText <- paste0("Growing conditions for ", speciesName, ", cultivar ", cultivar, ", are ", suitabilityLevel, "\n", "during the ", k, " scenario", ", period ", gsub("_", "-", yearSpan))
  if (suitabilityLevel == "good") {
    captionString <- "Note: Locations (green) not limited by temperature for %s, cultivar %s, include at least %s chill portions, fewer than %s days of spring frost risk, \na minimum of %s growing degree days and %s days of summer heat greater than %s°C. Gray shading indicates early 21st century area \nfor all %s varieties according to data from http://www.earthstat.org. Pink shading indicates early century non-limited areas."
    caption <- sprintf(paste(captionString, collapse = " "), speciesName, cultivar, CPfruit, frostRiskDays[2], gddsFruit, heatRiskDays[2], summerHeat, speciesName)
    suitcol <- "green"
  }
  
  if (suitabilityLevel == "acceptable") {
    captionString <- "Note: Acceptable growing conditions for %s, cultivar %s, include at least %s chill portions, %s - %s  days of spring frost risk, \na minimum of %s growing degree days and \n%s - %s days of summer heat greater than %s°C. Gray shading indicates early 21st century area for all %s varieties according to data from \nhttp://www.earthstat.org. Pink shading indicates suitable areas in the beginning of the century."
    caption <- sprintf(paste(captionString, collapse = " "), speciesName, cultivar, CPfruit, frostRiskDays[3], frostRiskDays[4], gddsFruit, heatRiskDays[3], heatRiskDays[4], summerHeat, speciesName)
    suitcol <- "yellow"
  }
  
  if (suitabilityLevel == "bad") {
    captionString <- "Note: Bad growing conditions for %s, cultivar %s, include at least %s chill portions, %s - %s  days of spring frost risk, \na minimum of %s growing degree days and \n%s - %s days of summer heat greater than %s°C. Gray shading indicates early 21st century area for all %s varieties according to data from \nhttp://www.earthstat.org. Pink shading indicates suitable areas in the beginning of the century."
    caption <- sprintf(paste(captionString, collapse = " "), speciesName, cultivar, CPfruit, frostRiskDays[5], frostRiskDays[6], gddsFruit, heatRiskDays[5], heatRiskDays[6], summerHeat, speciesName)
    suitcol <- "red"
  }
  
  # this file has 6 layers. The first one is the combined suitable locations. The list of layer names - "combinedSuit"      "extremeColdSuit"   "springFrostSuit"     "heatSuit"          "chillPortionsSuit  gddsSuit"
  fileName_in <- paste0("data/perennials/nonlimiting_all_", speciesName, "_", k, "_", suitabilityLevel, "_", yearSpan, ".tif")
  r_combined <- rast(fileName_in)
  r <- r_combined[[1]]
  fileName_hist_suit_in <- paste0("data/perennials/nonlimiting_all_", speciesName, "_", "historical", "_", suitabilityLevel, "_", "1991_2010", ".tif")
  r_combined_hist <- rast(fileName_hist_suit_in)
  r_hist <- r_combined_hist[[1]]
  suitableArea_historical <- project(r_hist, crsRob)
  suitableArea_historical_df <- as.data.frame(suitableArea_historical, xy = TRUE)
  names(suitableArea_historical_df) <- c("x", "y", "value")
  suitableArea_historical_df <- round(suitableArea_historical_df, 0)
  suitableArea_historical_df[suitableArea_historical_df == 0] <- NA
  suitableArea_historical_df$value <- as.factor(suitableArea_historical_df$value)
  
  # now get harvested area map
  harvestArea_earlyCent <- f_harvestArea(speciesName, minArea = 100)
  harvestArea_earlyCent[harvestArea_earlyCent > 0] <- 1
  # convert to vector for a border and project
  harvestArea_earlyCent_rob <- project(harvestArea_earlyCent, crsRob)
  r_df <- project(r, crsRob) |>
    as.data.frame(xy = TRUE) |>
    round(0)
  names(r_df) <- c("x", "y", "value")
  r_df$value <- as.factor(r_df$value)
  harvestArea_df <- as.data.frame(harvestArea_earlyCent, xy = TRUE)
  names(harvestArea_df) <- c("x", "y", "value_harvest")
  harvestArea_df <- round(harvestArea_df, 0)
  harvestArea_df[harvestArea_df == 0] <- NA
  harvestArea_df$value_harvest <- as.factor(harvestArea_df$value_harvest)
  outF <- paste0(path_graphics, speciesName, "_", suitabilityLevel, "_", k, "_", yearSpan, ".pdf")
  
  # do without legend
  g <- ggplot() +
    geom_tile(data = r_df, aes(x, y, fill = value), show.legend = FALSE) +
    scale_fill_manual(values = c("white", "green", "grey")) +
    labs(title = titleText, fill = legendTitle, x = "", y = "", caption = caption) +
    geom_sf(data = coastline_cropped_Rob_sf, color = "black", size = 0.2) +
    theme_custom() +
    geom_tile(
      data = dplyr::filter(harvestArea_df, !is.na(value_harvest)),
      aes(x = x, y = y), fill = "grey60", alpha = .2, show.legend = FALSE
    ) +
    geom_tile(
      data = dplyr::filter(suitableArea_historical_df, !is.na(value)),
      aes(x = x, y = y), fill = "chocolate1", alpha = .2, show.legend = FALSE
    ) +
    NULL
  
  print(g)
  ggsave(filename = outF, plot = g, width = 8, height = 4, units = "in", dpi = 300)
  knitr::plot_crop(outF) # gets rid of margins around the plot
  print(paste0("file name out: ", outF))
  g <- NULL
}
# prepare suitability facet map graphics -----
f_suitableLocsGraphics_clean_old <- function(speciesName) {
  l <- 1991
  yearSpan_early <- paste0(l, "_", l + yearRange)
  l <- 2041
  yearSpan_mid <- paste0(l, "_", l + yearRange)
  l <- 2081
  yearSpan_end <- paste0(l, "_", l + yearRange)
  suitabilityLevel <- "good"
  suitcol <- "green"
    
  #  legendTitle <- suitabilityLevel
  speciesName <- gsub(chillLevel, "", speciesName) # needed for the harvested area data
  k <- "ssp585"
  fileName_in_ssp585_mid <- fileName_in <- paste0(path_data, "nonlimiting_all_", speciesName, "_", k, "_", suitabilityLevel, "_", yearSpan_mid, ".tif")
  k <- "ssp126"
  fileName_in_ssp126_mid <- fileName_in <- paste0(path_data, "nonlimiting_all_", speciesName, "_", k, "_", suitabilityLevel, "_", yearSpan_mid, ".tif")
  
  k <- "ssp585"
  fileName_in_ssp585_end <- fileName_in <- paste0(path_data, "nonlimiting_all_", speciesName, "_", k, "_", suitabilityLevel, "_", yearSpan_mid, ".tif")
  k <- "ssp126"
  fileName_in_ssp126_end <- fileName_in <- paste0(path_data, "nonlimiting_all_", speciesName, "_", k, "_", suitabilityLevel, "_", yearSpan_mid, ".tif")
  
  r_combined_ssp585_mid <- rast(fileName_in_ssp585_mid, lyrs = 1)
  r_combined_ssp126_mid <- rast(fileName_in_ssp126_mid, lyrs = 1)
  r_combined_ssp585_end <- rast(fileName_in_ssp585_end, lyrs = 1)
  r_combined_ssp126_end <- rast(fileName_in_ssp126_end, lyrs = 1)
  
  fileName_hist_suit_in <- paste0(path_data, "nonlimiting_all_", speciesName, "_", "historical", "_", suitabilityLevel, "_", yearSpan_early, ".tif")
  
  r_combined_hist <- rast(fileName_hist_suit_in, lyrs = 1)
  suitableArea_historical_rob <- project(r_combined_hist, crsRob)
  
  r_combined_ssp585_mid_rob <- project(r_combined_ssp585_mid, crsRob)
  r_combined_ssp126_mid_rob <- project(r_combined_ssp126_mid, crsRob)
  r_combined_ssp585_end_rob <- project(r_combined_ssp585_end, crsRob)
  r_combined_ssp126_end_rob <- project(r_combined_ssp126_end, crsRob)
  
  suitableArea_historical_df <- as.data.frame(suitableArea_historical_rob, xy = TRUE)
  suitableArea_historical_df$period <- "Early century"
  r_combined_ssp126_mid_rob_df <- as.data.frame(r_combined_ssp126_mid_rob, xy = TRUE)
  r_combined_ssp126_mid_rob_df$period <- "Mid century, SSP1-2.6"
  r_combined_ssp585_mid_rob_df <- as.data.frame(r_combined_ssp585_mid_rob, xy = TRUE)
  r_combined_ssp585_mid_rob_df$period <- "Mid century, SSP5-8.5"
  
  r_combined_ssp126_end_rob_df <- as.data.frame(r_combined_ssp126_end_rob, xy = TRUE)
  r_combined_ssp126_end_rob_df$period <- "End century, SSP1-2.6"
  r_combined_ssp585_end_rob_df <- as.data.frame(r_combined_ssp585_end_rob, xy = TRUE)
  r_combined_ssp585_end_rob_df$period <- "End century, SSP5-8.5"
  
  r_combined_df <- rbind(suitableArea_historical_df, r_combined_ssp126_mid_rob_df, r_combined_ssp585_mid_rob_df, r_combined_ssp126_end_rob_df, r_combined_ssp585_end_rob_df)
  names(r_combined_df) <- c("x", "y", "value", "period")
  r_combined_df$period <- factor(r_combined_df$period, levels = c("Early century", "Mid century, SSP5-8.5", "Mid century, SSP1-2.6", "End century, SSP5-8.5", "End century, SSP1-2.6"))
  
  r_combined_df$period_new <- factor(r_combined_df$period, levels = sort(c("", " ", levels(r_combined_df$period))))
  r_combined_df$period_rev <- factor(r_combined_df$period, levels = sort(unique(r_combined_df$period), decreasing = TRUE))
  
  r_combined_df$value <- round(r_combined_df$value, 0)
  r_combined_df[r_combined_df == 0] <- NA
  r_combined_df$value <- as.factor(r_combined_df$value)
  # now get harvested area map
  harvestArea_earlyCent <- f_harvestArea(speciesName, minArea = 100)
  harvestArea_earlyCent[harvestArea_earlyCent > 0] <- 1
  
  # convert to vector for a border and project
  harvestArea_earlyCent_p_sf_rob <- as.polygons(harvestArea_earlyCent) |>
    st_as_sf() |>
    st_transform(crsRob)
  
  harvestArea_df <- project(harvestArea_earlyCent, crsRob) |>
    as.data.frame(xy = TRUE) |>
    round(0)
  names(harvestArea_df) <- c("x", "y", "value_harvest")
  harvestArea_df[harvestArea_df == 0] <- NA
  
  # do without legend, title or caption
  outF <- paste0(path_graphics, "facetMaps_", speciesName, "_", suitabilityLevel, ".pdf")
  g <- ggplot() +
    geom_tile(data = r_combined_df, aes(x, y, fill = value), show.legend = FALSE, stat = "identity", position = "identity") +
    scale_fill_manual(values = c("green"), na.value = "white") +
    labs(x = "", y = "") +
    geom_sf(data = coastline_cropped_Rob_sf, color = "black", size = 0.1) +
    theme_custom() +
    geom_tile(
      data = dplyr::filter(harvestArea_df, !is.na(value_harvest)),
      aes(x = x, y = y), fill = "grey24", alpha = .2, show.legend = FALSE
    ) +
    facet_wrap(~period_rev, ncol = 2, as.table = FALSE) +
    scale_x_continuous(sec.axis = sec_axis(~., name = speciesName, breaks = NULL, labels = NULL))
  print(g)
  ggsave(filename = outF, plot = g, width = 8, height = 8, units = "in", dpi = 300)
  knitr::plot_crop(outF) # gets rid of margins around the plot
  print(paste0("file name out: ", outF))
  g <- NULL
  
  # do three period version -----
  r_combined_3period_df <- r_combined_df[r_combined_df$period %in% c("Early century", "Mid century, SSP1-2.6", "End century, SSP5-8.5"), ]
  r_combined_3period_df$period <- as.factor(r_combined_3period_df$period)
  
  outF <- paste0(path_graphics, "facetMaps_", speciesName, "_", "suitabilityLevel_3periods", ".pdf")
  g <- ggplot() +
    geom_tile(data = r_combined_3period_df, aes(x, y, fill = value), show.legend = FALSE, stat = "identity", position = "identity") +
    scale_fill_manual(values = c("green"), na.value = "white") +
    labs(x = "", y = "") +
    geom_sf(data = coastline_cropped_Rob_sf, color = "black", size = 0.1) +
    geom_sf(data = r_combined_hist_p_s_rob, color = "black", size = 0.1) +
    theme_custom() +
    geom_tile(
      data = dplyr::filter(harvestArea_df, !is.na(value_harvest)),
      aes(x = x, y = y), fill = "grey24", alpha = .2, show.legend = FALSE
    ) +
    facet_wrap(~period, ncol = 3, as.table = FALSE) +
    scale_x_continuous(sec.axis = sec_axis(~., name = speciesName, breaks = NULL, labels = NULL))
  print(g)
  ggsave(filename = outF, plot = g, width = 8, height = 8, units = "in", dpi = 300)
  knitr::plot_crop(outF) # gets rid of margins around the plot
  print(paste0("file name out: ", outF))
  g <- NULL
}

f_graphics_demoSlides <- function(r, titleText, caption, outF, col) {
  r <- project(r, crsRob, method = "near")
  r_df <- as.data.frame(r, xy = TRUE)
  names(r_df) <- c("x", "y", "value")
  r_df$value <- as.factor(r_df$value)
  
  print(paste0("file name out: ", outF))
  g <- ggplot() +
    geom_tile(data = r_df, aes(x, y, fill = value), show.legend = FALSE) +
    scale_fill_manual(values = col) +
    labs(title = titleText, fill = legendTitle, x = "", y = "", caption = caption) +
    theme_custom() +
    geom_sf(
      data = coastline_cropped_Rob_sf,
      color = "black", size = 0.1, stat = "sf", fill = NA,
      position = "identity"
    ) +
    NULL
  print(g)
  ggsave(filename = outF, plot = g, width = 8, height = 4, units = "in", dpi = 300)
  knitr::plot_crop(outF) # gets rid of margins around the plot
}

#  Growing season of ‘frost free season’ for all crops is assumed to be the period from last spring frost (defined as -2C) to first autumn frost (Tmin ≤- 2°C).
# get a years worth of data
f_yearSubset <- function(l, yearRange, r, hem) {
  yearSpan <- paste0(l, "_", l + yearRange)
  startDate <- paste0(l, "-01-01")
  endDate <- paste0(l, "-12-31") # one year of data
  if (hem == "SH") startDate <- paste0(yearnum, "-07-01")
  endDate <- paste0(yearnum + 1, "-06-30")
  indices <- seq(as.Date(startDate), as.Date(endDate), 1)
  indices <- paste0("X", as.character(indices))
  yearLayers <- subset(r, indexList)
}
