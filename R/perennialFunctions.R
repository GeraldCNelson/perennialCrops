# functions for the perennials calculations
{
  #source("R/globallyUsed.R")
  require(terra)
  library(ggplot2)
  source("R/ISIMIPconstants.R")
  source("R/ISIMIPspatialConstants.R")
  source("R/perennialsPrep.R") # creates the data tables majorCropValues_main, majorCropValues_lo and majorCropValues_hi
  library(Rcpp)
  sourceCpp("R/cpp/gdd.cpp")
  
  options(warn = 1)
  # file locations -----
  #  locOfClimFiles <- "/Volumes/ExtremeSSD2/ISIMIP/cmip6/"
  locOfClimFiles <- "climdata/"
  locOfgddsFiles <- "data/cmip6/growingDegreeDays/"
  # constants, general
  
  # choose whether to do the base cps, or the lo or hi cp requirements varieties -----
  varietiesChoices <- c("varieties_lo", "varieties_main", "varieties_hi")
  suitabilityLevels <- c("good", "acceptable", "bad")
  yearRangeSH <- 18 # one less year because of 6 month offset
  
  varietiesChoiceSuffix <- "_main" # whether to use lo, main, or hi cp values
  varietiesChoice <- paste0("varieties", varietiesChoiceSuffix)
  cropVals <- get(paste0("majorCropValues", varietiesChoiceSuffix))
  
  # constants, perennials -----
  minimumGrwSeasonLength = 100
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
  # suitability day counts - first two numbers are good window, second two numbers are acceptable window; third two numbers are bad window; fourth 2 numbers are unsuitable range
  frostRiskDays <- c(0, 6, 7, 20, 21, 45, 46, springFrostLength) 
  heatRiskDays <- c(0, 12, 13, 30, 31, 45, 46, heatDamageLength) 
  
  
  # variables that define the runs parameters. This code is used to create the runs values in runsHemisphere.R Copying it here is a kludge. Need to think of a better way to deal with the need for this in f_gddSums
  {
    climVal <- -2 # changed from 0 because subject experts say plants can tolerate this
    test_logic <- paste0("x > ", climVal)
    logicDirection <- ">"
    if (logicDirection == ">") ldtext <-"gt"
    if (logicDirection == "<") ldtext <-"lt"
    runlengthChoices <- c(100) # at least 100 days of tmin > 0
    climateVariable <- "tasmin"
    runsParms <- c(climVal, test_logic, logicDirection, ldtext, runlengthChoices, climateVariable)
  }
  #test values, perennials -----
  speciesChoices <- unique(cropVals$cropName)
  speciesChoice <- "cherry_main"
  suitabilityLevel <- "good"
  climateVarChoice <- "tasmin"
  hem <- "NH"
  
  # functions -----
  f_range <- function(x, rangeVals) {
    # set locations with values outside the range to 999 as an indicator where land is for later adjustment
    # x[x < rangeVals[1]] <- NA
    # x[x > rangeVals[2]] <- NA
    x[x < rangeVals[1]] <- 999
    x[x > rangeVals[length(rangeVals)]] <- 999 # length needed because rangeVals can sometimes have more than 2 entries.
    return(x)
  }
  
  f_convertToSHdays <- function(yearNum, calIn) {
    if (lubridate::leap_year(yearNum)) {dayCt <- 366}else{dayCt <- 365}
    NHYear <- rep(1:dayCt, 1)
    SHyrStart <- rep(182:dayCt,1) # july 1
    SHyrmid <- rep(1:181,1)
    SHYear <- c(SHyrStart, SHyrmid)
    daysLookup <- data.frame(cbind(NHYear, SHYear))
    uniqueVals <- unique(calIn)
    for (val in 1:length(uniqueVals)) {
      calIn[calIn == uniqueVals[val]] <- daysLookup[uniqueVals[val], "SHYear"]
    }
    return(calIn)
  }
  
  f_readRast_count <- function(modelChoice, climateVarChoice, threshold, layersToKeep, probVal, k, l, hem) {
    print(climateVarChoice)
    print(threshold)
    modelChoice_lower <- tolower(modelChoice)
    yearSpan <- paste0(l, "_", l + yearRange)
    fileName_hem_in <- paste0(locOfClimFiles, modelChoice_lower, "_", climateVarChoice, "_", k, "_", hem, "_", yearSpan, ".tif")
    print(paste0("climate fileName_hem_in: ", fileName_hem_in))
    r <- rast(fileName_hem_in)
    #    print(r)
    # convert day numbers to calendar days to subset
    # indices values are from 1 to the number of years in the subsetted file, usually 20.
    # layersToKeep are the layer numbers to keep in each year.
    datesToKeep <- c()
    indices <- c()
    for (yearNum in l:(l + yearRange)) {
      temp <- as.Date(layersToKeep, origin = paste0(yearNum, "-01-01"))
      indices_yearNum <- rep(yearNum - l + 1, length(layersToKeep))
      indices <-c(indices, indices_yearNum)
      temp <- paste0("X", temp)
      datesToKeep <- c(datesToKeep, temp)
    }
    r <- subset(r, datesToKeep)
    # for spring frost damage
    if (climateVarChoice == "tasmin") f_ct <- function(x) (sum(x < threshold)) 
    if (climateVarChoice == "tasmax") f_ct <- function(x) (sum(x > threshold)) 
    # print(paste0("length r: ", length(names(r)), ", length indices: ", length(indices)))
    
    print(system.time(tempCt <- tapp(r, indices, f_ct)))
    print(tempCt)
    return(tempCt)
  }
  
  f_readRast_extreme <- function(modelChoice_lower, climateVarChoice, funDir, hem) {
    yearSpan <- paste0(l, "_", l + yearRange)
    fileName_hem_in <- paste0(locOfClimFiles, modelChoice_lower, "_", climateVarChoice, "_", k, "_", hem, "_", yearSpan, ".tif")
    # fileName_hem_in <- paste0(locOfClimFiles, "mean_daily/", "mean_daily", "_", climateVarChoice, "_", modelChoice_lower, "_", k, "_", hem, "_", yearSpan, ".tif") 
    print(paste0("climate fileName_hem_in in: ", fileName_hem_in))
    r <- rast(fileName_hem_in)
    print(r)
    indices <- format(as.Date(names(r), format = "X%Y-%m-%d"), format = "%Y") # %Y is year as 4 digit number
    indices <- as.numeric(indices)
    indices <- indices - l  + 1
    print(system.time(extremeTemp <- tapp(r, indices, funDir))) #indices are from 1 to the number of years in the input file. 365 1s, then 365 2s, etc. The funDir is the minimum or maximum value of the temp var. extremeTemp is the highest or lowest temp value in each year in each cell.
    print(extremeTemp)
    return(extremeTemp)
  }
  
  f_extremeCold <- function(k, l, speciesChoice, hem, varietiesChoiceSuffix) {
    yearSpan <- paste0(l, "_", l + yearRange)
    climateVarChoice <- "tasmin"
    funDir <- "min"
    probVal <- 0.80
    speciesName <- gsub(varietiesChoiceSuffix, "", speciesChoice) # needed for the harvested area data
    system.time(x <- lapply(modelChoices_lower, f_readRast_extreme, climateVarChoice, funDir, hem)) # read in tasmin for the relevant period and all ESMs
    r <- rast(x)
    r
    extremeColdCutoff <- cropVals[cropName == speciesChoice, low_temp_threshold] 
    
    # now do ensemble mean and cutoff
    print(paste0("Working on extreme cold for speciesChoice: ", speciesChoice, ", working on ssp: ", k, ", start year ", l, ", hemisphere ", hem))
    fileName_lo_out <- paste0(locOfDataFiles_perennials, "extremeCold_cutoff_", speciesName, "_lo", "_", k, "_", hem, "_", yearSpan, ".tif")
    fileName_main_out <- paste0(locOfDataFiles_perennials, "extremeCold_cutoff_", speciesName, "_main", "_", k, "_", hem, "_", yearSpan, ".tif")
    fileName_high_out <- paste0(locOfDataFiles_perennials, "extremeCold_cutoff_", speciesName, "_hi", "_", k, "_", hem, "_", yearSpan, ".tif")
    print(system.time(extremeCold_quant <- quantile(r, probs = probVal, na.rm = TRUE)))
    extremeCold_quant[extremeCold_quant < extremeColdCutoff] <- 0 #  extreme cold limited
    extremeCold_quant[extremeCold_quant >= extremeColdCutoff] <- 1 # not extreme cold limited
    # rangeVals <- c(extremeColdCutoff, 50)
    # extremeCold_quant <- f_range(extremeCold_quant, rangeVals) 
    # extremeCold_quant[extremeCold_quant >= extremeColdCutoff & extremeCold_quant < 999] <- 1 # use of 999 here and in f_range is a kludge to set non-temp-limited land areas to zero.
    # extremeCold_quant[extremeCold_quant == 999] <- 0
    # 
    print(system.time(writeRaster(extremeCold_quant, filename = fileName_lo_out, overwrite = TRUE, wopt = woptList)))
    print(paste0("fileName_out extreme cold: ", fileName_lo_out))
    file.copy(from = fileName_lo_out, to = fileName_main_out, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
    file.copy(from = fileName_lo_out, to = fileName_high_out, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
    
    plot(extremeCold_quant, main = speciesChoice)
    return(extremeCold_quant)
  }
  
  f_computeGDDs <- function(k, l, modelChoice, cropVals) {
    yearSpan <- paste0(l, "_", l + yearRange)
    print(paste0("start year: ", l, ", ssp: ", k,  " model: ", modelChoice, ", start year: ", l, ", hemisphere: ", hem))
    modelChoice_lower <- tolower(modelChoice)
    yearSpan <- paste0(l, "_", l + yearRange)
    fileName_tas_in <- paste0(locOfClimFiles, modelChoice_lower, "_", "tas", "_", k, "_", yearSpan, ".tif")
    tas <- rast(fileName_tas_in)
    print(tas)
    
    # split tas up into individual years and run gdd on those
    # with July 6, 2021 adjustments gddtb and GDD_opt are identical for two groups - almond, apple and cherry in one and olive and winegrape in the other
    # first do gdds for one from each group - choose almond and winegrape
    for (speciesChoice in c(paste0("almond", varietiesChoiceSuffix), paste0("winegrape", varietiesChoiceSuffix))) {
      fileName_out <- paste0(locOfgddsFiles, modelChoice_lower, "_", "gdd", "_", speciesChoice, "_", k, "_", yearSpan, ".tif")
      #      if (!fileName_out %in% gddFilesCompleted) { commented out temporarily to make sure all are redone
      print(paste0("Working on: ", fileName_out))
      topt_min <- cropVals[cropName == speciesChoice, gddtb]
      topt_max <- cropVals[cropName == speciesChoice, GDD_opt]
      print(paste0("crop: ", speciesChoice, " topt_min: ", topt_min, " topt_max: ", topt_max, " fileName_out: ", fileName_out))
      print(system.time(gdd <- app(tas, fun = f_gdd, topt_min = topt_min, topt_max = topt_max, cores = 1, filename = fileName_out, overwrite = TRUE, wopt = woptList)))
      print(paste0("gdd file out name: ", fileName_out))
      if (speciesChoice == paste0("almond", varietiesChoiceSuffix)) {
        print("copying almond to apple and cherry")
        for (i in c("apple", "cherry")) {
          speciesChoice_alt <- paste0(i, varietiesChoiceSuffix)
          fileName_out_alt <- paste0(locOfgddsFiles, modelChoice_lower, "_", "gdd", "_", speciesChoice_alt, "_", k, "_", yearSpan, ".tif")
          file.copy(fileName_out, fileName_out_alt, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
        }
      }
      
      if (speciesChoice == paste0("winegrape", varietiesChoiceSuffix)) {
        for (i in c("olive")) {
          print("copying winegrape to olive")
          speciesChoice_alt <- paste0(i, varietiesChoiceSuffix)
          fileName_out_alt <- paste0(locOfgddsFiles, modelChoice_lower, "_", "gdd", "_", speciesChoice_alt, "_", k, "_", yearSpan, ".tif")
          file.copy(fileName_out, fileName_out_alt, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
        }
      }
      #   
    }
    #    gdd <- NULL
    #    gc()
    # #  }else{
    # #    print(paste("This file has already been created: ", fileName_out))
    # #  }
    #  # copy gdd files to other crops with same gdd ranges
    #  gddFilesCompleted <- list.files(locOfgddsFiles,  full.names = TRUE)
    #  gddFilesCompleted <- gddFilesCompleted[!grepl("aux.xml", gddFilesCompleted, fixed = TRUE)]
    #  gddFilesCompleted <- gsub("//", "/", gddFilesCompleted)
    #  gddFilesCompleted <- gddFilesCompleted[!grepl("gddSum_mean_", gddFilesCompleted, fixed = TRUE)]
    #  gddFilesCompleted <- gddFilesCompleted[!grepl("ensemble_", gddFilesCompleted, fixed = TRUE)]
    #  gddFilesCompleted <- gddFilesCompleted[!grepl("gddSum_", gddFilesCompleted, fixed = TRUE)]
    #  
    #  gddFiles_almond <- gddFilesCompleted[grepl("almond", gddFilesCompleted, fixed = TRUE)]
    #  gddFiles_winegrape <- gddFilesCompleted[grepl("winegrape", gddFilesCompleted, fixed = TRUE)]
    #  
    #  for (i in gddFiles_almond) {
    #    print(i)
    #    gdd_cherry <- gsub("almond", "cherry", i)
    #    print(gdd_cherry)
    #    file.copy(from = i, to = gdd_cherry, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
    #    gdd_apple <- gsub("almond", "apple", i)
    #    file.copy(from = i, to = gdd_apple, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
    #  }
    #  gdd_olive <- gsub("winegrape", "olive", gddFiles_winegrape)
    #  file.copy(from = gddFiles_winegrape, to = gdd_olive, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
  }
  
  f_gddSums <- function(k, l, speciesChoice, hem, runsParms) {
    yearSpan <- paste0(l, "_", l + yearRange)
    logicDirection <- runsParms[3]
    climVal <- runsParms[1]
    climateVariable <- runsParms[6]
    ldtext <- runsParms[4]
    runlength <- runsParms[5]
    for (modelChoice in modelChoices) {
      modelChoice_lower <- tolower(modelChoice)
      fileName_gdd_in <- paste0(locOfgddsFiles, modelChoice_lower, "_", "gdd", "_", speciesChoice, "_", k, "_", yearSpan, ".tif")
      gdds <- rast(fileName_gdd_in)
      gdds <- crop(gdds, get(paste0("extent_", hem)))
      # gdds are daily for the 20 year period
      # if (hem == "SH")  {startDate <-  paste0(l, "-07-01"); endDate <- paste0(l + yearRange-1, "-06-30")} # in southern hemisphere search July 1 to June 30 of the next year. NH is just the calendar year
      # if (hem == "NH")  {startDate <-  paste0(l, "-01-01"); endDate <- paste0(l + yearRange, "-12-31")} # in southern hemisphere search July 1 to June 30 of the next year. NH is just the calendar year
      startDate <-  paste0(l, "-01-01"); endDate <- paste0(l + yearRange, "-12-31")
      indices <- seq(as.Date(startDate), as.Date(endDate), by = "days")
      indicesChar <- paste0("X", indices)
      # in case gdds doesn't have correct names
      names(gdds) <- indicesChar
      indicesYr <- unique(as.numeric(format(indices, "%Y")))
      if (hem == "SH") indicesYr <- indicesYr[1:yearRange]
      if (!nlyr(gdds) == length(indicesChar)) gdds <- subset(gdds, indicesChar) # if SH, gets rid of the first 1/2 year and  last 1/2 year. may not be necessary because I think gdds in sh file already have this done. If statement may capture this
      fileName_startDay1_in <- paste0(locOfRunsFiles, "startday_1_", climateVariable,"_",modelChoice_lower, "_run_", runlength, "_lim_", ldtext, climVal, "_", hem, "_", k, "_", yearSpan, ".tif")
      fileName_endDay1_in <- paste0(locOfRunsFiles, "endday_1_", climateVariable,"_",modelChoice_lower, "_run_", runlength, "_lim_", ldtext, climVal, "_", hem, "_", k, "_", yearSpan, ".tif")
      startDay <- rast(fileName_startDay1_in)
      endDay <- rast(fileName_endDay1_in)
      names(startDay) <-names(endDay) <- sort(unique(indicesYr))
      
      # now do calc by year
      for (yearNumber in 1:nlyr(startDay)) {
        print(paste0("Working on species choice: ", speciesChoice, ", ssp: " , k, ", startYear: ", l, ", hem: ", hem, ", model: ", modelChoice, ", yearNumber: ", yearNumber))
        startDay_yr <- subset(startDay, yearNumber)
        endDay_yr <- subset(endDay, yearNumber)
        startYear <- l + yearNumber - 1
        
        if (hem == "SH")  {
          startDate <-  paste0(startYear, "-07-01"); endDate <- paste0(startYear + 1, "-06-30")} # in southern hemisphere search July 1 to June 30 of the next year. NH is just the calendar year
        if (hem == "NH")  {
          startDate <-  paste0(startYear, "-01-01"); endDate <- paste0(startYear, "-12-31")
        }
        indices <- seq(as.Date(startDate), as.Date(endDate), by = "days")
        indicesChar <- paste0("X", indices)
        #            print(system.time(sum_gdds <- app(gdds_yr, f_sumVec, startDay_yr, endDay_yr)))
        if ((hem == "SH" & yearNumber < 20) | (hem == "NH")) {
          gdds_yr <- subset(gdds, indicesChar)
        }
        # rapp needs to have a start day that is greater than 0. Next two lines sets all 0 values to 1
        startDay_yr[startDay_yr == 0] <- 1
        endDay_yr[endDay_yr == 0] <- 1
        endDay_yr[endDay_yr > nlyr(gdds_yr)] <- nlyr(gdds_yr)  # needed because the end day can be # 367 if we're in a tropical region; end day is the frost day that ends the run
        print(system.time(sum_gdds <- rapp(gdds_yr, startDay_yr, endDay_yr, "sum")))
        plot(sum_gdds, main = paste0("Sum of gdds in a run of at least 100 days, scenario: ", k, ", period: ", yearSpan, ", year number: ", yearNumber, ", model: ", modelChoice_lower, ", crop: ", speciesChoice))
        if (yearNumber == 1 ) {
          period_sums <- sum_gdds
        } else {
          period_sums <- c(period_sums, sum_gdds)
        }
      }
      gc()
      period_sums
      fileName_gddSums_out <- paste0(locOfgddsFiles, "gddSum", "_", modelChoice_lower, "_", hem, "_",  speciesChoice, "_", k, "_", yearSpan, ".tif")
      print(system.time(writeRaster(period_sums, filename = fileName_gddSums_out,  overwrite = TRUE, wopt= woptList))); flush.console()
      print(paste0("fileName_gddSums_out: ", fileName_gddSums_out))
    }
  }
  
  f_gddSum_mean <- function (k, l, speciesChoices, modelChoices, hemispheres) {
    yearSpan <- paste0(l, "_", l + yearRange)
    for (speciesChoice in speciesChoices) {
      for (modelChoice in modelChoices) {
        modelChoice_lower <- tolower(modelChoice)
        for (hem in hemispheres) {
          fileName_gddSums_in <- paste0(locOfgddsFiles, "gddSum", "_", modelChoice_lower, "_", hem, "_",  speciesChoice, "_", k, "_", yearSpan, ".tif")
          r_in <- rast(fileName_gddSums_in)
          fileName_out = paste0(locOfgddsFiles, "gddSum_mean", "_", modelChoice_lower, "_", hem, "_",  speciesChoice, "_", k, "_", yearSpan, ".tif")
          test <- app(r_in, mean, filename = fileName_out, overwrite = TRUE, wopt = woptList)
          print(paste0("fileName_out: ", fileName_out))
        }
      }
    }
  }
  
  f_readRast_gddSum <- function(modelChoice, speciesChoice, k, l, hem) {
    yearSpan <- paste0(l, "_", l + yearRange)
    modelChoice_lower <- tolower(modelChoice)
    fileName_in = paste0(locOfgddsFiles, "gddSum_mean", "_", modelChoice_lower, "_", hem, "_",  speciesChoice, "_", k, "_", yearSpan, ".tif")
    print(paste0("speciesChoice: ", speciesChoice, ", k: ", k, ", modelChoice: ", modelChoice, ", fileName in: ", fileName_in))
    r <- rast(fileName_in)
  }
  
  f_ensemble_GDD_sum_mean <- function(k, l, yearRange, speciesChoices, cropVals) {
    yearSpan <- paste0(l, "_", l + yearRange)
    for (speciesChoice in speciesChoices) {
      gddsRequired <- cropVals[cropName == speciesChoice, gdd]
      for (hem in hemispheres) {
        x <- lapply(modelChoices, f_readRast_gddSum, speciesChoice, k, l, hem)
        r <- rast(x)
        indices_day <- rep(seq(1, nlyr(x[[1]]), 1), 5) # 5 is number of models; if omitted should get the same result
        fileName_out <- paste0(locOfgddsFiles, "ensemble_gddSum_mean", "_", hem, "_",  speciesChoice, "_", k, "_", yearSpan, ".tif")
        print(paste0("speciesChoice: ", speciesChoice, ", ensemble ssp: ", k, ", start year: ", l , ", fileName out: ", fileName_out))
        print(system.time(r_mean <- tapp(r, indices_day, fun = "mean", na.rm = TRUE, filename = fileName_out, overwrite = TRUE, wopt = woptList)))
        main <- paste0("Perennial: ", speciesChoice, ", GDDs required: ", gddsRequired, ", hemisphere: ", hem, ", ssp: ", k, ", period: ", yearSpan)
        plot(r_mean, main = main, axes = FALSE)
      }
    }
  }
  
  f_gddsSuitability <- function(k, l, speciesChoice, hem) {
    yearSpan <- paste0(l, "_", l + yearRange)
    # get gdds ensemble mean
    fileName_in <- paste0(locOfgddsFiles, "ensemble_gddSum_mean", "_", hem, "_",  speciesChoice, "_", k, "_", yearSpan, ".tif")
    gdds <- rast(fileName_in)                     
    # get gdd requirements
    gddsRequired <- cropVals[cropName == speciesChoice, gdd]
    gddsSuitable <- gdds
    gddsSuitable[gddsSuitable <  gddsRequired] <- 0 
    gddsSuitable[gddsSuitable >=  gddsRequired] <- 1 
    fileName_out <- paste0(locOfgddsFiles, "gdds_not_limiting", "_", hem, "_",  speciesChoice, "_", k, "_", yearSpan, ".tif")
    print(system.time(writeRaster(gddsSuitable, filename = fileName_out, overwrite = TRUE, wopt = woptList)))
    print(paste0("fileName_out suitable gdds: ", fileName_out))
    
    plot(gddsSuitable, main = paste0("Adequate GDDs for ", speciesChoice, ", minimum required ", gddsRequired, ", hemisphere: ", hem, ", scenario: ", k, ", period: ", yearSpan))
    return(gddsSuitable)
  }
  
  # chillPortions calcs done in chillPortions.R -----
  
  f_frostDamage <- function(k, l, speciesChoice, hem, varietiesChoice, suitabilityLevel) {
    yearSpan <- paste0(l, "_", l + yearRange)
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
    system.time(x <- lapply(modelChoices_lower, f_readRast_count, climateVarChoice, threshold, layersToKeep, probVal, k, l, hem))
    r <- rast(x)
    r[r <= frDays[2]] <- 1
    r[r > frDays[2]] <- 0
    #       fr <- f_range(r, frDays)
    system.time(fr <- quantile(r, probs = probVal, na.rm = TRUE)) # note: if all layers have the same value quantile returns that value
    fr[fr > 0] <- 1
    
    titleText_fr <- paste0("Green indicates locations where frost days are not limiting for ", strsplit(speciesChoice, "_")[[1]][1], " during the ", k, " scenario", ", ", gsub("_", "-", yearSpan), ". \nUnsuitable frost risk is more than ", frDays[2], " frost days (-2°C) during the spring frost window.")
    plot(fr, main = titleText_fr)
    fileName_fr_out <- paste0(locOfDataFiles_perennials, "frostDamage_", speciesChoice, "_", k, "_", suitabilityLevel, "_", hem, "_", yearSpan, ".tif")
    writeRaster(fr, fileName_fr_out, overwrite = TRUE, wopt = woptList)
    print(paste0(" frost damage fileName out: ", fileName_fr_out))
    print("--------------------------------------------")
  }
  
  f_heatDamage <- function(k, l, speciesChoice, hem, varietiesChoice, suitabilityLevel) {
    yearSpan <- paste0(l, "_", l + yearRange)
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
    system.time(x <- lapply(modelChoices_lower, f_readRast_count, climateVarChoice, threshold, layersToKeep, probVal, k, l, hem))
    r <- rast(x)
    r[r <= hdDays[2]] <- 1
    r[r > hdDays[2]] <- 0
    #        hd <- f_range(r, hdDays)
    system.time(hd <- quantile(r, probs = probVal, na.rm = TRUE)) # note: if all layers have the same value quantile returns that value
    hd[hd > 0] <- 1
    #  r[r > 0] <- 1
    titleText_hd <- paste0("1 indicates locations where extreme summer heat is not limiting for ", strsplit(speciesChoice, "_")[[1]][1], " during the ", k, " scenario", ", ", gsub("_", "-", yearSpan), ". \nUnsuitable heat is more than ", hdDays[2], " days above ", threshold, "°C during the summer window.")
    
    plot(hd, main = titleText_hd)
    fileName_hd_out <- paste0(locOfDataFiles_perennials, "heatDamage_", speciesChoice, "_", k, "_", suitabilityLevel, "_", hem, "_", yearSpan, ".tif")
    writeRaster(hd, fileName_hd_out, overwrite = TRUE, wopt = woptList)
    print(paste0(" heat damage fileName out: ", fileName_hd_out))
    print("--------------------------------------------")
  }
  
  f_combinedDamage <- function(k, l, speciesChoice, suitabilityLevel) {
    # combine suitability metrics from chill portions, extreme cold, spring frost, and summer heat; locations with value 1 is suitable
    yearSpan <- paste0(l, "_", l + yearRange)
    for (hem in hemispheres) {
      # read in all the rasters needed
      fileName_extremeColdCt_in <- paste0(locOfDataFiles_perennials, "extremeCold_cutoff_", speciesChoice, "_", k, "_", hem, "_", yearSpan, ".tif")
      extremeColdCt <- rast(fileName_extremeColdCt_in) 
      fileName_fr_in <- paste0(locOfDataFiles_perennials, "frostDamage_", speciesChoice, "_", k, "_", suitabilityLevel, "_", hem, "_", yearSpan, ".tif") 
      frostCt <- rast(fileName_fr_in) 
      fileName_hd_in <- paste0(locOfDataFiles_perennials, "heatDamage_", speciesChoice, "_", k, "_", suitabilityLevel, "_", hem, "_", yearSpan, ".tif") 
      heatCt <- rast(fileName_hd_in) 
      fileName_cp_in <- paste0("data/cmip6/chillPortions/chill_portions/", "ensemble_chill_cutoff_", speciesChoice, "_", k, "_", hem, "_", yearSpan, ".tif")
      chillPortionsCutoff <- rast(fileName_cp_in) 
      fileName_gdds_in <- paste0(locOfgddsFiles, "gdds_not_limiting", "_", hem, "_",  speciesChoice, "_", k, "_", yearSpan, ".tif")
      gdds_suitable <- rast(fileName_gdds_in)
      
      print(paste0("working on combined damage ", speciesChoice, " in hemisphere ", hem, ", year ", l, ", scenario ", k))
      
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
      fileName_nonlimiting_all_hem_out <- paste0(locOfDataFiles_perennials, "nonlimiting_all_", speciesChoice, "_",  k, "_", suitabilityLevel, "_", hem, "_", yearSpan, ".tif")
      fileName_nonlimiting_all_out <- paste0(locOfDataFiles_perennials, "nonlimiting_all_", speciesChoice, "_",  k, "_", suitabilityLevel, "_", yearSpan, ".tif")
      print(system.time(writeRaster(r_all, filename = fileName_nonlimiting_all_hem_out,  overwrite = TRUE,  wopt= woptList))); flush.console()
    }
    r_nonlimiting_globe <- merge(r_nonlimiting_NH, r_nonlimiting_SH) # just the combined nonlimiting value
    r_nonlimiting_all_globe <- merge(r_nonlimiting_all_NH, r_nonlimiting_all_SH)
    print(system.time(writeRaster(r_nonlimiting_all_globe, filename = fileName_nonlimiting_all_out,  overwrite = TRUE,  wopt= woptList))); flush.console()
    fileName_nonlimiting_all_df_out <- paste0(locOfDataFiles_perennials, "nonlimiting_all_", speciesChoice, "_",  k, "_", suitabilityLevel, "_", yearSpan, ".csv")
    r_nonlimiting_all_globe_df <- as.data.frame(r_nonlimiting_all_globe, xy = TRUE, na.rm = FALSE)
    fileName_nonlimiting_all_out <- paste0(locOfDataFiles_perennials, "nonlimiting_all_", speciesChoice, "_",  k, "_", suitabilityLevel, "_", yearSpan, ".tif")
    write.csv(r_nonlimiting_all_globe_df, fileName_nonlimiting_all_df_out)
    titleText <- paste0("Locations where suitability is ", suitabilityLevel, " for ", strsplit(speciesChoice, "_")[[1]][1], ", during the ", k, " scenario", ", period ", gsub("_", "-", yearSpan))
    pal <- colorRampPalette(c("red", "green"))
    
    plot(r_nonlimiting_globe, main = titleText, legend = FALSE, xlab = FALSE, axes=FALSE, col = pal(2))
    plot(coastline_cropped, add = TRUE)
  }
  
  f_suitableLocsPpt <- function(k, l, yearSpan, suitabilityLevel, speciesChoice) {
    yearSpan <- paste0(l, "_", l + yearRange)
    fileName_in <- paste0(lofOfGraphicsFiles, "perennials/", speciesChoice, "_", suitabilityLevel, "_", k, "_", yearSpan, ".png")
    extImg <- external_img(src = fileName_in, width = defaultWidth, height = defaultHeight)
    my_pres <- add_slide(x = my_pres, layout = 'Title Only', master = 'Office Theme')
    my_pres <- ph_with(x = my_pres, value = extImg, location = ph_location(left = defaultLeft, top = defaultTop, width = defaultWidth, height = defaultHeight - 0.5))
    return(my_pres)
  }
  
  f_suitableLocsGraphics <- function(k, l, speciesChoice, suitabilityLevel) {
    yearSpan <- paste0(l, "_", l + yearRange)
    legendTitle <- suitabilityLevel
    speciesName <- gsub(varietiesChoiceSuffix, "", speciesChoice) # needed for the harvested area data
    # the harvested area data is just for grapes so need to get rid of wine in the names
    if (speciesChoice == paste0("winegrape", varietiesChoiceSuffix)) speciesName <- "grape"
    
    # crop out areas where summer heat or spring frost are greater than the unsuitable values, either unsuitable_springFreezeDays or unsuitable_summerHotDays
    CPfruit <- cropVals[cropName == speciesChoice, chill_portions]
    summerHeat <- cropVals[cropName == speciesChoice, summer_heat_threshold]
    cultivar <-  cropVals[cropName == speciesChoice, cultivar]
    gddsFruit <- cropVals[cropName == speciesChoice, gdd]
    
    titleText <- paste0("Growing conditions for ", speciesName,", cultivar ", cultivar, ", are ", suitabilityLevel, "\n" , "during the ", k, " scenario", ", period ", gsub("_", "-", yearSpan))
    if (suitabilityLevel == "good") {
      # caption <- paste0("Note: Good growing conditions for ", speciesName, ", cultivar ", cultivar, " include chill portions of at least ", CPfruit, ", less than ", frostRiskDays[3], " days of spring frost risk and \nless than ", heatRiskDays[3], " days of summer heat greater than ", summerHeat )#, "°C. Gray chading indicates early 21st century area.")
      captionString <- "Note: Locations (green) not limited by temperature for %s, cultivar %s, include at least %s chill portions, fewer than %s days of spring frost risk, \na minimum of %s growing degree days and %s days of summer heat greater than %s°C. Gray shading indicates early 21st century area \nfor all %s varieties according to data from http://www.earthstat.org. Pink shading indicates early century non-limited areas."
      caption <- sprintf(paste(captionString, collapse = " ") , speciesName, cultivar, CPfruit, frostRiskDays[2], gddsFruit, heatRiskDays[2], summerHeat, speciesName)
      suitcol = "green"}
    
    if (suitabilityLevel == "acceptable") {
      captionString <-   "Note: Acceptable growing conditions for %s, cultivar %s, include at least %s chill portions, %s - %s  days of spring frost risk, \na minimum of %s growing degree days and \n%s - %s days of summer heat greater than %s°C. Gray shading indicates early 21st century area for all %s varieties according to data from \nhttp://www.earthstat.org. Pink shading indicates suitable areas in the beginning of the century."
      caption <- sprintf(paste(captionString, collapse = " ") , speciesName, cultivar, CPfruit, frostRiskDays[3], frostRiskDays[4], gddsFruit, heatRiskDays[3], heatRiskDays[4], summerHeat, speciesName)
      suitcol = "yellow"}
    
    if (suitabilityLevel == "bad") {
      captionString <-   "Note: Bad growing conditions for %s, cultivar %s, include at least %s chill portions, %s - %s  days of spring frost risk, \na minimum of %s growing degree days and \n%s - %s days of summer heat greater than %s°C. Gray shading indicates early 21st century area for all %s varieties according to data from \nhttp://www.earthstat.org. Pink shading indicates suitable areas in the beginning of the century."
      caption <- sprintf(paste(captionString, collapse = " ") , speciesName, cultivar, CPfruit, frostRiskDays[5], frostRiskDays[6], gddsFruit, heatRiskDays[5], heatRiskDays[6], summerHeat, speciesName)
      suitcol = "red"}
    
    # this file has 6 layers. The first one is the combined suitable locations. The list of layer names - "combinedSuit"      "extremeColdSuit"   "springFrostSuit"     "heatSuit"          "chillPortionsSuit  gddsSuit"
    fileName_in <- paste0("data/cmip6/perennials/nonlimiting_all_", speciesChoice, "_", k, "_", suitabilityLevel, "_", yearSpan, ".tif")
    r_combined <- rast(fileName_in)
    r <- r_combined[[1]]
    # if not historical, get historical suitability for use as background shading
    #    if (!k == "historical") {
    fileName_hist_suit_in <- paste0("data/cmip6/perennials/nonlimiting_all_", speciesChoice, "_", "historical", "_", suitabilityLevel, "_", "1991_2010", ".tif")
    r_combined_hist <- rast(fileName_hist_suit_in)
    r_hist <- r_combined_hist[[1]]
    suitableArea_historical <- project(r_hist, crsRob)
    suitableArea_historical_df <- as.data.frame(suitableArea_historical, xy = TRUE)
    names(suitableArea_historical_df) <- c("x", "y", "value")
    suitableArea_historical_df <- round(suitableArea_historical_df, 0)
    suitableArea_historical_df[suitableArea_historical_df == 0] <- NA
    suitableArea_historical_df$value = as.factor(suitableArea_historical_df$value)
    #    }
    # now get harvested area map
    rInArea <- rast(paste0(locOfHarvestDataFiles, speciesName,"/",  speciesName, "_HarvestedAreaHectares.tif"))
    harvestArea_earlyCent <- aggregate(rInArea, fact = 6, fun = "sum") # convert 5 arc minutes to 1/2 degrees
    maskMin <- switch(
      speciesName,
      "almond" = 100,
      "apple" = 100,
      "cherry" = 100,
      "grape" = 100,
      "olive" = 100
    )
    harvestArea_earlyCent[harvestArea_earlyCent < maskMin] <- 0 # set minimum area to be greater than 50 hectares per grid cell
    harvestArea_earlyCent[harvestArea_earlyCent > 0] <- 1
    r <- project(r, crsRob)
    harvestArea_earlyCent <- project(harvestArea_earlyCent, crsRob)
    r_df <- as.data.frame(r, xy = TRUE)
    names(r_df) <-   c("x", "y", "value")
    r_df <- round(r_df, 0)
    r_df$value = as.factor(r_df$value)
    harvestArea_df <- as.data.frame(harvestArea_earlyCent, xy = TRUE)
    names(harvestArea_df) <-   c("x", "y", "value_harvest")
    harvestArea_df <- round(harvestArea_df, 0)
    harvestArea_df[harvestArea_df == 0] <- NA
    harvestArea_df$value_harvest = as.factor(harvestArea_df$value_harvest)
    fileName_out <- paste0(lofOfGraphicsFiles, "perennials/", speciesChoice, "_", suitabilityLevel, "_", k, "_", yearSpan, ".png")
    
    # do without legend
    g <- ggplot() +
      geom_tile(data = r_df, aes(x, y, fill = value), show.legend = FALSE) +
      #      geom_raster(data = r_df, aes(x, y, fill = value)) +
      scale_fill_manual(values = c("white", "green", "grey")) +
      labs(title = titleText, fill = legendTitle, x = "", y = "", caption = caption) + 
      
      geom_sf(data = coastline_cropped_Rob_sf,  color="black", size = 0.2) +
      #     coord_sf(xlim=c(3,35), ylim=c(52,72)) + 
      theme_bw() +
      theme(
        legend.text.align = 1,
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 12, hjust = 0.5),
        plot.caption = element_text(hjust = 0, vjust = 7.0, size = 7)
      ) +
      #      geom_tile(data = harvestArea_df, aes(x, y), fill = "gray", alpha = .2, show.legend = FALSE)
      geom_tile(data = dplyr::filter(harvestArea_df, !is.na(value_harvest)), 
                aes(x = x, y = y), fill = "grey60", alpha = .2, show.legend = FALSE) +
      geom_tile(data = dplyr::filter(suitableArea_historical_df, !is.na(value)), 
                aes(x = x, y = y), fill = "chocolate1", alpha = .2, show.legend = FALSE) +
      NULL
    
    print(g)
    ggsave(filename = fileName_out, plot = g, width = 8, height = 4, units = "in", dpi = 300)
    knitr::plot_crop(fileName_out) # gets rid of margins around the plot
    print(paste0("file name out: ", fileName_out))
    g <- NULL
  }
  # prepare suitability facet map graphics -----
  f_suitableLocsGraphics_clean <- function(speciesChoice) {
    l <- 1991; yearSpan_early <- paste0(l, "_", l + yearRange)
    l <- 2041; yearSpan_mid <- paste0(l, "_", l + yearRange)
    l <- 2081; yearSpan_end  <- paste0(l, "_", l + yearRange)
    suitabilityLevel <- "good"
    suitcol = "green"
    
    #  legendTitle <- suitabilityLevel
    speciesName <- gsub(varietiesChoiceSuffix, "", speciesChoice) # needed for the harvested area data
    # the harvested area data is just for grapes so need to get rid of wine in the names
    if (speciesChoice == paste0("winegrape", varietiesChoiceSuffix)) speciesName <- "grape"
    
    k <- "ssp585"; fileName_in_ssp585_mid <- fileName_in <- paste0("data/cmip6/perennials/nonlimiting_all_", speciesChoice, "_", k, "_", suitabilityLevel, "_", yearSpan_mid, ".tif")
    k <- "ssp126"; fileName_in_ssp126_mid <- fileName_in <- paste0("data/cmip6/perennials/nonlimiting_all_", speciesChoice, "_", k, "_", suitabilityLevel, "_", yearSpan_mid, ".tif")
    
    k <- "ssp585"; fileName_in_ssp585_end <- fileName_in <- paste0("data/cmip6/perennials/nonlimiting_all_", speciesChoice, "_", k, "_", suitabilityLevel, "_", yearSpan_mid, ".tif")
    k <- "ssp126"; fileName_in_ssp126_end <- fileName_in <- paste0("data/cmip6/perennials/nonlimiting_all_", speciesChoice, "_", k, "_", suitabilityLevel, "_", yearSpan_mid, ".tif")
    
    r_combined_ssp585_mid <- rast(fileName_in_ssp585_mid)[[1]]
    r_combined_ssp126_mid <- rast(fileName_in_ssp126_mid)[[1]]
    r_combined_ssp585_end <- rast(fileName_in_ssp585_end)[[1]]
    r_combined_ssp126_end <- rast(fileName_in_ssp126_end)[[1]]
    
    fileName_hist_suit_in <- paste0("data/cmip6/perennials/nonlimiting_all_", speciesChoice, "_", "historical", "_", suitabilityLevel, "_", yearSpan_early, ".tif")
    r_combined_hist <- rast(fileName_hist_suit_in)[[1]]
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
    
    #   r_combined_df <- rbind(suitableArea_historical_df, r_combined_ssp126_mid_rob_df, r_combined_ssp585_mid_rob_df)
    r_combined_df <- rbind(suitableArea_historical_df, r_combined_ssp126_mid_rob_df, r_combined_ssp585_mid_rob_df, r_combined_ssp126_end_rob_df, r_combined_ssp585_end_rob_df)
    names(r_combined_df) <- c("x", "y", "value", "period")
   # r_combined_df$period <- factor(r_combined_df$period, levels = c("Early century", "Mid century, ssp126", "Mid century, ssp126", "End century, ssp585", "End century, ssp126"))
     r_combined_df$period <- factor(r_combined_df$period, levels = c("Early century", "Mid century, SSP5-8.5", "Mid century, SSP1-2.6", "End century, SSP5-8.5", "End century, SSP1-2.6"))
     
    r_combined_df$period_new <- factor(r_combined_df$period, levels=sort(c(""," ",levels(r_combined_df$period))))
    r_combined_df$period_rev <- factor(r_combined_df$period, levels=sort(unique(r_combined_df$period), decreasing=TRUE))
    
    r_combined_df$value <- round(r_combined_df$value, 0)
    r_combined_df[r_combined_df == 0] <- NA
    r_combined_df$value = as.factor(r_combined_df$value)
    # now get harvested area map
    rInArea <- rast(paste0(locOfHarvestDataFiles, speciesName,"/",  speciesName, "_HarvestedAreaHectares.tif"))
    harvestArea_earlyCent <- aggregate(rInArea, fact = 6, fun = "sum") # convert 5 arc minutes to 1/2 degrees
    maskMin <- switch(
      speciesName,
      "almond" = 100,
      "apple" = 100,
      "cherry" = 100,
      "grape" = 100,
      "olive" = 100
    )
    harvestArea_earlyCent[harvestArea_earlyCent < maskMin] <- 0 # set minimum area to be greater than 50 hectares per grid cell
    harvestArea_earlyCent[harvestArea_earlyCent > 0] <- 1
    harvestArea_earlyCent <- project(harvestArea_earlyCent, crsRob)
    harvestArea_df <- as.data.frame(harvestArea_earlyCent, xy = TRUE)
    names(harvestArea_df) <-   c("x", "y", "value_harvest")
    harvestArea_df <- round(harvestArea_df, 0)
    harvestArea_df[harvestArea_df == 0] <- NA
    
    # do without legend, title or caption
    fileName_out <- paste0(lofOfGraphicsFiles, "perennials/facetMaps_", speciesChoice, "_", suitabilityLevel, ".png")
    g <- ggplot() +
      geom_tile(data = r_combined_df, aes(x, y, fill = value), show.legend = FALSE, stat = "identity", position = "identity") +
      scale_fill_manual(values = c("green"), na.value = 'white') +
      #   labs(title = titleText, fill = legendTitle, x = "", y = "", caption = caption) + 
      labs( x = "", y = "") +
      geom_sf(data = coastline_cropped_Rob_sf,  color="black", size = 0.1) +
      theme_bw() +
      theme(
        legend.text.align = 1,
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 12, hjust = 0.5),
        plot.caption = element_text(hjust = 0, vjust = 7.0, size = 7)
      ) +
      #      geom_tile(data = harvestArea_df, aes(x, y), fill = "gray", alpha = .2, show.legend = FALSE)
      geom_tile(data = dplyr::filter(harvestArea_df, !is.na(value_harvest)), 
                aes(x = x, y = y), fill = "grey24", alpha = .2, show.legend = FALSE) +
      # geom_tile(data = dplyr::filter(suitableArea_historical_df, !is.na(value)), 
      #           aes(x = x, y = y), fill = "chocolate1", alpha = .2, show.legend = FALSE) +
      facet_wrap(~period_rev, ncol = 2,  as.table=FALSE) +
      scale_x_continuous(sec.axis = sec_axis(~ . , name = speciesName, breaks = NULL, labels = NULL))
    print(g)
    ggsave(filename = fileName_out, plot = g, width = 8, height = 8, units = "in", dpi = 300)
    knitr::plot_crop(fileName_out) # gets rid of margins around the plot
    print(paste0("file name out: ", fileName_out))
    g <- NULL
    
    # do three period version -----
    r_combined_3period_df <- r_combined_df[r_combined_df$period %in% c("Early century",  "Mid century, SSP1-2.6", "End century, SSP5-8.5"),]
    r_combined_3period_df$period <- as.factor(r_combined_3period_df$period)
    
    fileName_out <- paste0(lofOfGraphicsFiles, "perennials/facetMaps_", speciesChoice, "_", "suitabilityLevel_3periods", ".png")
    g <- ggplot() +
      geom_tile(data = r_combined_3period_df, aes(x, y, fill = value), show.legend = FALSE, stat = "identity", position = "identity") +
      scale_fill_manual(values = c("green"), na.value = 'white') +
      #   labs(title = titleText, fill = legendTitle, x = "", y = "", caption = caption) + 
      labs( x = "", y = "") +
      geom_sf(data = coastline_cropped_Rob_sf,  color="black", size = 0.1) +
      theme_bw() +
      theme(
        legend.text.align = 1,
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 12, hjust = 0.5),
        plot.caption = element_text(hjust = 0, vjust = 7.0, size = 7)
      ) +
      #      geom_tile(data = harvestArea_df, aes(x, y), fill = "gray", alpha = .2, show.legend = FALSE)
      geom_tile(data = dplyr::filter(harvestArea_df, !is.na(value_harvest)), 
                aes(x = x, y = y), fill = "grey24", alpha = .2, show.legend = FALSE) +
      # geom_tile(data = dplyr::filter(suitableArea_historical_df, !is.na(value)), 
      #           aes(x = x, y = y), fill = "chocolate1", alpha = .2, show.legend = FALSE) +
      facet_wrap(~period, ncol = 3,  as.table=FALSE) +
      scale_x_continuous(sec.axis = sec_axis(~ . , name = speciesName, breaks = NULL, labels = NULL))
    print(g)
    ggsave(filename = fileName_out, plot = g, width = 8, height = 8, units = "in", dpi = 300)
    knitr::plot_crop(fileName_out) # gets rid of margins around the plot
    print(paste0("file name out: ", fileName_out))
    g <- NULL
  }
  
  # used to prepare suitability facet map graphics
  f_suitableLocsGraphics_compare_chillPortions <- function(speciesChoice) {
    #  l <- 1991; yearSpan_early <- paste0(l, "_", l + yearRange)
    l <- 2081; yearSpan_end <- paste0(l, "_", l + yearRange)
    suitabilityLevel <- "good"
    suitcol = "green"
    
    #  legedTitle <- suitabilityLevel
    speciesName <- gsub(varietiesChoiceSuffix, "", speciesChoice) 
    # the harvested area data is just for grapes so need to get rid of wine in the names
    if (speciesChoice == paste0("winegrape", varietiesChoiceSuffix)) speciesName <- "grape"
    
    #main
    speciesChoice_main <- paste0(speciesName, "_main")
    k <- "ssp585"; fileName_in_ssp585_end <- fileName_in <- paste0("data/cmip6/perennials/nonlimiting_all_", speciesChoice_main, "_", k, "_", suitabilityLevel, "_", yearSpan_end, ".tif")
    k <- "ssp126"; fileName_in_126_end <- fileName_in <- paste0("data/cmip6/perennials/nonlimiting_all_", speciesChoice_main, "_", k, "_", suitabilityLevel, "_", yearSpan_end, ".tif")
    r_combined_ssp585_end_main <- rast(fileName_in_ssp585_end)[[1]]
    r_combined_ssp126_end_main <- rast(fileName_in_ssp585_end)[[1]]
    
    #lo
    speciesChoice_lo <- paste0(speciesName, "_lo")
    k <- "ssp585"; fileName_in_ssp585_end_lo <- fileName_in <- paste0("data/cmip6/perennials/nonlimiting_all_", speciesChoice_lo, "_", k, "_", suitabilityLevel, "_", yearSpan_end, ".tif")
    k <- "ssp126"; fileName_in_126_end_lo <- fileName_in <- paste0("data/cmip6/perennials/nonlimiting_all_", speciesChoice_lo, "_", k, "_", suitabilityLevel, "_", yearSpan_end, ".tif")
    r_combined_ssp585_end_lo <- rast(fileName_in_ssp585_end)[[1]]
    r_combined_ssp126_end_lo <- rast(fileName_in_ssp585_end)[[1]]
    
    r_combined_ssp585_end_main_rob <- project(r_combined_ssp585_end_main, crsRob)
    r_combined_ssp126_end_main_rob <- project(r_combined_ssp126_end_main, crsRob)
    r_combined_ssp126_end_main_rob_df <- as.data.frame(r_combined_ssp126_end_main_rob, xy = TRUE)
    r_combined_ssp126_end_main_rob_df$period <- "End century, SSP1-2.6"
    r_combined_ssp585_end_main_rob_df <- as.data.frame(r_combined_ssp585_end_main_rob, xy = TRUE)
    r_combined_ssp585_end_main_rob_df$period <- "End century, SSP5-8.5"
    
    r_combined_ssp585_end_lo_rob <- project(r_combined_ssp585_end_lo, crsRob)
    r_combined_ssp126_end_lo_rob <- project(r_combined_ssp126_end_lo, crsRob)
    r_combined_ssp126_end_lo_rob_df <- as.data.frame(r_combined_ssp126_end_lo_rob, xy = TRUE)
    r_combined_ssp126_end_lo_rob_df$period <- "End century, SSP1-2.6"
    r_combined_ssp585_end_lo_rob_df <- as.data.frame(r_combined_ssp585_end_lo_rob, xy = TRUE)
    r_combined_ssp585_end_lo_rob_df$period <- "End century, SSP5-8.5"
    
    r_combined_df <- rbind(r_combined_ssp126_end_main_rob_df, r_combined_ssp585_end_main_rob_df)
    names(r_combined_df) <- c("x", "y", "value", "period")
    r_combined_df$period <- as.factor(r_combined_df$period)
    r_combined_df$value <- round(r_combined_df$value, 0)
    r_combined_df[r_combined_df == 0] <- NA
    r_combined_df$value = as.factor(r_combined_df$value)
    
    # now get harvested area map
    rInArea <- rast(paste0(locOfHarvestDataFiles, speciesName,"/",  speciesName, "_HarvestedAreaHectares.tif"))
    harvestArea_earlyCent <- aggregate(rInArea, fact = 6, fun = "sum") # convert 5 arc minutes to 1/2 degrees
    maskMin <- switch(
      speciesName,
      "almond" = 100,
      "apple" = 100,
      "cherry" = 100,
      "grape" = 100,
      "olive" = 100
    )
    harvestArea_earlyCent[harvestArea_earlyCent < maskMin] <- 0 # set minimum area to be greater than 50 hectares per grid cell
    harvestArea_earlyCent[harvestArea_earlyCent > 0] <- 1
    harvestArea_earlyCent <- project(harvestArea_earlyCent, crsRob)
    harvestArea_df <- as.data.frame(harvestArea_earlyCent, xy = TRUE)
    names(harvestArea_df) <-   c("x", "y", "value_harvest")
    harvestArea_df <- round(harvestArea_df, 0)
    harvestArea_df[harvestArea_df == 0] <- NA
    harvestArea_df$value_harvest = as.factor(harvestArea_df$value_harvest)
    
    # do without legend, title or caption
    fileName_out <- paste0(lofOfGraphicsFiles, "perennials/", speciesChoice, "_", suitabilityLevel, "_", k, "_", yearSpan, ".png")
    g <- ggplot() +
      geom_tile(data = r_combined_df, aes(x, y, fill = value), show.legend = FALSE, stat = "identity", position = "identity") +
      scale_fill_manual(values = c("green"), na.value = 'white') +
      #   labs(title = titleText, fill = legendTitle, x = "", y = "", caption = caption) + 
      labs( x = "", y = "") +
      geom_sf(data = coastline_cropped_Rob_sf,  color="black", size = 0.2) +
      theme_bw() +
      theme(
        legend.text.align = 1,
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 12, hjust = 0.5),
        plot.caption = element_text(hjust = 0, vjust = 7.0, size = 7)
      ) +
      #      geom_tile(data = harvestArea_df, aes(x, y), fill = "gray", alpha = .2, show.legend = FALSE)
      geom_tile(data = dplyr::filter(harvestArea_df, !is.na(value_harvest)), 
                aes(x = x, y = y), fill = "grey24", alpha = .2, show.legend = FALSE) +
      # geom_tile(data = dplyr::filter(suitableArea_historical_df, !is.na(value)), 
      #           aes(x = x, y = y), fill = "chocolate1", alpha = .2, show.legend = FALSE) +
      facet_wrap(~period, ncol = 2) +
      scale_x_continuous(sec.axis = sec_axis(~ . , name = speciesName, breaks = NULL, labels = NULL))
    NULL
    
    print(g)
    ggsave(filename = fileName_out, plot = g, width = 8, height = 4, units = "in", dpi = 300)
    knitr::plot_crop(fileName_out) # gets rid of margins around the plot
    print(paste0("file name out: ", fileName_out))
    g <- NULL
  }
  
  f_graphics_demoSlides <- function(r, titleText, caption, fileName_out, col) {
    r <- project(r, crsRob, method = "near")
    r_df <- as.data.frame(r, xy = TRUE)
    names(r_df) <-   c("x", "y", "value")
    r_df$value <- as.factor(r_df$value)
    
    print(paste0("file name out: ", fileName_out))
    g <- ggplot() +
      geom_tile(data = r_df, aes(x, y,  fill = value), show.legend = FALSE) +
      #              scale_fill_discrete(colors = col, drop = FALSE, na.value = 'grey95') +
      scale_fill_manual(values = col) +
      # the na.value doesn't work yet. See https://stackoverflow.com/questions/45144630/scale-fill-manual-define-color-for-na-values/45147172 for a possible solution
      labs(title = titleText, fill = legendTitle, x = "", y = "", caption = caption) + 
      theme_bw()  +
      #      theme(plot.title = element_text(size = 12, hjust = 0.5), plot.caption = element_text(hjust = 0, size = 8)) +
      theme(
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 12, hjust = 0.5),
        plot.caption = element_text(hjust = 0, vjust = 7.0, size = 7)
      ) +
      geom_sf(data = coastline_cropped_Rob_sf , 
              color = "black", size = 0.1, stat = "sf", fill = NA,
              position = "identity") +
      NULL
    print(g)
    ggsave(filename = fileName_out, plot = g, width = 8, height = 4, units = "in", dpi = 300)
    knitr::plot_crop(fileName_out) # gets rid of margins around the plot
  }
  
  #  Growing season of ‘frost free season’ for all crops is assumed to be the period from last spring frost (defined as -2C) to first autumn frost (Tmin ≤- 2°C).
  # get a years worth of data
  f_yearSubset <- function(l, yearRange, r, hem) {
    yearSpan <- paste0(l, "_", l + yearRange)
    startDate <- paste0(l, "-01-01"); endDate <- paste0(l, "-12-31") # one year of data
    if (hem == "SH") startDate <-  paste0(yearnum, "-07-01"); endDate <- paste0(yearnum + 1, "-06-30")
    indices <- seq(as.Date(startDate), as.Date(endDate), 1)
    indices <- paste0("X", as.character(indices))
    yearLayers <- subset(r, indexList)
  }
}
