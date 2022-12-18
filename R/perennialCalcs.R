# code to do various perennial crop calculations

speciesChoices <- c("almond_main", "apple_main",  "cherry_main", "olive_main",  "grape_main")
speciesNames <- gsub("_main", "", speciesChoices)
chillLevels <- c("_lo", "_main")
source("R/perennialFunctions.R")
source("R/perennialsPrep.R")
suitabilityLevel <- "good"
# the code in {} below runs all the raw data crunching - extreme cold, heat, frost, GDDs - except CPs which is done in ; the code below that does aggregation - means by model and ensemble means. Don't run these unless the files need to be updated and you have some time to kill

# Note: that the difference between _lo and _main is just in the chill portions for the crops

{
  # extreme cold calcs -----
  #scenarios
  cropVals <- get(paste0("majorCropValues", "_main"))
  for (hem in hemispheres) {
    for (speciesName in speciesNames) {
      for (k in sspChoices) {
        for (l in startYearChoices) {
          f_extremeCold(k, l, speciesName, hem, modelChoices_lower, cropVals)
        }
      }
    
    #historical -----
    k <- "historical"
    l <- 1991
    f_extremeCold(k, l, speciesName, hem, modelChoices_lower, cropVals)
    }
  }
  
  # heat damage  -----
  # scenarios
  cropVals <- get(paste0("majorCropValues", "_main"))
  for (hem in hemispheres) {
    for (speciesName in speciesNames) {
      for (k in sspChoices) {
        for (l in startYearChoices) {
          f_heatDamage(k, l, speciesName, hem, suitabilityLevel, cropVals)
        }
      }
      #historical
      k <- "historical"
      l <- 1991
      f_heatDamage(k, l, speciesName, hem, suitabilityLevel, cropVals)
    }
  }
  
  # frost damage -----
  #scenarios
  for (hem in hemispheres) {
    for (speciesName in speciesNames) {
      for (k in sspChoices) {
        for (l in startYearChoices) {
          f_frostDamage(k, l, speciesName, hem, suitabilityLevel, cropVals)
        }
      }
    }
    #historical
    k <- "historical"
    l <- 1991
    f_frostDamage(k, l, speciesName, hem, suitabilityLevel, cropVals)
  }
  
  # gdd calcs ------
  # note that gdds don't change on the basis of species chill portions, so generate just with file names
  #scenarios
  for (speciesName in speciesNames) {
    speciesChoice <- paste0(speciesName, "_main") # needed to get a cropVals variable with the gdds, _lo would also work as well
    cropVals <- get(paste0("majorCropValues", "_main"))
    topt_min <- cropVals[cropName == speciesChoice, gddtb]
    topt_max <- cropVals[cropName == speciesChoice, GDD_opt]
    for (modelChoice in modelChoices) {
      for (k in sspChoices) {
        for (l in startYearChoices) {
          print(system.time(f_computeGDDs(k, l, speciesName, modelChoice, topt_min, topt_max)))
        }
      }
      
      # historical
      k <- "historical"
      l <- 1991
      print(system.time(f_computeGDDs(k, l, speciesName, modelChoice, topt_min, topt_max)))
    }
  }
}
# runs files -----
# runsParms created in perennialFunctions.R
#scenarios
for (k in sspChoices) {
  for (l in startYearChoices) {
    f_runsSetup(k, l, runsParms)
  }
}
#historical
k = "historical"
l = 1991
f_runsSetup(k, l, runsParms)

# gdd sums  ------ this takes forever becuase it works year by year.
#scenarios
for (speciesName in speciesNames) {
  for (modelChoice in modelChoices) {
    for (hem in hemispheres) {
      for (k in sspChoices) {
        for (l in startYearChoices) {
          f_gddSums(k, l, speciesName, hem, runsParms, modelChoice)
        }
      }
      #historical
      k <- "historical"
      l <- 1991
      f_gddSums(k, l, speciesName, hem, runsParms, modelChoice)
    }
  }
}

# aggregation code -----
# runs only used in GDD calcs so no need to do means
{
  # GDD sum, means by model -----
  #  scenarios
  cropVals <- get(paste0("majorCropValues", "_main"))
  for (k in sspChoices) {
    for (l in startYearChoices) {
      f_gddSum_mean(k, l, speciesNames, modelChoices, hemispheres)
    }
  }
  #historical -----
  k <- "historical"
  l <- 1991
  f_gddSum_mean(k, l, speciesNames, modelChoices, hemispheres)
  
  
  # ensemble GDD sum -----
  #scenarios 
  for (k in sspChoices) {
    for (l in startYearChoices) {
      yearSpan <- paste0(l, "_", l + yearRange)
      f_ensemble_GDD_sum_mean(k, l, yearRange, speciesNames, cropVals)
    }
  }
  #historical
  k <- "historical"
  l <- 1991
  f_ensemble_GDD_sum_mean(k, l, yearRange, speciesNames, cropVals)
}

# gdds suitable-----

for (hem in hemispheres) {
  for (speciesName in speciesNames) {
    #scenarios
    for (k in sspChoices) {
      for (l in startYearChoices) {
        print(paste0("speciesName: ", speciesName, ", ssp choice: ", k, ", start year: ", l))
        f_gddsSuitability(k, l, speciesName, hem, cropVals) 
      }
    }
    #historical
    k <- "historical"
    l <- 1991
    print(paste0("speciesName: ", speciesName, ", ssp choice: ", k, ", start year: ", l))
    f_gddsSuitability(k, l, speciesName, hem, cropVals) 
  }
}

# combined damage, scenarios -----
# code to read in 1/0 metrics files and produce 1/0 tifs where the crop is potentially growable. The chill portions files are created in the chillPortions.R script
#Important note: The chillPortions.R script must be run before combined damage whenever any chill portion value is changed.
{
  # combined damage -----
  for (chillLevel in chillLevels) {
    for (speciesName in speciesNames) {
      speciesName <- f_speciesName(speciesName)
      for (k in sspChoices) {
        for (l in startYearChoices) {
          print(paste0("speciesName: ", speciesName, ", chilllevel: ", chillLevel,  ", ssp choice: ", k, ", start year: ", l))
          f_combinedDamage(k, l, speciesName, suitabilityLevel, chillLevel) 
        }
      }
      #historical
      k <- "historical"
      l <- 1991
      print(paste0("speciesName: ", speciesName, ", chilllevel: ", chillLevel,  ", ssp choice: ", k, ", start year: ", l))
      f_combinedDamage(k, l, speciesName, suitabilityLevel, chillLevel) 
    } 
  }
}

#suitable locations graphics, historical -----
for (hem in hemispheres) {
  for (speciesName in speciesNames) {
    k <- "historical"
    l <- 1991
    print(paste0("Suitability level:  ", suitabilityLevel, ", speciesName: ", speciesName, ", ssp choice: ", k, ", start year: ", l))
    f_suitableLocsGraphics(k, l, speciesName, suitabilityLevel)
    
    #suitable locations graphics, scenarios -----
    for (k in sspChoices) {
      for (l in startYearChoices) {
        print(paste0("Suitability level:  ", suitabilityLevel, ", speciesName: ", speciesName, ", ssp choice: ", k, ", start year: ", l))
        f_suitableLocsGraphics(k, l, speciesName, suitabilityLevel)
      }
    }
  }
}

# moved to facet maps code to GCBFig2.R

# #create facet maps for good suitability -----
# for (speciesName in speciesNames) {
#   f_suitableLocsGraphics_clean(speciesName)
# }

# demo slides ------
k = "historical"
l = 1991
yearSpan <- paste0(l, "_", l + yearRange)
speciesName <- "cherry_main"
legendTitle <- "Suitable"
defaultColor <- c("green", "red")
CPfruit <- cropVals[cropName == speciesName, chill_portions]
gddsFruit <- cropVals[cropName == speciesName, gdd]
excessiveHeat <- cropVals[cropName == speciesName, summer_heat_threshold]

caption <- paste0("Note: ", suitabilityLevel, " growing conditions for ", strsplit(speciesName, "_")[[1]][1], ", ", cropVals[cropName == speciesName, cultivar], " variety", " include chill portions of at least ", CPfruit, ", at least ", gddsFruit, " growing degree days, \nno more than ", frostRiskDays[2], " days of frost risk during spring and no more than ", heatRiskDays[2], " days of excessive heat (above ", excessiveHeat,  "°C) in summer.")

fileName_nonlimiting_all_NH_in <- paste0(locOfDataFiles_perennials, "nonlimiting_all_", speciesName, "_",  k, "_", suitabilityLevel, "_", "NH", "_", yearSpan, ".tif")
suitable_combined_NH <- rast(fileName_nonlimiting_all_NH_in)
suitable_NH <- subset(suitable_combined_NH, "combinedSuit") 
extremeCold_NH <- subset(suitable_combined_NH, "extremeColdSuit") 
unsuitableFrost_NH <- subset(suitable_combined_NH, "springFrostSuit")
unsuitableHeat_NH <- subset(suitable_combined_NH, "heatSuit")
chillPortionsCutoff_NH <- subset(suitable_combined_NH, "chillPortionsSuit")
gddsSuitable_NH <- subset(suitable_combined_NH, "gddsSuit") 

fileName_nonlimiting_all_SH_in <- paste0(locOfDataFiles_perennials, "nonlimiting_all_", speciesName, "_",  k, "_", suitabilityLevel, "_", "SH", "_", yearSpan, ".tif")
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
tminExtremeVal <- cropVals[cropName == speciesName, low_temp_threshold]
titleText <- paste0("Temperatures below ", tminExtremeVal, "°C for at least one day can kill ", strsplit(speciesName, "_")[[1]][1], ". \nResults are for the ", k, " scenario, ", "period ", gsub("_", "-", yearSpan), ". \nGreen areas are not limited by this metric.")
fileName_out <- paste0(lofOfGraphicsFiles, "perennials/demoSlide_extremeColdMask_", speciesName, "_", suitabilityLevel, "_", k, "_", yearSpan, ".png")
f_graphics_demoSlides(r, titleText, caption, fileName_out, col)

r <- chillPortionsCutoff # 1 is where chill portions is adequate 
col <- rev(defaultColor)
titleText <- paste0("Locations where chill portions are adequate for ", strsplit(speciesName, "_")[[1]][1], ", ", cropVals[cropName == speciesName, cultivar], " variety", ". \nResults are for the ", k, " scenario, ", "period ", gsub("_", "-", yearSpan))
fileName_out <- paste0(lofOfGraphicsFiles, "perennials/demoSlide_chillPortionsCutoff_", speciesName, "_", suitabilityLevel, "_", k, "_", yearSpan, ".png")
f_graphics_demoSlides(r, titleText, caption, fileName_out, col)

{
  # these next lines of code must stay in this order.
  # titleText <- "Locations where the chill portions are not affected by extreme cold"
  # fileName_out <- paste0(lofOfGraphicsFiles, "perennials/demoSlide_chillPortionsGood", speciesName, "_", suitabilityLevel, "_", k, "_", yearSpan, ".png")
  # f_graphics_demoSlides(r, titleText, caption, fileName_out, col)
  
  r <- unsuitableFrost
  col <- rev(defaultColor)
  titleText = paste0("Locations where the frost days during spring are not limiting for ", strsplit(speciesName, "_")[[1]][1], ". \nResults are for the ", k, " scenario, ", "period ", gsub("_", "-", yearSpan))
  fileName_out <- paste0(lofOfGraphicsFiles, "perennials/demoSlide_springFrostRisk_", speciesName, "_", suitabilityLevel, "_", k, "_", yearSpan, ".png")
  f_graphics_demoSlides(r, titleText, caption, fileName_out, col)
  
  r <- unsuitableHeat # 1 is where heat Ct is above the unsuitable level
  col <- rev(defaultColor)
  titleText = paste0("Locations where summer hot days for ",  strsplit(speciesName, "_")[[1]][1], " are not limiting.", "\nResults are for the ", k, " scenario, ", "period ", gsub("_", "-", yearSpan))
  fileName_out <- paste0(lofOfGraphicsFiles, "perennials/demoSlide_summerHeat_", speciesName, "_", suitabilityLevel, "_", k, "_", yearSpan, ".png")
  f_graphics_demoSlides(r, titleText, caption, fileName_out, col)
  
  r <- gddsSuitable # 1 is where heat Ct is above the unsuitable level
  col <- rev(defaultColor)
  titleText = paste0("Locations where the growing degree days requirement for ",  strsplit(speciesName, "_")[[1]][1], " is met.", "\nResults are for the ", k, " scenario, ", "period ", gsub("_", "-", yearSpan))
  fileName_out <- paste0(lofOfGraphicsFiles, "perennials/demoSlide_gdds_", speciesName, "_", suitabilityLevel, "_", k, "_", yearSpan, ".png")
  f_graphics_demoSlides(r, titleText, caption, fileName_out, col)
}

# r <- heatCt
# col <- defaultColor
# titleText = "Locations in red are where the number of summer hot days are not suitable."
# fileName_out <- paste0(lofOfGraphicsFiles, "perennials/demoSlide_summerHeat_", speciesName, "_", suitabilityLevel, "_", k, "_", yearSpan, ".png")
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
  for (speciesName in speciesNames) {
    ensembleTitle <- paste("Non-temperature-limited locations for ", strsplit(speciesName, "_")[[1]][1]) #, ", Locations quality: ", suitabilityLevel
    my_pres <- add_slide(x = my_pres, layout = 'Section Header', master = 'Office Theme')
    my_pres <- ph_with(x = my_pres, value = ensembleTitle, location = ph_location_type(type = "title"))
    
    k <- "historical"
    l <- 1991
    my_pres <- f_suitableLocsPpt(k, l, speciesName, suitabilityLevel)
    
    for (k in sspChoices) {
      for (l in startYearChoices) {
        my_pres <- f_suitableLocsPpt(k, l, speciesName, suitabilityLevel)
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
speciesName <- "cherry_main"
suitabilityLevel = "good"
cropLocations <- as.data.table(read_excel("data-raw/perennials/cropLocations.xlsx"))
#locsVector <- cropLocations[, location := NULL]
#speciesName <- gsub("_main", "", speciesName)
speciesName <- f_speciesName(speciesName)
#if (speciesName == "winegrape_main") speciesName <- "wine grape"
locsVector <- cropLocations[species == speciesName,]
setnames(locsVector, old = c("lon", "lat"), new = c("x", "y"))
fileName_nonlimiting_all_in <- paste0(locOfDataFiles_perennials, "nonlimiting_all_", speciesName, "_",  k, "_", suitabilityLevel, "_", yearSpan, ".tif")
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
text(subset(labs, labs$islon), lab = parse(text = labs$lab[labs$islon]), pos = 1, cex = 0.75, srt = 90)
text(subset(labs, !labs$islon), lab = parse(text = labs$lab[!labs$islon]), pos = 2, cex = 0.75)
par(op)
points(x = locsVector$x, y = locsVector$y)


# old code -----
# the special code commented out just below is (hopefully replaced by running f_gddSums over all crops)
# added because almond, apple, and cherry have the same gdd requirements and grape and olive also have the same
# first get all the gdd file names and pull out just the ones for almond and grape
# gddFilesCompleted <- list.files(locOfgddsFiles,  full.names = TRUE)
# gddFilesCompleted <- gddFilesCompleted[!grepl("aux.xml", gddFilesCompleted, fixed = TRUE)]
# gddFilesCompleted <- gsub("//", "/", gddFilesCompleted)
# gddSumFilesCompleted <- gddFilesCompleted[!grepl("gddSum_mean_", gddFilesCompleted, fixed = TRUE)]
# gddSumFilesCompleted <- gddSumFilesCompleted[!grepl("ensemble_", gddSumFilesCompleted, fixed = TRUE)]
# gddSumFilesCompleted <- gddSumFilesCompleted[grepl("gddSum_", gddSumFilesCompleted, fixed = TRUE)]
# 
# gddSumFiles_almond <- gddSumFilesCompleted[grepl("almond", gddSumFilesCompleted, fixed = TRUE)]
# gddSumFiles_grape <- gddSumFilesCompleted[grepl("grape", gddSumFilesCompleted, fixed = TRUE)]
# 
# # because the gdds are the same for main and lo, the next bit of code just creates a new  filename and copies the main to the new
# 
# for (i in gddSumFiles_almond) {
#   gdd_cherry <- gsub("almond", "cherry", i)
#   
#   file.copy(from = i, to = gdd_cherry, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
#   gdd_apple <- gsub("almond", "apple", i)
#   file.copy(from = i, to = gdd_apple, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
#   print(gdd_cherry)
#   print(i)
#   al <- rast(i)
#   plot(al, 1)
#   al <- rast(gdd_cherry)
#   plot(al, 1)
#   
# }
# for (i in gddSumFiles_grape) {
#   gdd_olive <- gsub("grape", "olive", i)
#   file.copy(from = i, to = gdd_olive, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
# }
# 
# for ( i in gddSumFilesCompleted) {
#   i_lo <- gsub("_main", "_lo", i) 
#   file.copy(from = i, to = i_lo, overwrite = TRUE, copy.mode = TRUE)
# }
