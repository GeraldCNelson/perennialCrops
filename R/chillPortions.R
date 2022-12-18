# calculate 1-0 values for locations where chill portions are enough, by fruit
{
  require(terra)
  require(sf)
  library("crayon") # for red color in cat
  require(ggplot2)
  source("R/ISIMIPconstants.R")
  source("R/ISIMIPspatialConstants.R")
  source("R/perennialsPrep.R") # get the latest chill portions data
  # choose whether to do the base chill portions, or the lo chill protion requirements
  varietyChoices <- c("varieties_main")
  chillLevels <- c("_lo", "_main")  
  #test values
  modelChoice <- "UKESM1-0-LL"
  k <- "ssp585"
  l <- 2041
  midYear <- 2080
  hem <- "NH"
  speciesChoice <- "cherry_lo"
  chillLevel <- "_lo"
  
  # choice for varietyChoice in next line-----------
  cropVals <- get(paste0("majorCropValues", chillLevel))
  speciesChoices <- sort(unique(cropVals$cropName))
  
  f_readRast_quantile <- function(modelChoice, k, l, hem) {
    modelChoice_lower <- tolower(modelChoice)
    midYear <- l + 9
    if (hem == "NH") hem_full <- "north"
    if (hem == "SH") hem_full <- "south"
    fileName_in_hem <- paste0(locOfCPFiles, k,"/", modelChoice_lower, "/", k, "_", modelChoice_lower, "_", midYear, "_chill_portions_", hem_full, ".tif")
    print(paste0("fileName in: ", fileName_in_hem))
    r <- rast(fileName_in_hem)
    if (hem == "NH") r <- subset(r, 1:19) # because NH chill season crosses end of calendar year
    system.time(chillPortion <- quantile(r, probs = 0.2, na.rm = TRUE))
    print(chillPortion)
    return(chillPortion)
  }
  
  f_chillPortions <- function(k, l, speciesChoice, hem, chillLevel) {
    cropVals <- get(paste0("majorCropValues", chillLevel))
    midYear <- l + 9
    yearSpan <- paste0(l, "_", l + yearRange)
    ext_hem <- get(paste0("extent_", hem))
    system.time(x <- lapply(modelChoices_lower, f_readRast_quantile, k, l, hem))
    r <- rast(x)
    r
    print(system.time(r.mean <- app(r, fun = "mean", na.rm = TRUE)))
    print(r.mean)
    #    browser()
    fileNameCP_out <- paste0(locOfCPFiles, "ensemble_chill_portions", "_", k, "_", hem, "_", yearSpan, ".tif")
    writeRaster(r.mean, filename = fileNameCP_out, overwrite = TRUE, wopt = woptList)
    # now do ensemble 1-0 calcs
    r.mean_copy <-r.mean
    cplimit <- cropVals[cropName == speciesChoice, chill_portions]
    print(paste0("working on ssp: ", k, ", start year ", l, ", hemisphere ", hem, ", crop ", speciesChoice, ", cplimit: ", cplimit))
    fileName_out <- paste0(locOfCPFiles, "ensemble_chill_cutoff_", speciesChoice, "_", k, "_", hem, "_", yearSpan, ".tif")
    r.mean_copy[r.mean_copy < cplimit] <- 0 # not suitable
    r.mean_copy[r.mean_copy >= cplimit] <- 1 # suitable
    print(r.mean_copy)
    titleText <- paste0("ensemble_chill_cutoff_", speciesChoice, "_", k, "_", hem, "_", yearSpan)
    plot(r.mean_copy, main = titleText)
    print(system.time(writeRaster(r.mean_copy, filename = fileName_out, overwrite = TRUE, wopt = woptList)))
    print(paste0("fileName out: ", fileName_out))
    maxVal <- round(max(minmax(r)), 2)
    minVal <- round(min(minmax(r)), 2)
    cat(paste0(red("species: ", speciesChoice, ", ensemble ssp: ", k, ", start year: ", l, ", minVal ", minVal,  ", maxVal ", maxVal, ", fileName out: ", fileName_out), "\n\n"))
  }
  
  f_chillPortionsGraphs <- function(k, l) {
    gc()
    midYear <- l + 9
    yearSpan <- paste0(l, "_", l + yearRange)
    
    # graphics for total chill portions 
    fileNameCP_NH_in <- paste0(locOfCPFiles, "ensemble_chill_portions", "_", k, "_", "NH", "_", yearSpan, ".tif")
    fileNameCP_SH_in <- paste0(locOfCPFiles, "ensemble_chill_portions", "_", k, "_", "SH", "_", yearSpan, ".tif")
    r_CP_NH <- rast(fileNameCP_NH_in)
    r_CP_SH <- rast(fileNameCP_SH_in)
    r_CP <- merge(r_CP_NH, r_CP_SH)
    #    r_CP <- crop(r_CP, extent_noAntarctica)
    r_CP_Rob <- project(r_CP, crsRob)
    r_CP_Rob_df <- as.data.frame(r_CP_Rob, xy = TRUE)
    names(r_CP_Rob_df) <- c("x", "y", "value")
    
    # code for the chill portions values
    titleText <- paste0("Total chill portions, scenario:  ", k,  ", year span: ", gsub("_", "-", yearSpan))
    legendTitle <- "Chill portions"
    #        colorList <- (RColorBrewer::brewer.pal(2, "YlOrRd"))
    #    colorList <- c( "yellow", "green")
    g <- ggplot() +
      geom_tile(data = r_CP_Rob_df, aes(x, y, fill = value), show.legend = TRUE) +
      labs(title = titleText, fill = legendTitle) + theme(plot.title = element_text(size = 12, hjust = 0.5)) +
      labs(x = "", y = "") +
      scale_fill_viridis_c(option = "D",  direction = -1, limits = c(1, 100), na.value = "white") +
      geom_sf(color = "ghostwhite", lwd = 0.2) +
      theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), ) +
      theme(panel.background = element_rect(fill = "aliceblue")) +
      NULL
    print(g)
      fileName_out <- paste0(lofOfGraphicsFiles, "chillPortions/ChillPortions", "_", k, "_", yearSpan, ".png")
    ggsave(filename = fileName_out, plot = g, width = 6, height = 6, units = "in", dpi = 300)
    knitr::plot_crop(fileName_out)
    print(paste0("file name out: ", fileName_out))
    g <- NULL
    
    # graphics for species-specific cutoff
    for (speciesChoice in speciesChoices) {
      cropVals <- get(paste0("majorCropValues", chillLevel))
      cplimit <- cropVals[cropName == speciesChoice, chill_portions]
      cultivar <-  cropVals[cropName == speciesChoice, cultivar]
      fileNameCP_cutoff_NH_in <- paste0(locOfCPFiles, "ensemble_chill_cutoff_", speciesChoice, "_", k, "_", "NH", "_", yearSpan, ".tif")
      fileNameCP_cutoff_SH_in <- paste0(locOfCPFiles, "ensemble_chill_cutoff_", speciesChoice, "_", k, "_", "SH", "_", yearSpan, ".tif")
      
      r_cutoff_NH <- rast(fileNameCP_cutoff_NH_in)
      r_cutoff_SH <- rast(fileNameCP_cutoff_SH_in)
      r_cutoff <- merge(r_cutoff_NH, r_cutoff_SH)
      r_cutoff_Rob <- project(r_cutoff, crsRob)
      r_cutoff_Rob_df <- as.data.frame(r_cutoff_Rob, xy = TRUE)
      names(r_cutoff_Rob_df) <- c("x", "y", "value")
      r_cutoff_Rob_df <- round(r_cutoff_Rob_df, 0) # get the value back to 0s and 1s
      
      # code for the cutoff
      titleText <- paste0("Minimum chill portion requirement met, species: ", speciesChoice, ", \nscenario:  ", k,  ", year span: ", gsub("_", "-", yearSpan))
      legendTitle <- "Adequate chill portions"
      caption <- paste0("Chill portion requirement for ", speciesChoice, " (cultivar ", cultivar, ") is ", cplimit, ".")
      #        colorList <- (RColorBrewer::brewer.pal(2, "YlOrRd"))
      colorList <- c("white", "green")
      #    custom_bins = c(0, 1)
      g <- ggplot() +
        labs(title = titleText, fill = legendTitle) + theme(plot.title = element_text(size = 12, hjust = 0.5)) +
        labs(x = "", y = "", caption = caption) +
        geom_tile(data = r_cutoff_Rob_df, aes(x, y, fill = value), show.legend = FALSE) +
        scale_fill_gradientn(colours=c("white","green")) +
        geom_sf(color = "ghostwhite", lwd = 0.2) +
        theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), ) +
        theme(panel.background = element_rect(fill = "aliceblue"), 
              plot.caption = element_text(hjust = 0, vjust = 7.0, size = 8)
        ) +
        NULL
      print(g)
      fileName_out <- paste0(lofOfGraphicsFiles, "chillPortions/adeqChillPortions_", speciesChoice, "_", k, "_", yearSpan, ".png")
      ggsave(filename = fileName_out, plot = g, width = 6, height = 6, units = "in", dpi = 300)
      knitr::plot_crop(fileName_out)
      print(paste0("file name out: ", fileName_out))
      g <- NULL
    }
  }
}

# chill portions -----
#scenarios
for (chillLevel in chillLevels) {
  for (hem in hemispheres) {
    for (k in sspChoices) {
      #  k = "ssp585"
      for (l in startYearChoices) {
        for (speciesChoice in speciesChoices) {
          print(system.time(f_chillPortions(k, l, speciesChoice, hem, chillLevel)))
        }
      }
    }
    #historical
    k <- "historical"
    l <- 1991
    for (speciesChoice in speciesChoices) {
      print(system.time(f_chillPortions(k, l, speciesChoice, hem, chillLevel)))
    }
  }
}

# graphics -----

# chill portions graphics, scenarios -----
for (k in sspChoices) {
  for (l in startYearChoices) {
    f_chillPortionsGraphs(k, l)
  }
}

# chill portions graphics, historical -----
{
  k = "historical"
  l <- 1991
  f_chillPortionsGraphs(k, l)
}

# chillPortions ppt -----
library(officer)
library(flextable)
library(magrittr)

defaultWidth <- 10
defaultHeight <- 5
defaultLeft <- 0
defaultTop <- 1
# defaultTopSH <- 4

f_chillportionsPpt <- function(fruit) {
  fileNameStart <- paste0("adeqChillPortions_")
  fileName_in <- paste0(lofOfGraphicsFiles, "chillPortions/", fileNameStart, fruit, "_", k, "_", yearSpan, ".png")
  extImg_cp <- external_img(src = fileName_in, width = defaultWidth, height = defaultHeight)
  my_pres <- add_slide(x = my_pres, layout = 'Title Only', master = 'Office Theme')
  my_pres <- ph_with(x = my_pres, value = extImg_cp, location = ph_location(left = defaultLeft, top = defaultTop, width = defaultWidth, height = defaultHeight - 0.5) )
  return(my_pres)
}

colsToDelete <- names(cropVals)[!names(cropVals) %in% c("cropName", "cultivar", "chill_portions")]
CPs <- cropVals[, (colsToDelete) := NULL ]

# presentation intro -----
titleString <- paste0("Adequate Chill Portions by Fruit, Time Period, and Scenario")
contentString <- paste0("Powerpoint produced on ", Sys.Date())

introText1 <- "These slides display locations where chill portions are adequate 90% of the time of various perennial fruits. "
introText2 <- "The table below shows the crops we’re considering and the chill portion requirements used in the following graphs. "

#dataText5 <- "Crop area is based on the SAGE cropping calendar data set, based in the early 2000s. Areas that are not cropped have an NA value and are displayed in white. Areas in gray have chilling less than the minimum requirements. Areas in yellow have chilling hours between the lower and upper range of the requirements."
dataText <- c(dataText1, dataText2, dataText3, dataText4) #, dataText5)

fp_1 <- fp_text(bold = TRUE, color = "pink", font.size = 0)
fp_2 <- fp_text(bold = FALSE, font.size = 12)
fp_3 <- fp_text(italic = TRUE, color = "black", font.size = 14)

blIntro <- block_list(
  fpar(
    ftext(introText1, fp_2),
    ftext(introText2, fp_2)
  ))

my_pres <- read_pptx()
my_pres <- add_slide(x = my_pres, layout = 'Title Slide', master = 'Office Theme')
my_pres <- ph_with(x = my_pres, value = titleString, location = ph_location_type(type = "ctrTitle"))
my_pres <- ph_with(x = my_pres, value = contentString, location = ph_location_type(type = "subTitle"))

my_pres <- add_slide(my_pres, layout = "Title and Content", master = "Office Theme")
my_pres <-  ph_with(x = my_pres, value = "Introduction", location = ph_location_type(type = "title"))
my_pres <- ph_with(my_pres, head(CPs), location = ph_location(left = 2.5, top = 2.5, width = 4, height = 3) )

my_pres <- ph_with(x = my_pres, value = blIntro, location = ph_location_type(type = "body") )

# presentation for loop -----
#browser()
for (fruit in speciesChoices) {
  ensembleTitle <- paste("Adequate Chill Portions for ", fruit)
  my_pres <- add_slide(x = my_pres, layout = 'Section Header', master = 'Office Theme')
  my_pres <- ph_with(x = my_pres, value = ensembleTitle, location = ph_location_type(type = "title"))
  
  # do historical first, then ssps and future periods
  k <- "historical"
  l <- 1991
  yearSpan <- paste0(l, "_", l + yearRange)
  f_chillportionsPpt(fruit)
  
  for (k in sspChoices) {
    for (l in startYearChoices) {
      yearSpan <- paste0(l, "_", l + yearRange)
      f_chillportionsPpt(fruit)
    }
  }
}

my_pres <- add_slide(my_pres, layout = "Title and Content", master = "Office Theme")
my_pres <-  ph_with(x = my_pres, value = "Data Source", location = ph_location_type(type = "title"))
my_pres <- f_addDataSlide()

print(my_pres, target = "presentations/cmip6/perennials/adequateChillPortions.pptx") %>% browseURL()


