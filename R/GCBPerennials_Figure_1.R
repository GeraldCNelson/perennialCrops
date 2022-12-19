library(terra)
library(viridis)
path_data <- "data/perennials/"
speciesChoices <- c("almond_main", "apple_main",  "cherry_main",  "grape_main", "olive_main")
minArea = 1
crsRob <-  "+proj=robin"
globe <- geodata::world(path = "data-raw/gadm") |> crop(ext(-180, 180, -60, 90)) |> project(crsRob)
e <- ext(globe)
xm <- e[1]; xx <- e[2]; ym <- e[3]; yx <- e[4]
f_harvestArea <- function(speciesName, minArea) {
  r_area <- rast(paste0("data-raw/crops/HarvestedAreaYield175Crops_Geotiff/GeoTiff/", speciesName,"/",  speciesName, "_HarvestedAreaHectares.tif"))
  harvestArea <- aggregate(r_area, fact = 6, fun = "sum") # convert 5 arc minutes to 1/2 degrees
  maskMin <- switch(
    speciesName,
    "almond" = minArea,
    "apple" = minArea,
    "cherry" = minArea,
    "grape" = minArea,
    "olive" = minArea
  )
  harvestArea[harvestArea < maskMin] <- 0 # set minimum area to be greater than minarea hectares per grid cell
  harvestArea[harvestArea > 1] <- 1
  harvestArea[harvestArea < 1] <- NA
  harvestArea <- project(harvestArea, crsRob) |> crop(e)
  
  return(harvestArea)
}

# test
speciesChoice <- "olive_main"
graphCounter = 1
for (speciesChoice in speciesChoices) {
  speciesName <- gsub("_main", "", speciesChoice) # needed for the harvested area data
  harvestArea_earlyCent_rob <- f_harvestArea(speciesName, minArea = 1) 
  
  inF_hist_suit <- paste0("data/perennials/nonlimiting_all_", speciesChoice, "_", "historical", "_", "good", "_", "1991_2010", ".tif")
  suitableArea_historical_rob <- rast(inF_hist_suit, lyrs = 1)  |> project(crsRob) 
  suitableArea_historical_rob[suitableArea_historical_rob > 1] <- 1
  suitableArea_historical_rob[suitableArea_historical_rob < 1] <- NA 
  inF_end_suit_ssp585 <- paste0("data/perennials/nonlimiting_all_", speciesChoice, "_", "ssp585", "_", "good", "_", "2081_2100", ".tif")
  suitableArea_end_ssp585_rob <- rast(inF_end_suit_ssp585, lyrs = 1)  |> project(crsRob) 
  suitableArea_end_ssp585_rob[suitableArea_end_ssp585_rob > 1] <- 1
  suitableArea_end_ssp585_rob[suitableArea_end_ssp585_rob < 1] <- NA
  grat <- graticule(30, 30, crs = crsRob) |> crop(e)
  outf <- paste0("graphics/figure_1_", letters[graphCounter], "_", speciesName, ".pdf")
  
  pdf(outf)
  layout(matrix(c(1,2), 1, 2, byrow = TRUE))
  plot(grat, col = "gray", background = "azure", lty = 2,  lwd = .5, mar = c(.1,.1,.1,.1), labels = FALSE)
  polys(globe, col=gray(.99), lwd = .1, alpha = 1)
  plot(harvestArea_earlyCent_rob, add = TRUE, axes = FALSE, col = "gray", legend = FALSE)
  plot(suitableArea_historical_rob, add = TRUE, axes = FALSE, col = "green3", legend = FALSE, alpha = .3)
  text(xm-xm/6, ym, cex = .8, (bquote(paste((bold(.(letters[graphCounter])))*' ', .(speciesName)))))
  if (speciesName == "olive") {
    ltext = c("Suitable (Recent historical\n1991-2010)", "Early century\narea")
    legend(xm/2, ym, fill = c("green3", "gray"), legend = ltext, text.width = strwidth(ltext)[4]/1.9, horiz = TRUE,  cex = .5, x.intersp = .0025, lwd = 0, bty = "n", xjust = 0)
  }
  
  graphCounter <- graphCounter + 1
  # future 
  plot(grat, col = "gray", background = "azure", lty = 2,  lwd = .5, mar = c(.1,.1,.1,.1), labels = FALSE)
  polys(globe, col=gray(.99), lwd = .1, alpha = 1)
  plot(harvestArea_earlyCent_rob, add = TRUE, axes = FALSE, col = "gray", legend = FALSE)
  plot(suitableArea_end_ssp585_rob, add = TRUE, axes = FALSE, col = "green3", legend = FALSE, alpha = .3)
  text(xm-xm/8, ym, cex = .8, (bquote(paste((bold(.(letters[graphCounter])))))))
  if (speciesName == "olive") {
    ltext = c("Suitable, end century\n(2081-2100), SSP5-8.5", "Early century\narea")
    legend(xm/2, ym, fill = c("green3", "gray"), legend = ltext, text.width = strwidth(ltext)[4]/1.9, horiz = TRUE,  cex = .5, x.intersp = .0025, lwd = 0, bty = "n", xjust = 0)
  }
  dev.off()
  graphCounter <- graphCounter + 1
  system2('pdfcrop', c(outf, outf)) # gets rid of white space around the figure in the pdf
}
