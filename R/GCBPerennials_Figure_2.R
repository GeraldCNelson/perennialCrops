library(terra)
library(viridis)
path_data <- "data/perennials/"
speciesChoices <- c("almond_main", "apple_main",  "cherry_main", "olive_main",  "grape_main")
minArea = 1
crsRob <-  "+proj=robin"
globe_ll <- geodata::world(path = "data-raw/gadm") |> crop(ext(-180, 180, -60, 90))
globe <- globe_ll |> project(crsRob)
e <- ext(globe)
xm <- e[1]; xx <- e[2]; ym <- e[3]; yx <- e[4]
m <- c(1, Inf, 1,
       -Inf, 1, NA)
rclmat <- matrix(m, ncol=3, byrow=TRUE)

# test
speciesChoice <- "grape_main"

graphCounter = 1
for (speciesChoice in speciesChoices) {
  speciesName <- gsub("_main", "", speciesChoice) # needed for the harvested area data
  colList <- c("#D55E00", "#56B4E9", "#009E73", "#F0E442") # orange, blue, green, yellow
  
  # suitable locs in each period, note that these files have the whole planet, not by hemisphere
  suitable_early <- rast(paste0(path_data, "nonlimiting_all_", speciesChoice, "_",  "historical", "_", "good", "_", "1991_2010", ".tif"), lyrs = "combinedSuit") 
  suitable_end <-rast( paste0(path_data, "nonlimiting_all_", speciesChoice, "_",                  "ssp585", "_", "good", "_", "2081_2100", ".tif"), lyrs = "combinedSuit") 
  suitable_end_lo <- rast(paste0(path_data, "nonlimiting_all_", paste0(speciesName, "_lo"), "_",  "ssp585", "_", "good", "_", "2081_2100", ".tif"), lyrs = "combinedSuit") 
  
  suitable_endNearly <- suitable_early * suitable_end 
  suitable_end_new <- (suitable_end - suitable_early) # has values of -1, 0, and 1. One is new suitable land end century, -1 is loss of suitable land
  suitable_early_lost <- suitable_end_new
  suitable_early_lost[suitable_early_lost >= 0] <- NA
  suitable_early_lost[suitable_early_lost== -1] <- 1
  
  suitable_end_new[suitable_end_new == -1] <- NA
  suitable_end_new_lo <- (suitable_end_lo - suitable_end) 
  suitable_end_new_lo[suitable_end_new_lo == -1] <- NA
  suitable_early_recovered <- suitable_early * suitable_end_lo
  combined <- c(suitable_early, suitable_end, suitable_end_lo, suitable_endNearly, suitable_end_new, suitable_end_new_lo, suitable_early_recovered)
  combined[combined == 0 ] <- NA
  combined <- project(combined, crsRob)
  names(combined) <- c("suitable_early_lost", "suitable_end", "suitable_end_lo", "suitable_endNearly", "suitable_end_new", "suitable_end_new_lo", "suitable_early_recovered")
  grat <- graticule(30, 30, crs = crsRob) |> crop(e)
  
  
  outf <- paste0("graphics/figure_2_", letters[graphCounter], "_", speciesName, ".pdf")
  pdf(outf)
  plot(grat, col = "gray", background = "azure", lty = 2, mar = c(.1,.1,.1,.1), labels = FALSE)
  polys(globe, col=gray(.99), lwd = .5, alpha = 1)
  plot(combined[[1]], add = TRUE, axes = FALSE, col = colList[[1]], alpha = 1, legend = FALSE) # suitable early
  plot(combined[[4]], add = TRUE, axes = FALSE, col = colList[[2]], legend = FALSE)
  plot(combined[[5]], add = TRUE, axes = FALSE, col = colList[[3]], legend = FALSE)
  plot(combined[[6]], add = TRUE, axes = FALSE, col = colList[[4]], alpha = 1, legend = FALSE) 
  text(xm-xm/6, ym, cex = 1, (bquote(paste((bold(.(letters[graphCounter])))*' ', .(speciesName)))))
  if (speciesName == "grape") {
    ltext = c("Recent historical \nsuitability lost", "Recent historical\nsuitability retained", "End century\nnew suitable areas", "End century suitable\nareas, low chill portions")
    legend(xm - xm/4, ym, fill = colList, legend = ltext, text.width = strwidth(ltext)[4]/1.9, horiz = TRUE,  cex = .7, x.intersp = .0025, lwd = 0, bty = "n", xjust = 0)
  }
  graphCounter <- graphCounter + 1
  dev.off()
  system2('pdfcrop', c(outf, outf)) # gets rid of white space around the figure in the pdf
}

