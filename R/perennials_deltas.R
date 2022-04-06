source("R/perennialFunctions.R")
# source("R/perennialsPrep.R") # creates the data tables majorCropValues_main, majorCropValues_lo and majorCropValues_hi. Included in perennialFunctions.R
woptList <- list(gdal=c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL = 6", "NUM_THREADS=ALL_CPUS"))
terraOptions(memfrac = 0.9, copies = 1, progress = 10, tempdir="data/ISIMIP", verbose = TRUE) # memfrac = 2, 
extent_noAntarctica <- ext(-180, 180, -60, 90) #-60 gets rid of Antarctica for global
crsRob <-  "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
#coastline data
coastline <- vect("data-raw/regionInformation/ne_50m_coastline/ne_50m_coastline.shp")
coastline_cropped <- crop(coastline, extent_noAntarctica )
coastline_cropped_Rob <- project(coastline_cropped, crsRob)
coastline_cropped_Rob_sf <- sf::st_as_sf(coastline_cropped_Rob)
colList <- c("#e66101","#fdae61", "#abd9e9", "#2c7bb6")

# import regions to zoom in on 
library(readxl)
library(PROJ)
library(data.table)
library(ggforce)
library(dplyr)
library(ggplot2)
regions <- as.data.table(read_excel("data-raw/regionInformation/regions.xlsx"))
test_loc <- "west_coast_US"
ext_west_coast_US <- ext(as.numeric(regions[location == test_loc,2:5]))
xmin <- as.numeric(regions[location == test_loc,2])
xmax <- as.numeric(regions[location == test_loc,3])
ymin <- as.numeric(regions[location == test_loc,4])
ymax <- as.numeric(regions[location == test_loc,5])
xymin_rob <- as.numeric(proj_trans(cbind(xmin, ymin), target = RobinsonProj, source = "epsg:4326"))
xymax_rob <- as.numeric(proj_trans(cbind(xmax, ymax), target = RobinsonProj, source = "epsg:4326"))
#ext_rob <- ext(xminmax_rob, yminmax_rob)

yearRange <- 19
locOfClimFiles <- "climdata"

f_convert_graphics <- function(rawIn, desc) {
  rawIn_rob <- project(rawIn, crsRob)
  rawIn_rob_df <- as.data.frame(rawIn_rob, xy = TRUE)
  names(rawIn_rob_df) <- c("x", "y", "value")
  rawIn_rob_df$type = desc
  rawIn_rob_df$value <- round(rawIn_rob_df$value, 0)
  rawIn_rob_df$value[rawIn_rob_df$value == 0] <- NA
  return(rawIn_rob_df)
}

f_getArea <- function(r, layer) {
  r_sub <- subset(r, layer)
  r_sub[r_sub == 0] <- NA
  r_sub_area <- (expanse(r_sub, unit = "km"))/1000
  r_sub_area <- round(r_sub_area, 0)
}

yearSpan_early <- "1991_2010"
yearSpan_mid <- "2041_2060"
yearSpan_end <- "2081_2100"
speciesChoices <- unique(cropVals$cropName)

#test values, perennials -----
speciesChoice <- "cherry_main"
i <- "cherry"
# _lo data
# varietiesChoiceSuffix <- "_main" # whether to use lo, main, or hi cp values. moved to perennialsFunctions.R
speciesNames <- gsub(varietiesChoiceSuffix, "", speciesChoices) 

# graphing ------
# empty table for area data -----
areaDat <- data.table(species = character(), earlyC_cp = character(), endC_main_cp = character(),  endC_lo_cp = character(), earlyC_heatSuit = character(), endC_heatSuit = character(),  endC_main_combinedSuit = character(), endC_lo_combinedSuit = character())
k = "ssp585"
if (k == "ssp585") scenarioName = "SSP5-8.5"
if (k == "ssp126") scenarioName = "SSP1-2.6"

for (i in speciesNames) {
  summerHeat <- majorCropValues_main[cropName == paste0(i, "_main"), summer_heat_threshold]
  chillPortions_main <- majorCropValues_main[cropName == paste0(i, "_main"), chill_portions]
  chillPortions_lo <- majorCropValues_lo[cropName == paste0(i, "_main"), chill_portions]
  r_early_in <- paste0(locOfDataFiles_perennials, "nonlimiting_all_", paste0(i, "_main"), "_",  "historical", "_", suitabilityLevel, "_", yearSpan_early, ".tif")
  r_mid_in <- paste0(locOfDataFiles_perennials, "nonlimiting_all_", paste0(i, "_main"), "_",  k, "_", suitabilityLevel, "_", yearSpan_mid, ".tif")
  r_end_in <- paste0(locOfDataFiles_perennials, "nonlimiting_all_", paste0(i, "_main"), "_",  k, "_", suitabilityLevel, "_", yearSpan_end, ".tif")
  r_end_in_lo <- paste0(locOfDataFiles_perennials, "nonlimiting_all_", paste0(i, "_lo"), "_",  k, "_", suitabilityLevel, "_", yearSpan_end, ".tif")
  
  r_early <- rast(r_early_in)
  r_mid <- rast(r_mid_in)
  r_end <- rast(r_end_in)
  r_end_lo <- rast(r_end_in_lo)
  # suitable locs in each period
  suitable_early <- r_early$combinedSuit
  # suitable_early[suitable_early == 0] <- NA
  suitable_mid <- r_mid$combinedSuit
  suitable_end <- r_end$combinedSuit
  suitable_end_lo <- r_end_lo$combinedSuit
  # common suitability in early and mid period
  suitable_mid_early <- suitable_early * suitable_mid # locations where both period values are 1 remain 1
  # common suitability in early and end period; will assume it is also common with mid
  suitable_end_early <- suitable_early * suitable_end # locations where both period values are 1 remain 1
  early_end_loss <- suitable_end_early - suitable_early
  suitable_end_early <- suitable_end_early * 2
#  test <- suitable_early + suitable_mid + suitable_end
  suitable_end_new <- (suitable_end - suitable_early) * 3
  suitable_end_new <- app(suitable_end_new, fun=function(x){ x[x < 3] <- 0; return(x)} )
  suitable_loss_early <- suitable_mid - suitable_early # 0 = both period values are 1 or both are 0; 1 = mid suitable, early not, -1 = suitable early, not in mid
  lo_benefits <- (suitable_end_lo - suitable_end) * 4 # new areas where a low chill portion requirement would become suitable
#  test <- c(suitable_early, suitable_end_early)
  
  # get areas ------
  earlyCent_main_cp <- f_getArea(r_early, "chillPortionsSuit")
  earlyCent_main_heatSuit <- f_getArea(r_early, "heatSuit")
  
  endCent_main_cp <- f_getArea(r_end, "chillPortionsSuit")
  endCent_lo_cp <- f_getArea(r_end_lo, "chillPortionsSuit")
  endCent_main_heatSuit <- f_getArea(r_end, "heatSuit")
  endCent_main_combinedSuit <- f_getArea(r_end, "combinedSuit")
  endCent_lo_combinedSuit <- f_getArea(r_end_lo, "combinedSuit")
  area_combined <- list(i, earlyCent_main_cp, endCent_main_cp, endCent_lo_cp, earlyCent_main_heatSuit, endCent_main_heatSuit, endCent_main_combinedSuit, endCent_lo_combinedSuit)
  areaDat <- rbind(areaDat, area_combined)
  
  # graphics, convert to dfs -----
   suitable_early_rob_df <- f_convert_graphics(suitable_early, "suitable early") 
   suitable_end_rob_df <- f_convert_graphics(suitable_end, "suitable end") 
   suitable_end_early_rob_df <- f_convert_graphics(suitable_end_early, "suitable, early & end century") 
   suitable_end_new_rob_df <- f_convert_graphics(suitable_end_new, "new suitable, end century")
   lo_benefits_rob_df <- f_convert_graphics(lo_benefits, "added end-century low-chill suitability") 
 
  r_combined_df <- rbind(suitable_early_rob_df, suitable_end_early_rob_df, suitable_end_new_rob_df, lo_benefits_rob_df)
  #  r_combined_df$type <- factor(r_combined_df$type, levels = c("suitable early", "suitable, early & end century", "added end-century low-chill suitability"))
  r_combined_df$value <- round(r_combined_df$value, 0)
#  r_combined_df[r_combined_df == 0] <- NA
 r_combined_df$value = as.factor(r_combined_df$value)
  #forcats::fct_rev(r_combined_df$value)
  # now get harvested area map
 if (i == "winegrape") i <- "grape"
  rInArea <- rast(paste0(locOfHarvestDataFiles, i,"/",  i, "_HarvestedAreaHectares.tif"))
  harvestArea_earlyCent <- aggregate(rInArea, fact = 6, fun = "sum") # convert 5 arc minutes to 1/2 degrees
  maskMin <- switch(
    i,
    "almond" = 100,
    "apple" = 100,
    "cherry" = 100,
    "grape" = 100,
    "olive" = 100
  )
  if (i == "grape") i <- "winegrape"
  harvestArea_earlyCent[harvestArea_earlyCent < maskMin] <- 0 # set minimum area to be greater than 0 hectares per grid cell
  harvestArea_earlyCent[harvestArea_earlyCent > 0] <- 1
  harvestArea_earlyCent <- project(harvestArea_earlyCent, crsRob)
  harvestArea_df <- as.data.frame(harvestArea_earlyCent, xy = TRUE)
  names(harvestArea_df) <-   c("x", "y", "value_harvest")
  harvestArea_df <- round(harvestArea_df, 0)
  harvestArea_df[harvestArea_df == 0] <- NA
 # browser()
  fileName_out <- paste0(lofOfGraphicsFiles, "perennials/delta_suitability_", i, ".png")
  
  titleText <- paste0("   Changing suitability for ", i, "\nfrom early to end 21st century, climate scenario ", scenarioName)
  legendTitle <- "Suitability"
  cropNameVal <- paste0(i, varietiesChoiceSuffix)
  CPfruit <- majorCropValues_main[cropName == cropNameVal, chill_portions]
  CPfruit_lo <- majorCropValues_lo[cropName == paste0(i, "_lo"), chill_portions]
  summerHeat <- majorCropValues_main[cropName == cropNameVal, summer_heat_threshold]
  cultivar <-  majorCropValues_main[cropName == cropNameVal, cultivar]
  gddsFruit <- majorCropValues_main[cropName == cropNameVal, gdd]
  captionString <- "Note: Suitable locations are where the %s cultivar, %s, is not limited by temperature, has at least %s chill portions, fewer than %s days of spring frost risk, \na minimum of %s growing degree days and fewer than %s days of summer heat greater than %s°C. Locations in red are where the early century suitability is lost by the end of the century. \nGreen areas become suitable by end century. Blue areas become suitable with a variety with a reduced chill requirement of %s. Light gray shading indicates early 21st century \narea for all %s varieties according to data from http://www.earthstat.org."
  caption <- sprintf(paste(captionString, collapse = " ") , i, cultivar, CPfruit, frostRiskDays[2], gddsFruit, heatRiskDays[2], summerHeat, CPfruit_lo, i)
  g <- ggplot() +
    geom_tile(data = dplyr::filter(r_combined_df,  !is.na(value)), aes(x, y, fill = value), show.legend = TRUE, stat = "identity", position = "identity") +
    scale_fill_manual(labels = c("Suitability lost", "Suitability retained", "New suitable areas", "New areas, end century\n with low chill portions"), values = colList, na.value = 'white') +
    
#    scale_fill_manual(labels = c("Suitability lost", "Suitability retained", "New suitable areas", "New areas, end century\n with low chill portions"), values = c("red", "yellow", "mediumseagreen", "blue"), na.value = 'white') +
    labs(title = titleText, fill = legendTitle, x = "", y = "") + #, caption = caption
    geom_sf(data = coastline_cropped_Rob_sf,  color="black", size = 0.1) +
    theme_bw() +
    theme(
      legend.text.align = 0,
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(size = 12, hjust = 0.5),
      plot.caption = element_text(hjust = 0, vjust = 7.0, size = 7)
    ) +
    #      geom_tile(data = harvestArea_df, aes(x, y), fill = "gray", alpha = .2, show.legend = FALSE)
    geom_tile(data = dplyr::filter(harvestArea_df, !is.na(value_harvest)), 
              aes(x = x, y = y), fill = "grey24", alpha = .2, show.legend = FALSE) # +
  # geom_tile(data = dplyr::filter(suitableArea_historical_df, !is.na(value)), 
  #           aes(x = x, y = y), fill = "chocolate1", alpha = .2, show.legend = FALSE) +
  #  facet_wrap(~period, nrow = 3) +
  scale_x_continuous(sec.axis = sec_axis(~ . , name = i, breaks = NULL, labels = NULL))
  
  print(g)
  ggsave(filename = fileName_out, plot = g, width = 8, height = 8, units = "in", dpi = 300)
  knitr::plot_crop(fileName_out) # gets rid of margins around the plot
  print(paste0("file name out: ", fileName_out))
  g <- NULL
}

# code to show different suitabilities
suit_heat_lo_end <- r_end_lo$heatSuit
suit_heat_main_early <- r_early$heatSuit
suit_chill_lo_end <- r_end_lo$chillPortionsSuit
suit_chill_main_end <- r_end$chillPortionsSuit
suit_all_main_end <- r_end$combinedSuit
suit_all_lo_end <- r_end_lo$combinedSuit
suit_combined_endCent_heat <- r_end_lo$heatSuit + r_end_lo$chillPortionsSuit

plot(suit_heat_lo_end, main = paste0("crop ", i , ", Summer heat suitability, ", yearSpan_end, ", scenario ", k))
plot(suit_heat_main_early, main = paste0("crop ", i , ", Summer heat suitability, early ", yearSpan_early))
plot(suit_chill_lo_end, main = paste0("crop ", i , ", chill portions suitability, (low value), \nperiod ", yearSpan_end, ", scenario ", k))
plot(suit_chill_main_end, main = paste0("crop ", i , ", chill portions suitability, (main value), \nperiod ", yearSpan_end, ", scenario ", k))
plot(suit_combined_endCent_heat, main = paste0("Chill portions, heat stress and combined suitability,  \nperiod ", yearSpan_end, ", scenario ", k))
plot(suit_all_main_end, main = paste0("Combined suitability, all, \nperiod ", yearSpan_end, ", scenario ", k))
plot(suit_all_lo_end, main = paste0("Combined suitability, all, lo chill, \nperiod ", yearSpan_end, ", scenario ", k))

#suit_heat_lo_end_rob_df <- f_convert_graphics(suit_heat_lo_end, "suitable_heat_low_end") 
suit_chill_lo_end_rob_df <- f_convert_graphics(suit_chill_lo_end, "suitable_chill_low_end") 
suit_heat_main_end_rob_df <- f_convert_graphics(suit_heat_lo_end, "suitable_heat_main_end") 
suit_chill_main_end_rob_df <- f_convert_graphics(suit_chill_main_end, "suitable_chill_main_end") 

# combine different suitables -----
r_suit_combined_df <- rbind(suit_heat_main_end_rob_df, suit_chill_main_end_rob_df, suit_chill_lo_end_rob_df)
r_suit_combined_df$value <- round(r_suit_combined_df$value, 0)
#  r_combined_df[r_combined_df == 0] <- NA
r_suit_combined_df$value = as.factor(r_suit_combined_df$value)
titleText <- paste0("Suitability typesfor ", i, " end 21st century, climate scenario ", scenarioName)

g <- ggplot() +
  geom_tile(data = dplyr::filter (r_suit_combined_df, !is.na(value)), aes(x, y, fill = value), show.legend = TRUE, stat = "identity", position = "identity") +
  # scale_fill_manual(labels = c("Suitability lost", "Suitability retained", "New suitable areas", "New areas, end century\n with low chill portions"), values = c("red", "yellow", "mediumseagreen", "blue"), na.value = 'white') +
  scale_fill_manual(labels = c("Suitability lost", "Suitability retained", "New suitable areas", "New areas, end century\n with low chill portions"), values = colList, na.value = 'white') +
  
  labs(title = titleText, fill = legendTitle, x = "", y = "", caption = caption) + 
  #   labs( x = "", y = "") +
  geom_sf(data = coastline_cropped_Rob_sf,  color="black", size = 0.1) +
  theme_bw() +
  theme(
    legend.text.align = 0,
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 12, hjust = 0.5),
    plot.caption = element_text(hjust = 0, vjust = 7.0, size = 7)
  ) +
  #      geom_tile(data = harvestArea_df, aes(x, y), fill = "gray", alpha = .2, show.legend = FALSE)
  geom_tile(data = dplyr::filter(harvestArea_df, !is.na(value_harvest)), 
            aes(x = x, y = y), fill = "grey24", alpha = .2, show.legend = FALSE) # +
# geom_tile(data = dplyr::filter(suitableArea_historical_df, !is.na(value)), 
#           aes(x = x, y = y), fill = "chocolate1", alpha = .2, show.legend = FALSE) +
#  facet_wrap(~period, nrow = 3) +
scale_x_continuous(sec.axis = sec_axis(~ . , name = i, breaks = NULL, labels = NULL))
NULL

print(g)
ggsave(filename = fileName_out, plot = g, width = 8, height = 8, units = "in", dpi = 300)
knitr::plot_crop(fileName_out) # gets rid of margins around the plot
print(paste0("file name out: ", fileName_out))

# create table of area suitabilities -----
library(flextable)
library(officer)
areaTable <- areaDat
setcolorder(areaTable, neworder = c("species", "earlyC_cp", "earlyC_heatSuit", "endC_main_cp", "endC_lo_cp",  "endC_heatSuit", "endC_main_combinedSuit", "endC_lo_combinedSuit"))
areaTable[, species := stringr::str_to_title(species)][species == "Winegrape", species := "Wine grape"]
areaTable_flex <- flextable(areaTable)
typology_all <- data.frame(
  col_keys_all = names(areaTable),
  what = c("Species", #1
           "Early Century Suitability (1000 sq. km)", "Early Century Suitability (1000 sq. km)", "End Century Suitability (1000 sq. km)", #3
           "End Century Suitability (1000 sq. km)", "End Century Suitability (1000 sq. km)", "End Century Suitability (1000 sq. km)", "End Century Suitability (1000 sq. km)"), #4
  measure = c("Species", #1
              "Chill portions", "Summer heat", "Chill portions", "Chill portions, low variety", "Summer heat", "Combined suitability", "Combined suitability, low chill portion variety"), # 8
  stringsAsFactors = FALSE )

typology_footer <- data.frame(col_keys = "Species",
                              notes = "some text", stringsAsFactors = FALSE) 
areaTable_flex <- set_header_df(areaTable_flex, mapping = typology_all, key = "col_keys_all")
areaTable_flex <- align(areaTable_flex, align = "right", i = 1:2, j = 1:4, part = "all")
areaTable_flex <- merge_h(areaTable_flex, part = "header")
areaTable_flex <- align(areaTable_flex, align = "justify", i = 1, j = 1, part = "header")
areaTable_flex <- merge_v(areaTable_flex, j = c("species"), part = "header")
areaTable_flex <- theme_vanilla(areaTable_flex)
areaTable_flex <- fix_border_issues(areaTable_flex)
areaTable_flex <- autofit(areaTable_flex)
areaTable_flex <- colformat_num(areaTable_flex,i = 2:7, j = 3:8, big.mark=",", decimal.mark = ".", na_str = "N/A")
areaTable_flex <- width(areaTable_flex, j = 8, width = 1.5)
# areaTable_flex <- align(areaTable_flex, align = "center", part = "header")
#areaTable_flex <- align(areaTable_flex, align = "right", part = "body")
areaTable_flex <- set_footer_df(areaTable_flex, mapping = typology_footer, key = "col_keys")
# areaTable_flex <- add_footer(areaTable_flex, Species =  "Notes to be added.")
# areaTable_flex <- merge_at(areaTable_flex, part = "footer")
areaTable_flex <- height(areaTable_flex, height = .25, part = "body")
areaTable_flex

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
save_as_docx(areaTable_flex, values = NULL, path = "results/perennials_suitability.docx", pr_section = prsect_area)

# inset map ------
# not currently used
regions_latlon <- as.data.table(read_excel("data-raw/regionInformation/regions.xlsx"))
test_loc <- "west_coast_US"
# get the xy coordinates for the area you want to pull out
ext_west_coast_US <- ext(as.numeric(regions_latlon[location == test_loc,2:5]))
xmin <- as.numeric(regions[location == test_loc,2])
xmax <- as.numeric(regions[location == test_loc,3])
ymin <- as.numeric(regions[location == test_loc,4])
ymax <- as.numeric(regions[location == test_loc,5])
# convert to Robinson projection
xymin_rob <- as.numeric(proj_trans(cbind(xmin, ymin), target = RobinsonProj, source = "epsg:4326"))
xymax_rob <- as.numeric(proj_trans(cbind(xmax, ymax), target = RobinsonProj, source = "epsg:4326"))
xmin_rob <- xymin_rob[1]
ymin_rob <- xymin_rob[2]
xmax_rob <- xymax_rob[1]
ymax_rob <- xymax_rob[2]

# add rectngle around location to clip
boxCoords <- data.frame(lon = c(xmin_rob, xmax_rob, xmax_rob, xmin_rob, xmin_rob), lat = c(ymin_rob, ymin_rob, ymax_rob, ymax_rob, ymin_rob))
g_test <- g + geom_path(data = boxCoords, aes(x = lon, y = lat))
g_small <- g + 
  scale_x_continuous(limits = c(xmin_rob, xmax_rob)) +
  scale_y_continuous(limits = c(ymin_rob, ymax_rob))+
  guides(fill=FALSE) + #remove legend
  labs(title = "") # remove title
  
grobin <- ggplotGrob(g_small)

# where to put the grob -----
# South Pacific
ll_lower_left <- c(-180, -80)
ll_lower_left_rob <- as.numeric(proj_trans(cbind(ll_lower_left[1], ll_lower_left[2]), target = RobinsonProj, source = "epsg:4326"))

mag <- 8
xmin_gr = ll_lower_left_rob[1] # - 5000000
xmax_gr <- xmin_gr + (xmax_rob - xmin_rob) * mag
ymin_gr = ll_lower_left_rob[2]
ymax_gr <- ymin_gr + (ymax_rob - ymin_rob) * mag

g + annotation_custom(grobin, xmin = xmin_gr, xmax =  xmax_gr, 
                      ymin = ymin_gr, ymax = ymax_gr) +
  geom_path(data = )

# connect the dots -----
grobPoints <- c(xmin_gr, xmax_gr, ymin_gr, ymax_gr)
mapPoints <- c(xmin_rob, xmax_rob, ymin_rob, ymax_rob)

grobPointsX <- c(xmin_gr, xmin_rob)
grobPointsY <- c(ymin_gr, ymin_rob)

pointsToConnect_ll <- as.data.frame(cbind(grobPointsX, grobPointsY))

g + annotation_custom(grobin, xmin = xmin_gr, xmax =  xmax_gr, 
                      ymin = ymin_gr, ymax = ymax_gr) +
  geom_path(data = pointsToConnect_ll, aes(x = grobPointsX,
                                      y = grobPointsY))
# source: https://geocompr.github.io/post/2019/ggplot2-inset-maps/
library(cowplot)
library(rcartocolor)
pointsToConnect_ll <- as.data.frame(c(0))
gg_inset_map1 = ggdraw() +
  draw_plot(g) +
  draw_plot(grobin, x = 0.0, y = 0.32, width = 0.2, height = 0.2) +
  draw_line(x = c(0,0.2), y = c(0.1, 0.3))
                                           

#source: https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
scales::show_col(colorBlindGrey8)

# alternate source
#https://colorbrewer2.org/#type=diverging&scheme=RdYlBu&n=4

