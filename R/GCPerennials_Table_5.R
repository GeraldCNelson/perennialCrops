# area calculations for perennials paper with Tables 5 and 6 -----
chillLevel <- "_main" # kludge
path_data <- "data/perennials/"
source("R/perennialFunctions.R")
{  
  suitabilityLevels <- "good"
  cropVals <- get(paste0("majorCropValues", chillLevel)) # just need one to get to the crop names
  speciesChoices <- unique(cropVals$cropName)
  speciesNames <- gsub(chillLevel, "", speciesChoices) 
  
  for (chillLevel in c("_lo", "_main")) {
    namelist <- c("species", "cultivar", "chillPortions", "hemisphere", "quality", "ssp", "yearSpan",  "area_suitable", "rasterName")
    dt_area <- data.table(1)[,`:=`((namelist),NA)][,V1:=NULL][.0]
    cropVals <- get(paste0("majorCropValues", chillLevel))
    #scenarios 
    for (speciesName in speciesNames) {
      speciesChoice <- paste0(speciesName, chillLevel)
      CPfruit <- cropVals[cropName == speciesChoice, chill_portions]
      cultivar <-  cropVals[cropName == speciesChoice, cultivar]
      print(paste0("cultivar: ", cultivar))
      for (hem in hemispheres) {
        for (suitabilityLevel in suitabilityLevels) {
          for (k in sspChoices) {
            for (l in startYearChoices) {
              yearSpan <- paste0(l, "_", l + yearRange)
              fileName_in <- paste0(path_data, "nonlimiting_all_", speciesChoice, "_", k, "_", suitabilityLevel, "_", hem, "_", yearSpan, ".tif")
              rastName <- paste0("suitable_", speciesChoice, "_", hem, "_", k, "_", suitabilityLevel, "_", yearSpan)
              # r has all 5 suitability layers, get just the combined one
              r_combSuit <- rast(fileName_in, lyrs = 1 )
              r_combSuit[r_combSuit == 0] <- NA # only want to get area of cells that have a value of 1; ie, suitable
              r_area <- f_getArea(r_combSuit, layer = 1)
              print(paste0("speciesChoice: ", speciesChoice, ", ssp: ", k, ",start year: ", l, ", hemisphere: ", hem, ", chillVal: ", chillLevel, ", r_area: ", r_area))
              dt_area <- rbind(dt_area, list(speciesChoice, cultivar, CPfruit, hem, suitabilityLevel, k, yearSpan,  r_area, rastName))
            }
          }
          #historical
          k <- "historical"
          l <- 1991
          yearSpan <- paste0(l, "_", l + yearRange)
          fileName_in <- paste0(path_data, "nonlimiting_all_", speciesChoice, "_", k, "_", suitabilityLevel, "_", hem, "_", yearSpan, ".tif")
          rastName <- paste0("suitable_", speciesChoice, "_", hem, "_", k, "_", suitabilityLevel, "_", yearSpan)
          # r has all 5 suitability layers, get just the combined one
          r_combSuit <- rast(fileName_in, lyrs = 1 )
          r_combSuit[r_combSuit == 0] <- NA # only want to get area of cells that have a value of 1; ie, suitable
          r_area <- f_getArea(r_combSuit, layer = 1)
          print(paste0("speciesChoice: ", speciesChoice, ", ssp: ", k, ",start year: ", l, ", hemisphere: ", hem, ", chillVal: ", chillLevel, ", r_area: ", r_area))
          dt_area <- rbind(dt_area, list(speciesChoice, cultivar, CPfruit, hem, suitabilityLevel, k, yearSpan,r_area, rastName))
        }
      }
    }
    fileName_out <- paste0(path_data, "areaCalcs", chillLevel, ".csv")
    write.csv(dt_area, file = fileName_out, row.names = FALSE)
    print(paste0("fileName out: ", fileName_out)) # input into table 5
  }
}

# code below comes from https://stackoverflow.com/questions/37376398/how-to-create-an-empty-datatable-with-columns-names-and-then-append-datatables-t
{
  for (chillLevel in c("_lo", "_main"))
    namelist <- c("species", "cultivar", "chillPortions", "hemisphere", "quality", "ssp", "yearSpan", "area_base", "area_ssp", "area_both", "area_hist_loss", " area_ssp_gain", "area_ssp_gain_lo")
  dt_area_delta <- data.table(1)[,`:=`((namelist),NA)][,V1:=NULL][.0]
  cropVals <- get(paste0("majorCropValues", chillLevel))
  
  for (speciesName in speciesNames) {
    speciesChoice <- paste0(speciesName, chillLevel)
    CPfruit <- cropVals[cropName == speciesChoice, chill_portions]
    cultivar <-  cropVals[cropName == speciesChoice, cultivar]
    print(paste0("cropname: ", speciesName, ", cultivar: ", cultivar))
    
    for (hem in hemispheres) {
      suitabilityLevels <- "good"
      for (suitabilityLevel in suitabilityLevels) {
        filename_r_historical_in <- paste0(path_data, "nonlimiting_all_", speciesChoice, "_", "historical", "_", suitabilityLevel, "_", hem, "_", "1991_2010.tif")
        r_historical_combinedSuit <- rast(filename_r_historical_in, lyrs = "combinedSuit") # combined suitability means all of the temp metrics are satisfied (+1)
        r_base <- r_historical_combinedSuit; r_base[r_base < 1] <- NA
        area_base <- f_getArea(r_base, layer = 1)
        # do scenarios and deltas
        # use just end of century and ssp585
        startYearChoices <- 2081; sspChoices <- "ssp585"
        for (l in startYearChoices) {
          yearSpan <- paste0(l, "_", l + yearRange)
          for (k in sspChoices) {
            filename_r_ssp_in <- paste0(path_data, "nonlimiting_all_", speciesChoice, "_",  k, "_", suitabilityLevel, "_", hem, "_", yearSpan, ".tif")
            r_ssp_combinedSuit <- rast(filename_r_ssp_in, lyrs = "combinedSuit")
            r_sum <- r_ssp_combinedSuit + r_historical_combinedSuit # 2 - in both periods, 1 in 1 of the periods, 0 not in either
            r_both <- r_sum; r_both[!r_both == 2] <- NA # get area only for pixels with value of 2
            r_both <- r_both/2 # convert to 1/NA
            r_ssp <- r_ssp_combinedSuit; r_ssp[!r_ssp == 1] <- NA
            r_delta <- r_ssp_combinedSuit - r_historical_combinedSuit # 1 in 1 of the periods, not the other, 0 not in either or in both, -1 in historical but not in the future
            r_hist_loss <- r_delta; r_hist_loss[!r_hist_loss == -1] <- NA # get pixels only with the value of -1, lost historical area
            r_hist_loss <- r_hist_loss * -1 # convert the -1s to 1s
            r_ssp_gain <- r_delta; r_ssp_gain[!r_ssp_gain == 1 ] <- NA # new pixels from SSP; remove historical pixels lost (-1) and pixels where either both are 1 or both are 0
            
            fileName_r_end_lo_in <- paste0(path_data, "nonlimiting_all_", paste0(speciesName, "_lo"), "_",  k, "_", suitabilityLevel, "_", hem, "_", yearSpan, ".tif")
            r_ssp_combinedSuit_lo <- rast(fileName_r_end_lo_in, lyrs = "combinedSuit")
            r_ssp_lo <- r_ssp_combinedSuit_lo; r_ssp_lo[!r_ssp_lo == 1 ] <- NA # new pixels from SSP; remove historical pixels lost (-1) and pixels where either both are 1 or both are 0
            r_ssp_delta_lo <- r_ssp_combinedSuit_lo - r_ssp_combinedSuit 
            r_ssp_gain_lo <- r_ssp_delta_lo; r_ssp_gain_lo[!r_ssp_gain_lo == 1 ] <- NA # new pixels from SSP; remove historical pixels lost (-1) and pixels where either both are 1 or both are 0
            
            #area calcs, scenarios
            area_both <- f_getArea(r_both, layer = 1)
            area_hist_loss <- f_getArea(r_hist_loss, layer = 1)
            area_ssp <- f_getArea(r_ssp, layer = 1)
            area_ssp_gain <- f_getArea(r_ssp_gain, layer = 1)
            area_ssp_lo <- f_getArea(r_ssp_lo, layer = 1)
            area_ssp_gain_lo <- f_getArea(r_ssp_gain_lo, layer = 1)
            dt_area_delta <- rbind(dt_area_delta, list(speciesChoice, cultivar, CPfruit, hem, suitabilityLevel, k, yearSpan, area_base, area_ssp, area_both, area_hist_loss, area_ssp_gain, area_ssp_gain_lo))
          }
          test <- c(r_base, r_ssp, r_hist_loss, r_ssp_gain, r_ssp_lo, r_ssp_gain_lo)
          names(test) <- c("r_base", "r_ssp", "r_hist_loss", "r_ssp_gain", "r_ssp_lo", "r_ssp_gain_lo")
          fileName_various_out <- paste0(path_data, speciesName, "_", "areaSuitabilityCombos", "_", hem, ".tif")
          writeRaster(test, filename = fileName_various_out, overwrite = TRUE, wopt = woptList)
          print(paste0("fileName_various_out: ", fileName_various_out))
        }
      }
    }
  }
  
  fileName_out <- paste0(path_data, "areaCalcs_delta", chillLevel, ".csv")
  write.csv(dt_area_delta, file = fileName_out, row.names = FALSE)
  print(paste0("fileName out: ", fileName_out))
}

# harvest and suitability area calcs ------
# area common to mid and end century -----
for (chillLevel in c("_lo", "_main")) {
  namelist <- c("species", "cultivar", "chillPortions", "hemisphere",  "ssp", "yearSpan", "area_common")
  dt_area_common_mid <- data.table(1)[,`:=`((namelist),NA)][,V1 := NULL][.0]
  dt_area_common_mid[, area_common := as.numeric(area_common)]
  dt_area_common_end <-  dt_area_common_mid
  for (speciesName in speciesNames) {
    cropVals <- get(paste0("majorCropValues", chillLevel))
    speciesChoice <- paste0(speciesName, chillLevel)
    
    CPfruit <- cropVals[cropName == speciesChoice, chill_portions]
    cultivar <-  cropVals[cropName == speciesChoice, cultivar]
    
  #suitable areas -----
    # suitable area end century
    fileName_start <- paste0(path_data, "nonlimiting_all_", speciesChoice, "_")
    
    # suitable historical area
    fileName_in_historical_SH <- paste0(fileName_start, "historical", "_", suitabilityLevel, "_", "SH", "_", "1991_2010", ".tif")
    fileName_in_historical_NH <- paste0(fileName_start, "historical", "_", suitabilityLevel, "_", "NH", "_", "1991_2010", ".tif")
    suitableArea_historical_SH <- rast(fileName_in_historical_SH, lyrs = "combinedSuit")
    suitableArea_historical_NH <- rast(fileName_in_historical_NH, lyrs = "combinedSuit")
    
    fileName_in_ssp126_SH_end <- paste0(fileName_start, "ssp126", "_", suitabilityLevel, "_", "SH", "_", "2081_2100", ".tif")
    fileName_in_ssp126_NH_end <- paste0(fileName_start, "ssp126", "_", suitabilityLevel, "_", "NH", "_", "2081_2100", ".tif")
    fileName_in_ssp585_SH_end <- paste0(fileName_start, "ssp585", "_", suitabilityLevel, "_", "SH", "_", "2081_2100", ".tif")
    fileName_in_ssp585_NH_end <- paste0(fileName_start, "ssp585", "_", suitabilityLevel, "_", "NH", "_", "2081_2100", ".tif")
    
    suitableArea_ssp126_SH_end <- rast(fileName_in_ssp126_SH_end, lyrs = "combinedSuit")
    suitableArea_ssp126_NH_end <- rast(fileName_in_ssp126_NH_end, lyrs = "combinedSuit")
    suitableArea_ssp585_SH_end <- rast(fileName_in_ssp585_SH_end, lyrs = "combinedSuit")
    suitableArea_ssp585_NH_end <- rast(fileName_in_ssp585_NH_end, lyrs = "combinedSuit")
    
    # suitable area mid century
    fileName_in_ssp126_SH_mid <- paste0(fileName_start, "ssp126", "_", suitabilityLevel, "_", "SH", "_", "2041_2060", ".tif")
    fileName_in_ssp126_NH_mid <- paste0(fileName_start, "ssp126", "_", suitabilityLevel, "_", "NH", "_", "2041_2060", ".tif")
    fileName_in_ssp585_SH_mid <- paste0(fileName_start, "ssp585", "_", suitabilityLevel, "_", "SH", "_", "2041_2060", ".tif")
    fileName_in_ssp585_NH_mid <- paste0(fileName_start, "ssp585", "_", suitabilityLevel, "_", "NH", "_", "2041_2060", ".tif")
    
    suitableArea_ssp126_SH_mid <- rast(fileName_in_ssp126_SH_mid, lyrs = "combinedSuit")
    suitableArea_ssp126_NH_mid <- rast(fileName_in_ssp126_NH_mid, lyrs = "combinedSuit")
    suitableArea_ssp585_SH_mid <- rast(fileName_in_ssp585_SH_mid, lyrs = "combinedSuit")
    suitableArea_ssp585_NH_mid <- rast(fileName_in_ssp585_NH_mid, lyrs = "combinedSuit")
    
    #combined
    suitableArea_ssp126_mid <- merge(suitableArea_ssp126_SH_mid, suitableArea_ssp126_NH_mid)
    suitableArea_ssp585_mid <- merge(suitableArea_ssp585_SH_mid, suitableArea_ssp585_NH_mid)
    
    suitableArea_ssp126_end <- merge(suitableArea_ssp126_SH_end, suitableArea_ssp126_NH_end)
    suitableArea_ssp585_end <- merge(suitableArea_ssp585_SH_end, suitableArea_ssp585_NH_end)
    
   commonArea_ssp126_NH_mid <- suitableArea_ssp126_NH_mid * suitableArea_historical_NH 
    commonArea_ssp126_SH_mid <- suitableArea_ssp126_SH_mid * suitableArea_historical_SH
    commonArea_ssp585_NH_mid <- suitableArea_ssp585_NH_mid * suitableArea_historical_NH 
    commonArea_ssp585_SH_mid <- suitableArea_ssp585_SH_mid * suitableArea_historical_SH
    
    commonArea_ssp126_NH_end <- suitableArea_ssp126_NH_end * suitableArea_historical_NH 
    commonArea_ssp126_SH_end <- suitableArea_ssp126_SH_end * suitableArea_historical_SH
    commonArea_ssp585_NH_end <- suitableArea_ssp585_NH_end * suitableArea_historical_NH 
    commonArea_ssp585_SH_end <- suitableArea_ssp585_SH_end * suitableArea_historical_SH
    
    list_initial <- list(speciesChoice, cultivar, CPfruit)
    dt_area_common_mid <- rbind(dt_area_common_mid, append(list_initial, list("NH", "ssp126",  "2041_2060", f_getArea(commonArea_ssp126_NH_mid, layer = 1))))
    dt_area_common_mid <- rbind(dt_area_common_mid, append(list_initial, list("SH", "ssp126",  "2041_2060", f_getArea(commonArea_ssp126_SH_mid, layer = 1))))
    dt_area_common_mid <- rbind(dt_area_common_mid, append(list_initial, list("NH", "ssp585",  "2041_2060", f_getArea(commonArea_ssp585_NH_mid, layer = 1))))
    dt_area_common_mid <- rbind(dt_area_common_mid, append(list_initial, list("SH", "ssp585",  "2041_2060", f_getArea(commonArea_ssp585_SH_mid, layer = 1))))
    
    dt_area_common_end <- rbind(dt_area_common_end, append(list_initial, list("NH", "ssp126",  "2081_2100", f_getArea(commonArea_ssp126_NH_end, layer = 1))))
    dt_area_common_end <- rbind(dt_area_common_end, append(list_initial, list("SH", "ssp126",  "2081_2100", f_getArea(commonArea_ssp126_SH_end, layer = 1))))
    dt_area_common_end <- rbind(dt_area_common_end, append(list_initial, list("NH", "ssp585",  "2081_2100", f_getArea(commonArea_ssp585_NH_end, layer = 1))))
    dt_area_common_end <- rbind(dt_area_common_end, append(list_initial, list("SH", "ssp585",  "2081_2100", f_getArea(commonArea_ssp585_SH_end, layer = 1))))
  } 
  
  dt_area_common <- rbind(dt_area_common_mid, dt_area_common_end)
  dt_area_common[, ssp_year := paste0(ssp, "_", yearSpan)]
  dt_area_common[, yearSpan := NULL]

  # create table of area changes ------
  dt_area <- as.data.table(read.csv(file = paste0(path_data, "areaCalcs", chillLevel, ".csv")))
  dt_area[,ssp_year := paste0(ssp, "_", yearSpan)]
  dt_area[, c("quality", "rasterName",  "yearSpan") := NULL] #"ssp",
  dt_area_wide <-        dcast(dt_area,        species + cultivar + chillPortions + hemisphere ~ ssp_year, value.var = "area_suitable")
  dt_area_common_wide <- dcast(dt_area_common, species + cultivar + chillPortions + hemisphere ~ ssp_year, value.var = "area_common")
  sspNames <- c("ssp126_2041_2060", "ssp126_2081_2100", "ssp585_2041_2060", "ssp585_2081_2100")
  sspNames_new = paste0(sspNames, "_common")
  setnames(dt_area_common_wide, old = sspNames, new = sspNames_new)
  combined <- merge(dt_area_wide, dt_area_common_wide, all = TRUE)
  
  rInArea <- rast(fileName_in)
  harvestArea_earlyCent <- aggregate(rInArea, fact = 6, fun = "sum") # convert 5 arc minutes to 1/2 degrees
  
  combined[, ratioMid2Early_126 := 100 * (-1 + ssp126_2041_2060/historical_1991_2010)]
  combined[, ratioEnd2Early_126 := 100 * (-1 + ssp126_2081_2100/historical_1991_2010)]
  combined[, ratioMid2Early_585 := 100 * (-1 + ssp585_2041_2060/historical_1991_2010)]
  combined[, ratioEnd2Early_585 := 100 * (-1 + ssp585_2081_2100/historical_1991_2010)]
  combined[, ratioMidCommon2Early_126 := 100 * (ssp126_2041_2060_common/historical_1991_2010)]
  combined[, ratioEndCommon2Early_126 := 100 * (ssp126_2081_2100_common/historical_1991_2010)]
  combined[, ratioMidCommon2Early_585 := 100 * (ssp585_2041_2060_common/historical_1991_2010)]
  combined[, ratioEndCommon2Early_585 := 100 * (ssp585_2081_2100_common/historical_1991_2010)]
  combined[, ratioLossMid2Early_126 := 100 * ((historical_1991_2010 - ssp126_2041_2060)/historical_1991_2010)]
  combined[, ratioLossEnd2Early_126 := 100 * ((historical_1991_2010 - ssp126_2081_2100)/historical_1991_2010)]
  combined[, ratioLossMid2Early_585 := 100 * ((historical_1991_2010 - ssp585_2041_2060)/historical_1991_2010)]
  combined[, ratioLossEnd2Early_585 := 100 * ((historical_1991_2010 - ssp585_2081_2100)/historical_1991_2010)]
  ratioColumns <- c("ratioMid2Early_126", "ratioEnd2Early_126", "ratioMid2Early_585", "ratioEnd2Early_585", "ratioMidCommon2Early_126", "ratioEndCommon2Early_126", "ratioMidCommon2Early_585",  "ratioEndCommon2Early_585", "ratioLossMid2Early_126", "ratioLossEnd2Early_126", "ratioLossMid2Early_585", "ratioLossEnd2Early_585")
  sumColumns <- c("historical_1991_2010", "ssp126_2041_2060",  "ssp126_2081_2100", "ssp585_2041_2060", "ssp585_2081_2100", "ssp126_2041_2060_common", "ssp126_2081_2100_common", "ssp585_2041_2060_common", "ssp585_2081_2100_common")
  ssp126Columns <- c("ratioMid2Early_126", "ratioEnd2Early_126", "ratioMidCommon2Early_126", "ratioEndCommon2Early_126", "ratioLossMid2Early_126", "ratioLossEnd2Early_126") 
  combined[, (ratioColumns) := round(.SD, 1), .SDcols = ratioColumns]
  combined[, (sumColumns) := round(.SD, 0), .SDcols = sumColumns]
  combined[, chillPortions := NULL] # column not needed for presentation table
  
  write.csv(combined, file = paste0(path_data, "sumTable", chillLevel, ".csv"), row.names = FALSE)
}

# prepare summary table ------
# main varieties -----
library(flextable)
library(officer)
sumTable_main <- as.data.table(read.csv(file = paste0(path_data, "sumTable", "_main", ".csv")))
sumTable_main[, species := gsub("_main", "", species)]


sumTable_main_flex <- flextable(sumTable_main)

typology_all <- data.frame(
  col_keys_all = names(sumTable_main),
  what = c("Species", "Cultivar", "Hemi-\nsphere", #3
           "Area (sq. km)", "Area (sq. km)","Area (sq. km)","Area (sq. km)","Area (sq. km)","Area (sq. km)","Area (sq. km)","Area (sq. km)","Area (sq. km)","Area (sq. km)","Area (sq. km)","Area (sq. km)","Area (sq. km)", #9
           "Area ratio", "Area ratio", "Area ratio", "Area ratio", "Area ratio", "Area ratio", "Area ratio", "Area ratio"), #15
  measure = c("Species", "Cultivar", "Hemi-\nsphere", #3
              "Historical", "mid century, SSP1-2.6", "end century, SSP1-2.6", "mid century, SSP5-8.5", "end century, SSP5-8.5", #5
              "mid century, SSP1-2.6, common", "end century, SSP1-2.6, common", "mid century, SSP5-8.5, common", "end century, SSP5-8.5, common", #4
              "Mid to early recent historical, SSP1-2.6", "Mid to early recent historical, SSP5-8.5", "end to recent historical century, SSP1-2.6", "end to recent historical century, SSP5-8.5", #4
              "Mid common to early recent historical, SSP1-2.6", "Mid common to early recent historical, SSP5-8.5", "End common to early recent historical, SSP1-2.6", "End common to early recent historical, SSP5-8.5", "Loss of mid to early, 126", "Loss of mid to early, 585", "Loss of end to recent historical, 126", "Loss of end to recent historical, 585"), # 6
  stringsAsFactors = FALSE )

sumTable_main_flex <- set_header_df(sumTable_main_flex, mapping = typology_all, key = "col_keys_all")
sumTable_main_flex <- merge_h(sumTable_main_flex, part = "header")
sumTable_main_flex <- merge_v(sumTable_main_flex, j = c("species", "cultivar", "hemisphere", "ratioMid2Early_585"), part = "header")
sumTable_main_flex <- theme_vanilla(sumTable_main_flex)
sumTable_main_flex <- fix_border_issues(sumTable_main_flex)
sumTable_main_flex <- autofit(sumTable_main_flex)
sumTable_main_flex <- width(sumTable_main_flex, j = 7, width = 1.5)
sumTable_main_flex <- width(sumTable_main_flex, j = 6, width = 1.0)
sumTable_main_flex <- width(sumTable_main_flex, j = 3, width = 0.75)
sumTable_main_flex <- align(sumTable_main_flex, align = "center", i = 1, j = 4, part = "header")
sumTable_main_flex <- height(sumTable_main_flex, height = .25, part = "body")

#area only flextable -----
col_keys_area = names(sumTable_main)[!grepl("ratio", names(sumTable_main), fixed = TRUE)]
sumTable_main_flex_area <- flextable(sumTable_main, col_keys = col_keys_area)
typology_area <- data.frame(
  col_keys = col_keys_area,
  measure = c("Species", "Cultivar", "Hemi-\nsphere", #3
              "Recent historical", "mid century, SSP1-2.6", "end century, SSP1-2.6", "mid century, SSP5-8.5", "end century, SSP5-8.5", 
              "mid century, SSP1-2.6, common", "end century, SSP1-2.6, common", "mid century, SSP5-8.5, common", "end century, SSP5-8.5, common"), # 12
  stringsAsFactors = FALSE )

sumTable_main_flex_area <- set_header_df(sumTable_main_flex_area, mapping = typology_area)
sumTable_main_flex_area <- merge_h(sumTable_main_flex_area, part = "header")
# sumTable_main_flex_area <- merge_v(sumTable_main_flex_area, j = c("species", "cultivar", "hemisphere", "ratioMid2Early_585"), part = "header")
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
sumTable_main_flex_area
save_as_docx(sumTable_main_flex_area, values = NULL, path = "results/flextable_area_main.docx", pr_section = prsect_area)

#ratios only flextable -----
col_keys_ratios = c("species", "cultivar", "hemisphere", names(sumTable_main)[grepl("ratio", names(sumTable_main), fixed = TRUE)])
sumTable_main_flex_ratios <- flextable(sumTable_main, col_keys = col_keys_ratios)
typology_ratios <- data.frame(
  col_keys = col_keys_ratios,
  measure = c("Species", "Cultivar", "Hemi-\nsphere", #3
              "Mid to early recent historical, SSP1-2.6", "end to recent historical century, SSP1-2.6", "Mid to early recent historical, SSP5-8.5", "end to recent historical century, SSP5-8.5", 
              "Mid common to early recent historical, SSP1-2.6", "End common to early recent historical, SSP1-2.6", "Mid common to early recent historical, SSP5-8.5", "End common to early recent historical, SSP5-8.5", 
              "Mid loss of early to SSP1-2.6", "Mid loss of early to SSP5-8.5", "End loss of early to SSP1-2.6", "End loss of early to SSP5-8.5"), #6
  stringsAsFactors = FALSE )

sumTable_main_flex_ratios <- set_header_df(sumTable_main_flex_ratios, mapping = typology_ratios)
sumTable_main_flex_ratios <- merge_h(sumTable_main_flex_ratios, part = "header")
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
sumTable_main_flex_ratios
save_as_docx(sumTable_main_flex_ratios, values = NULL, path = "results/flextable_ratios_main.docx", pr_section = prsect_ratios)

# results for lo chill portions cultivars -----
sumTable_lo <- as.data.table(read.csv(file = paste0(path_data, "sumTable", "_lo", ".csv")))
sumTable_lo[, species := gsub("_lo", "", species)]
sumTable_lo_flex <- flextable(sumTable_lo)

typology_area <- data.frame(
  col_keys = col_keys_area,
  # what = c("Area (sq. km)", "Area (sq. km)","Area (sq. km)", #3
  #          "Area (sq. km)", "Area (sq. km)","Area (sq. km)","Area (sq. km)","Area (sq. km)","Area (sq. km)","Area (sq. km)","Area (sq. km)","Area (sq. km)"), #9
  measure = c("Species", "Cultivar", "Hemi-\nsphere", #3
              "Recent historical", "mid century, SSP1-2.6", "end century, SSP1-2.6", "mid century, SSP5-8.5", "end century, SSP5-8.5", 
              "mid century, SSP1-2.6, common", "end century, SSP1-2.6, common", "mid century, SSP5-8.5, common", "end century, SSP5-8.5, common"), # 12
  stringsAsFactors = FALSE )

sumTable_lo_flex <- set_header_df(sumTable_lo_flex, mapping = typology_all, key = "col_keys_all")
sumTable_lo_flex <- merge_h(sumTable_lo_flex, part = "header")
sumTable_lo_flex <- merge_v(sumTable_lo_flex, j = c("species", "cultivar", "hemisphere", "ratioMid2Early_585"), part = "header")
sumTable_lo_flex <- theme_vanilla(sumTable_lo_flex)
sumTable_lo_flex <- fix_border_issues(sumTable_lo_flex)
sumTable_lo_flex <- autofit(sumTable_lo_flex)
sumTable_lo_flex <- fit_to_width(sumTable_lo_flex, max_width = 10, inc = 1L, max_iter = 20, unit = "in")

sumTable_lo_flex <- width(sumTable_lo_flex, j = 7, width = 1.5)
sumTable_lo_flex <- width(sumTable_lo_flex, j = 6, width = 1.0)
sumTable_lo_flex <- width(sumTable_lo_flex, j = 3, width = 0.75)

sumTable_lo_flex <- align(sumTable_lo_flex, align = "center", i = 1, j = 4, part = "header")
sumTable_lo_flex <- height(sumTable_lo_flex, height = .25, part = "body")

#area only flextable -----
col_keys_area = names(sumTable_lo)[!grepl("ratio", names(sumTable_lo), fixed = TRUE)]
sumTable_lo_flex_area <- flextable(sumTable_lo, col_keys = col_keys_area)
typology_area <- data.frame(
  col_keys = col_keys_area,
  measure = c("Species", "Cultivar", "Hemi-\nsphere", #3
              "Recent historical", "mid century, SSP1-2.6", "end century, SSP1-2.6", "mid century, SSP5-8.5", "end century, SSP5-8.5", 
              "mid century, SSP1-2.6, common", "end century, SSP1-2.6, common", "mid century, SSP5-8.5, common", "end century, SSP5-8.5, common"), # 12
  stringsAsFactors = FALSE )

sumTable_lo_flex_area <- set_header_df(sumTable_lo_flex_area, mapping = typology_area)
sumTable_lo_flex_area <- merge_h(sumTable_lo_flex_area, part = "header")
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
sumTable_lo_flex_area

save_as_docx(sumTable_lo_flex_area, values = NULL, path = "results/flextable_area_lo.docx", pr_section = prsect_area)

#ratios only flextable -----
col_keys_ratios = c("species", "cultivar", "hemisphere", names(sumTable_lo)[grepl("ratio", names(sumTable_lo), fixed = TRUE)])
sumTable_lo_flex_ratios <- flextable(sumTable_lo, col_keys = col_keys_ratios)
typology_ratios <- data.frame(
  col_keys = col_keys_ratios,
  measure = c("Species", "Cultivar", "Hemi-\nsphere", #3
              "Mid to recent historical, SSP1-2.6", "End to recent historical, SSP1-2.6", "Mid to recent historical, SSP5-8.5", "End to recent historical, SSP5-8.5", 
              "Mid common to recent historical, SSP1-2.6", "End common to recent historical, SSP1-2.6", "Mid common to recent historical, SSP5-8.5", "End common to recent historical, SSP5-8.5", 
              "Mid loss of early to SSP1-2.6", "Mid loss of early to SSP5-8.5", "End loss of early to SSP1-2.6", "End loss of early to SSP5-8.5"), #6
  stringsAsFactors = FALSE )

sumTable_lo_flex_ratios <- set_header_df(sumTable_lo_flex_ratios, mapping = typology_ratios)
sumTable_lo_flex_ratios <- merge_h(sumTable_lo_flex_ratios, part = "header")
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
sumTable_lo_flex_ratios
save_as_docx(sumTable_lo_flex_ratios, values = NULL, path = "results/flextable_ratios_lo.docx", pr_section = prsect_ratios)
