library(terra)
library(data.table)
path_data <- "data/perennials/"
speciesChoices <- c("almond_main", "apple_main",  "cherry_main", "olive_main",  "grape_main")
speciesNames <- gsub("_main", "", speciesChoices)
dominantVarieties <- c("Gala", "Nonpareil", "Lapins", "Picual", "Northern highbush", "Chardonnay") # most widely grown
loChillVarieties <- c("Achak", "Eva", "17GE580", "Chemlali", "No_chill")
extent_noAntarctica <- ext(-180, 180, -60, 90)
hemispheres <- c("NH", "SH")
yearRange <- 19
f_getArea <- function(r, layer) {
  r_sub <- subset(r, layer)
  r_sub[r_sub < 1 | r_sub > 1] <- NA # keep only locations with value of 1
  r_sub_area <- (expanse(r_sub, unit = "km"))/1000 # convert to 1000 sq km
  r_sub_area <- round(r_sub_area, 0)
}

{
  thermal_coeffs <- as.data.table(readxl::read_excel("data-raw/perennials/supp_materials_2021_06_22.xlsx", 
                                                     col_types = c("text", "text", "numeric", 
                                                                   "skip", "numeric", "numeric", "skip", 
                                                                   "numeric", "numeric", "numeric", 
                                                                   "numeric", "skip", "skip", "skip", 
                                                                   "skip")))
  setnames(thermal_coeffs, old = names(thermal_coeffs), new = c("cropName", "cultivar", "chill_requirement", "chill_threshold", "low_temp_threshold", "summer_heat_threshold", "GDD", "GDD_Tb", "GDD_topt"))
  thermal_coeffs[, cropName := tolower(cropName)]
  thermal_coeffs <- thermal_coeffs[cropName %in% speciesNames,]
  
  thermal_coeffs[, CR_cultivar_mean := round(mean(chill_requirement), 1), by = "cultivar"]
  thermal_coeffs[, CR_cultivar_mean := round(mean(chill_requirement), 1), by = "cultivar"]
  thermal_coeffs[, CR_cultivar_min := min(chill_requirement, na.rm = TRUE), by = "cultivar"]
  thermal_coeffs[, CR_cultivar_max := max(chill_requirement, na.rm = TRUE), by = "cultivar"]
  thermal_coeffs[, CR_crop_mean := round(mean(chill_requirement, na.rm = TRUE), 1), by = "cropName"]
  thermal_coeffs[, CR_crop_min := round(min(chill_requirement, na.rm = TRUE), 1), by = "cropName"]
  thermal_coeffs[, CR_crop_max := round(max(chill_requirement, na.rm = TRUE), 1), by = "cropName"]
  thermal_coeffs[, chill_requirement := NULL]
  thermal_coeffs <- unique(thermal_coeffs)
  majorCropValues_main <- data.table::copy(thermal_coeffs)
  majorCropValues_main <- majorCropValues_main[cultivar %in% dominantVarieties,]
  majorCropValues_main[, chill_portions := CR_cultivar_mean]
  majorCropValues_main <- unique(majorCropValues_main)
  majorCropValues_lo <- data.table::copy(thermal_coeffs)
  majorCropValues_lo <- majorCropValues_lo[cultivar %in% loChillVarieties,]
  majorCropValues_lo[, chill_portions := CR_crop_min]
  majorCropValues_lo <- unique(majorCropValues_lo)
  
  for (chillLevel in c("_lo", "_main")) { 
    namelist <- c("species", "cultivar", "chillPortions", "hemisphere", "quality", "ssp", "yearSpan",  "area_suitable", "rasterName")
    dt_area <- data.table(1)[,`:=`((namelist),NA)][,V1:=NULL][.0]
    cropVals <- get(paste0("majorCropValues", chillLevel))
    #scenarios 
    for (speciesName in speciesNames) { # this construction to allow speciesChoice to include _lo
      speciesChoice <- paste0(speciesName, chillLevel)
      CPfruit <- mean(cropVals[cropName == speciesName, chill_portions])
      cultivar <-  cropVals[cropName == speciesName, cultivar]
      for (hem in hemispheres) {
        for (k in c("ssp126", "ssp585")) {
          l <- 2081
          yearSpan <- paste0(l, "_", l + yearRange)
          inF <- paste0(path, "nonlimiting_all_", speciesChoice, "_", k, "_", "good", "_", hem, "_", yearSpan, ".tif")
          rastName <- paste0("suitable_", speciesChoice, "_", hem, "_", k, "_", "good", "_", yearSpan)
          # r has all 5 suitability layers, get just the combined one
          r_combSuit <- rast(inF, lyrs = 1 )
          r_combSuit[r_combSuit == 0] <- NA # only want to get area of cells that have a value of 1; ie, suitable
          r_area <- f_getArea(r_combSuit, layer = 1)
          print(paste0("speciesChoice: ", speciesChoice, ", cultivar: ", cultivar, ", ssp: ", k, ",start year: ", l, ", hemisphere: ", hem, ", chillVal: ", chillLevel, ", r_area: ", r_area))
          dt_area <- rbind(dt_area, list(speciesChoice, cultivar, CPfruit, hem, "good", k, yearSpan,  r_area, rastName))
        }
        #historical
        k <- "historical"
        l <- 1991
        yearSpan <- paste0(l, "_", l + yearRange)
        inF <- paste0(path, "nonlimiting_all_", speciesChoice, "_", k, "_", "good", "_", hem, "_", yearSpan, ".tif")
        rastName <- paste0("suitable_", speciesChoice, "_", hem, "_", k, "_", "good", "_", yearSpan)
        # r has all 5 suitability layers, get just the combined one
        r_combSuit <- rast(inF, lyrs = 1 )
        r_combSuit[r_combSuit == 0] <- NA # only want to get area of cells that have a value of 1; ie, suitable
        r_area <- f_getArea(r_combSuit, layer = 1)
        print(paste0("speciesChoice: ", speciesChoice, ", cultivar: ", cultivar, ", ssp: ", k, ",start year: ", l, ", hemisphere: ", hem, ", chillVal: ", chillLevel, ", r_area: ", r_area))
        dt_area <- rbind(dt_area, list(speciesChoice, cultivar, CPfruit, hem, "good", k, yearSpan, r_area, rastName))
      }
    }
    outf <- paste0(path, "areaCalcs", chillLevel, ".csv")
    write.csv(dt_area, file = outf, row.names = FALSE)
    print(paste0("fileName out: ", outf)) # input into table 7
  }
  
  namelist <- c("species", "cultivar", "chillPortions", "hemisphere", "quality", "ssp", "yearSpan", "area_base", "area_ssp", "area_both", "area_hist_loss", " area_ssp_gain", "area_ssp_gain_lo")
  dt_area_delta <- data.table(1)[,`:=`((namelist),NA)][,V1:=NULL][.0]
  for (chillLevel in c("_lo", "_main")) {
    cropVals <- get(paste0("majorCropValues", chillLevel))
    for (speciesName in speciesNames) { # this construction to allow speciesChoise to include _lo
      speciesChoice <- paste0(speciesName, chillLevel)
      CPfruit <- cropVals[cropName == speciesName, chill_portions]
      cultivar <-  cropVals[cropName == speciesName, cultivar]
      print(paste0("cropname: ", speciesName, ", cultivar: ", cultivar))
      
      for (hem in hemispheres) {
        filename_r_historical_in <- paste0(path, "nonlimiting_all_", paste0(speciesName, chillLevel), "_", "historical", "_", "good", "_", hem, "_", "1991_2010.tif")
        r_historical_combinedSuit <- rast(filename_r_historical_in, lyrs = "combinedSuit") # combined suitability means all of the temp metrics are satisfied (+1)
        r_base <- r_historical_combinedSuit; r_base[r_base < 1] <- NA
        area_base <- f_getArea(r_base, layer = 1)
        # do scenarios and deltas
        # use just end of century and ssp585
        l <- 2081; 
        k <- "ssp585"
        yearSpan <- paste0(l, "_", l + yearRange)
        filename_r_ssp_in <- paste0(path, "nonlimiting_all_", paste0(speciesName, chillLevel), "_",  k, "_", "good", "_", hem, "_", "2081_2100", ".tif")
        r_ssp_combinedSuit <- rast(filename_r_ssp_in, lyrs = "combinedSuit")
        r_sum <- r_ssp_combinedSuit + r_historical_combinedSuit # 2 - in both periods, 1 in 1 of the periods, 0 not in either
        r_both <- r_sum; r_both[!r_both == 2] <- NA # get area only for pixels with value of 2
        r_both <- r_both/2 # convert to 1/NA
        r_ssp <- r_ssp_combinedSuit; r_ssp[!r_ssp == 1] <- NA
        r_delta <- r_ssp_combinedSuit - r_historical_combinedSuit # 1 in 1 of the periods, not the other, 0 not in either or in both, -1 in historical but not in the future
        r_hist_loss <- r_delta; r_hist_loss[!r_hist_loss == -1] <- NA # get pixels only with the value of -1, lost historical area
        r_hist_loss <- r_hist_loss * -1 # convert the -1s to 1s
        r_ssp_gain <- r_delta; r_ssp_gain[!r_ssp_gain == 1 ] <- NA # new pixels from SSP; remove historical pixels lost (-1) and pixels where either both are 1 or both are 0
        
        fileName_r_end_lo_in <- paste0(path, "nonlimiting_all_", paste0(speciesName, "_lo"), "_",  k, "_", "good", "_", hem, "_", "2081_2100", ".tif")
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
        
        dt_area_delta <- rbind(dt_area_delta, list(speciesChoice, cultivar, CPfruit, hem, "good", k, yearSpan, area_base, area_ssp, area_both, area_hist_loss, area_ssp_gain, area_ssp_gain_lo))
      }
    }
  }
  outf <- paste0(path, "areaCalcs_delta", ".csv")
  write.csv(dt_area_delta, file = outf, row.names = FALSE)
  print(paste0("fileName out: ", outf))
}

# harvest and suitability area calcs ------
# area common to mid and end century -----
for (chillLevel in c("_lo", "_main")) { 
  namelist <- c("species", "cultivar", "chillPortions", "hemisphere",  "ssp", "yearSpan", "area_common")
  dt_area_common_end <-  data.table(1)[,`:=`((namelist),NA)][,V1 := NULL][.0]
  dt_area_common_end[, area_common := as.numeric(area_common)]
  
  for (speciesName in speciesNames) { # this construction to allow speciesChoice to include _lo
    speciesChoice <- paste0(speciesName, chillLevel)
    cropVals <- get(paste0("majorCropValues", chillLevel))
    CPfruit <- mean(cropVals[cropName == speciesName, chill_portions])
    cultivar <-  cropVals[cropName == speciesName, cultivar]
    
    #suitable areas -----
    fileName_start <- paste0(path, "nonlimiting_all_", speciesChoice, "_")
    
    # suitable historical area
    inF_historical_SH <- paste0(fileName_start, "historical", "_", "good", "_", "SH", "_", "1991_2010", ".tif")
    inF_historical_NH <- paste0(fileName_start, "historical", "_", "good", "_", "NH", "_", "1991_2010", ".tif")
    suitableArea_historical_SH <- rast(inF_historical_SH, lyrs = "combinedSuit")
    suitableArea_historical_NH <- rast(inF_historical_NH, lyrs = "combinedSuit")
    
    inF_ssp126_SH_end <- paste0(fileName_start, "ssp126", "_", "good", "_", "SH", "_", "2081_2100", ".tif")
    inF_ssp126_NH_end <- paste0(fileName_start, "ssp126", "_", "good", "_", "NH", "_", "2081_2100", ".tif")
    inF_ssp585_SH_end <- paste0(fileName_start, "ssp585", "_", "good", "_", "SH", "_", "2081_2100", ".tif")
    inF_ssp585_NH_end <- paste0(fileName_start, "ssp585", "_", "good", "_", "NH", "_", "2081_2100", ".tif")
    
    suitableArea_ssp126_SH_end <- rast(inF_ssp126_SH_end, lyrs = "combinedSuit")
    suitableArea_ssp126_NH_end <- rast(inF_ssp126_NH_end, lyrs = "combinedSuit")
    suitableArea_ssp585_SH_end <- rast(inF_ssp585_SH_end, lyrs = "combinedSuit")
    suitableArea_ssp585_NH_end <- rast(inF_ssp585_NH_end, lyrs = "combinedSuit")
    
     commonArea_ssp126_NH_end <- suitableArea_ssp126_NH_end * suitableArea_historical_NH 
    commonArea_ssp126_SH_end <- suitableArea_ssp126_SH_end * suitableArea_historical_SH
    commonArea_ssp585_NH_end <- suitableArea_ssp585_NH_end * suitableArea_historical_NH 
    commonArea_ssp585_SH_end <- suitableArea_ssp585_SH_end * suitableArea_historical_SH
    
    list_initial <- list(speciesChoice, cultivar, CPfruit)
    dt_area_common_end <- rbind(dt_area_common_end, append(list_initial, list("NH", "ssp126",  "2081_2100", f_getArea(commonArea_ssp126_NH_end, layer = 1))))
    dt_area_common_end <- rbind(dt_area_common_end, append(list_initial, list("SH", "ssp126",  "2081_2100", f_getArea(commonArea_ssp126_SH_end, layer = 1))))
    dt_area_common_end <- rbind(dt_area_common_end, append(list_initial, list("NH", "ssp585",  "2081_2100", f_getArea(commonArea_ssp585_NH_end, layer = 1))))
    dt_area_common_end <- rbind(dt_area_common_end, append(list_initial, list("SH", "ssp585",  "2081_2100", f_getArea(commonArea_ssp585_SH_end, layer = 1))))
  } 
  
  dt_area_common <- dt_area_common_end
  dt_area_common[, ssp_year := paste0(ssp, "_", yearSpan)]
  dt_area_common[, yearSpan := NULL]
  
  # create table of area changes ------
  dt_area <- as.data.table(read.csv(file = paste0(path, "areaCalcs", chillLevel, ".csv")))
  dt_area[,ssp_year := paste0(ssp, "_", yearSpan)]
  dt_area[, c("quality", "rasterName",  "yearSpan") := NULL]
  dt_area_wide <-        dcast(dt_area,        species + cultivar + chillPortions + hemisphere  ~ ssp_year, value.var = "area_suitable")
  dt_area_common_wide <- dcast(dt_area_common, species + cultivar + chillPortions + hemisphere ~ ssp_year, value.var = "area_common")
  sspNames <- c("ssp126_2081_2100", "ssp585_2081_2100")
  sspNames_new = paste0(sspNames, "_common")
  setnames(dt_area_common_wide, old = sspNames, new = sspNames_new)
  combined <- merge(dt_area_wide, dt_area_common_wide, all = TRUE)
  
  combined[, ratioEnd2Early_126 := 100 * (-1 + ssp126_2081_2100/historical_1991_2010)]
  combined[, ratioEnd2Early_585 := 100 * (-1 + ssp585_2081_2100/historical_1991_2010)]
  combined[, ratioEndCommon2Early_126 := 100 * (ssp126_2081_2100_common/historical_1991_2010)]
  combined[, ratioEndCommon2Early_585 := 100 * (ssp585_2081_2100_common/historical_1991_2010)]
  combined[, ratioLossEnd2Early_126 := 100 * ((historical_1991_2010 - ssp126_2081_2100)/historical_1991_2010)]
  combined[, ratioLossEnd2Early_585 := 100 * ((historical_1991_2010 - ssp585_2081_2100)/historical_1991_2010)]
    sumColumns <- c("historical_1991_2010",  "ssp126_2081_2100", "ssp585_2081_2100", "ssp126_2081_2100_common", "ssp585_2081_2100_common")
  ssp126Columns <- c("ratioEnd2Early_126", "ratioEndCommon2Early_126", "ratioLossEnd2Early_126") 
  combined[, (sumColumns) := round(.SD, 0), .SDcols = sumColumns]
  out_f <- paste0(path, "sumTable", chillLevel, ".csv")
  write.csv(combined, file = out_f, row.names = FALSE)
  print(paste0("out file: ", out_f))
}

# prepare summary table ------
# main varieties -----
library(flextable)
library(officer)
sumTable_main <- as.data.table(read.csv(file = paste0(path, "sumTable", "_main", ".csv")))
sumTable_lo <- as.data.table(read.csv(file = paste0(path, "sumTable", "_lo", ".csv")))
sumTable <- rbind(sumTable_main, sumTable_lo)
sumTable[, species := gsub("_main", "", species)][, species := gsub("_lo", "", species)]
setorderv(sumTable, c("species", "hemisphere"))

#area only flextable -----
col_keys_area = names(sumTable_main)[!grepl("ratio", names(sumTable_main), fixed = TRUE)]
sumTable_main_flex_area <- flextable(sumTable_main, col_keys = col_keys_area)
typology_area <- data.frame(
  col_keys = col_keys_area,
  measure = c("Species", "Cultivar", "Chill portions", "Hemi-\nsphere", #3
              "Recent historical", "End century, SSP1-2.6", "End century, SSP5-8.5", 
              "End century, \nSSP1-2.6, \ncommon", "End century, \nSSP5-8.5, \ncommon"), 
  stringsAsFactors = FALSE )

sumTable_main_flex_area <- set_header_df(sumTable_main_flex_area, mapping = typology_area)
sumTable_main_flex_area <- merge_h(sumTable_main_flex_area, part = "header")
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

