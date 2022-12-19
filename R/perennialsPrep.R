#  creates the data tables majorCropValues_main, majorCropValues_lo and majorCropValues_hi
library(readxl)
library(data.table)
supp_materials_chill_portions <- as.data.table(read_excel("data-raw/perennials/supp_materials_2021_06_22.xlsx")) #, col_types = c("text", "text", "numeric", "text", "text")))
setnames(supp_materials_chill_portions, old = names(supp_materials_chill_portions), 
         new = c("cropName", "cultivar", "chill_requirement", "comment", "frost_threshold", "low_temp_threshold", "chill_hours", "summer_heat_threshold", "gdd", "gddtb", "GDD_opt", "reference_chill_portions", "reference_other_information", "reference_gdd", "other_comments"))
supp_materials_chill_portions[, cropName := tolower(cropName)]

# get rid of blueberries row
supp_materials_chill_portions <- supp_materials_chill_portions[!cropName == "blueberries"]
test <- data.table::copy(supp_materials_chill_portions)
test <- test[, CR_cultivar_mean := round(mean(chill_requirement), 1), by = "cultivar"]
test <- test[, CR_cultivar_min := min(chill_requirement, na.rm = TRUE), by = "cultivar"]
test <- test[, CR_cultivar_max := max(chill_requirement, na.rm = TRUE), by = "cultivar"]
test <- test[, CR_crop_mean := round(mean(chill_requirement, na.rm = TRUE), 1), by = "cropName"]
test <- test[, CR_crop_min := round(min(chill_requirement, na.rm = TRUE), 1), by = "cropName"]
test <- test[, CR_crop_max := round(max(chill_requirement, na.rm = TRUE), 1), by = "cropName"]
test[CR_crop_max =="Inf" | CR_crop_max =="-Inf", CR_crop_max := NA]
test[CR_crop_min =="Inf" | CR_crop_min =="-Inf", CR_crop_min := NA]

# added July 6, 2021
test[, frost_threshold := -2] # a general value for all plants. They can survive some period of temps below 0
# adjust GDD_opt to be the same for all crops
test[, GDD_opt := 25]

#adjust gddtb
test[cropName == "almond", gddtb := 4.5]
test[cropName == "apple", gddtb := 5]
test[cropName == "cherry", gddtb := 5]
test[cropName == "olive", gddtb := 10]
#test[cropName == "winegrape", gddtb := 10]
test[cropName == "grape", gddtb := 10]
# ---- end of adjustments of July 6, 2021
# use lower gdd value of 1100 for Chardonnay grapes and the artificial no-chill variety
test[cultivar %in% c("Chardonnay", "No_chill"), gdd := 1100]

# remove extraneous columns
test [, c("chill_requirement", "comment", "reference_chill_portions", "reference_other_information", "other_comments", "chill_hours", "reference_gdd") := NULL]
test <- unique(test)
dominantVarieties <- c("Gala", "Nonpareil", "Lapins", "Picual", "Northern highbush", "Chardonnay") # most widely grown
majorCropValues_main <- data.table::copy(test)
majorCropValues_main <- majorCropValues_main[cultivar %in% dominantVarieties,]
majorCropValues_main[, chill_portions := CR_cultivar_mean]
majorCropValues_lo <- data.table::copy(test)
majorCropValues_lo[, chill_portions := CR_crop_min]
majorCropValues_hi <-data.table::copy(test)
majorCropValues_hi[, chill_portions := CR_crop_max]

majorCropValues_hi[, cropName := lapply(.SD, paste0, "_hi"), .SDcols = "cropName"]
majorCropValues_lo[, cropName := lapply(.SD, paste0, "_lo"), .SDcols = "cropName"]
majorCropValues_main[, cropName := lapply(.SD, paste0, "_main"), .SDcols = "cropName"]
majorCropValues_lo <- majorCropValues_lo[CR_cultivar_min == chill_portions, ]
majorCropValues_main <- majorCropValues_main[CR_cultivar_mean == chill_portions, ]
majorCropValues_hi <- majorCropValues_hi[CR_cultivar_max == chill_portions, ]

deleteListCol <- c("CR_cultivar_min", "CR_cultivar_max", "CR_crop_mean", "CR_crop_min", "CR_crop_max", "CR_cultivar_mean") 
majorCropValues_main[, (deleteListCol) := NULL]
majorCropValues_lo[, (deleteListCol) := NULL]
majorCropValues_hi[, (deleteListCol) := NULL]
write.csv(majorCropValues_main, "data/perennials/majorCropValues_main.csv", row.names = FALSE)
write.csv(majorCropValues_lo, "data/perennials/majorCropValues_lo.csv", row.names = FALSE)
write.csv(majorCropValues_hi, "data/perennials/majorCropValues_hi.csv", row.names = FALSE)

