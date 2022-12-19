library(readxl)
library(data.table)
supp_materials_2021_06_22 <- as.data.table(read_excel("data-raw/perennials/supp_materials_2021_06_22.xlsx"))

chillPortions <- supp_materials_2021_06_22[, c("Crop", "Cultivar", "Chill requirement")]
setnames(chillPortions, old = c("Crop", "Cultivar", "Chill requirement"), new = c("crop", "cultivar", "chillRequirement"))
chillPortions <- chillPortions[!crop == "blueberries", ]
chillPortions[, minValByCrop := min(chillRequirement), by = crop]
chillPortions[, maxValByCrop := max(chillRequirement), by = crop]
chillPortions[, c("cultivar", "chillRequirement") := NULL]
chillPortions <- unique(chillPortions)
write.csv(chillPortions, "data/chillPortions/chillPortionRange.csv")
