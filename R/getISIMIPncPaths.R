this <- system('hostname', TRUE)
if (this == "LAPTOP-IVSPBGCA") {
  setwd("G:/.shortcut-targets-by-id/1mfeEftF_LgRcxOT98CBIaBbYN4ZHkBr_/share/pwc")
} else {
  setwd('/Users/gcn/Google Drive/My Drive/pwc')
}

library(data.table)

path <- "data-raw/ISIMIP/filelists/"

holder <- data.table(p = character())
temp <- list.files(path = path)
for (i in 1:length(temp)) {
  test <- read.csv(paste0(path,temp[i]))
  print(temp[i])
  names(test) <- "p"
  holder <- rbind(holder, test)
}

temp <- copy(holder)
temp <- temp[p %like% "1991|2001|2041|2051|2181|2091"]
temp <- temp[p %like% "_tas_|_tasmax_|tasmin"] # need to update for the PWC paper
temp <- temp[p %like% "_historical_|_ssp126_|ssp585"]
write.csv(temp, "/Users/gcn/Documents/workspace/perennialCrops/data-raw/ISIMIPNCfiles.csv")
