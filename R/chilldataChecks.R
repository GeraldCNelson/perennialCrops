# file data correct range?
fileLoc_out <- '../data/tasmin/'
badFileList <- character()
files_completed <- list.files(fileLoc_out, pattern = '.tif', recursive = TRUE, full.names = TRUE)
files_completed <- files_completed[!grepl("aux.xml", files_completed, fixed = TRUE)]
files_completed <- gsub("//", "/", files_completed)
files_completed <- unique(files_completed)
for (i in files_completed) {
  message(i)
  r <- rast(i)
  testVal_min <- min(minmax(r))
  testVal_max <- max(minmax(r))
  
  if(testVal_min < -75 | testVal_max > 100) {   
    print(paste0("min: ", testVal_min, ", max: ", testVal_max, ", file: ", i))    
    badfileList <- c(badFileList, i)   
  }  
}
# file exists?
yearRange <- 19
startYearChoices <- c(1991, 2041, 2081)
climateVars <- c("tasmin", "tasmax")
sspChoices <- c("ssp585", "ssp126")
modelChoices <- c( "GFDL-ESM4", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL", "IPSL-CM6A-LR")
modelChoices_lower <- tolower(modelChoices)
for (k in sspChoices) {
  
  for (l in startYearChoices) {    
    yearSpan <- paste0(l, "_", l + yearRange)   
    for (climateVar in climateVars) {
      for (modelChoice_lower in modelChoices_lower) {
        fileName_in <- paste0("../data/", climateVar, '/', modelChoice_lower, "_", climateVar, "_", k, "_", yearSpan, ".tif")
        if (!file.exists(fileName_in)) print(fileName_in)
        #        r <- rast(fileName_in)
      }     
    }   
  }
}
