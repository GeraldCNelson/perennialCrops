#  chilling hours calculations
source("R/globallyUsed.R")
#library(doParallel) #Foreach Parallel Adaptor 
# library(foreach) #Provides foreach looping construct, called with doParallel

locOfFiles <- locOfCMIP6ncFiles
sspChoices <- c("ssp585") #"ssp126", 
modelChoices <- c("IPSL-CM6A-LR")# 

modelChoices <- c("IPSL-CM6A-LR", "GFDL-ESM4", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL") 

startYearChoices <-  c(2021, 2051, 2091) 

yearRange <- 9

# f.chillhrs <- function(tmin, tmax) {
#   ch <- (7 - tmin)/(tmax - tmin)
#   ch[tmin > 7] <- 0
#   ch[tmax < 7 & tmin <= 7] <- 24
#   ch
# }


#test values
i <- "IPSL-CM6A-LR"
k <- "ssp585"
l <- 2091
northernHemWinter <- c("Nov", "Dec", "Jan", "Feb") #, "Mar", "Apr")
northernHemSummer <- c("May", "Jun", "Jul", "Aug", "Sep", "Oct")
southernHemWinter <- c("May", "Jun", "Jul", "Aug") #, "Sep", "Oct")
southernHemSummer <- c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr")
# useCores <- detectCores() - 2 # max number of cores
# useCores <- 2 # better for memory intensive activities


library(Rcpp)
cppFunction('std::vector<double> chill(std::vector<double> tmin, std::vector<double> tmax) {
size_t n = tmin.size();
std::vector<double> out(n);
for (size_t i=0; i<n; i++) {
if (tmax[i] < 7) {
out[i] = 24;
} else if (tmin[i] < 7) {
out[i] = (7 - tmin[i])/(tmax[i] - tmin[i]);
}
}
return out;
}')

# cl <- clusterSetup(varList, libList, useCores) # function created in globallyUsed.R
# foreach(l = startYearChoices) %:%
#   foreach(i = modelChoices) %:%
#   #  foreach(j = variableChoices) %:%
#   foreach(k = sspChoices) %dopar% {
for (l in startYearChoices) {
  for (i in modelChoices) {
    for (k in sspChoices) {
      
      print(paste0("start year: ", l, " ssp: ", k, " pid: ", Sys.getpid(), " systime: ", Sys.time()))
      modelName.lower <- tolower(i)
      startTime <-  Sys.time()
      yearSpan <- paste0(l, "_", l + yearRange)
      fileName_in <- paste(modelName.lower, k, "tasmax", "global_daily", yearSpan, sep = "_", ".tif")
      tmaxIn <- paste0(locOfFiles, k,"/", i, "/", fileName_in)
      
      fileName_in <- paste(modelName.lower, k, "tasmin", "global_daily", yearSpan, sep = "_", ".tif")
      tminIn <- paste0(locOfFiles, k,"/", i, "/", fileName_in)
      tmaxTminIn(tmaxIn, tminIn) # function to read in tmax and tmin with rast
      terra:::.mem_info(tmax, 1)
      tmp <- sds(tmin, tmax)
      
      print(system.time(chillHrs <- lapp(tmp, chill)))
      # tbase <- 7
      # system.time(chillHrs <- clamp(tmin, lower = tbase, upper = tbase_max, values = TRUE))
      names(chillHrs) <- names(tmax) # put the date info back into the names
      # 
      # do several count days in a month
      # first days with temp below zero
      print(paste0("Done with chillHrs function", " pid: ", Sys.getpid()))
      startDate <- paste0(l, "-01-01")
      endDate <- paste0(l + yearRange, "-12-31")
      indices <- seq(as.Date(startDate), as.Date(endDate), 1)
      indices <- format(indices, format = "%m")
      indices <- as.numeric(indices)
      print(system.time(monthZeroCount <- tapp(tmin, indices, fun = function(x, ...){sum(x <= 0)}, na.rm = TRUE)))
      names(monthZeroCount) <- month.abb
      fileName_outZero <- paste0("belowZeroCount_", modelName.lower, "_", k, "_", yearSpan, ".tif")
      print(paste0("Writing out ", fileName_outZero))
      writeRaster(monthZeroCount, filename = paste0("data/belowZero/", fileName_outZero),  overwrite = TRUE)
      
      #   rm(list = c("tmax", "tmin"))
      print(system.time(chillHrs.sumMonth <- tapp(chillHrs, indices, fun = sum, na.rm = TRUE)))
      chillHrs.sumMonth <- chillHrs.sumMonth/10 # to get to the monthly average over 10 years
      names(chillHrs.sumMonth) <- month.abb
      chillHrsNorthernHem <- subset(chillHrs.sumMonth, southernHemWinter) # note dropping layers for southern hemisphere winter
      chillHrsSouthernHem <- subset(chillHrs.sumMonth, northernHemWinter) # note dropping layers for northern hemisphere winter
      chillHrsNorthernHem <- sum(chillHrsNorthernHem)
      chillHrsSouthernHem <- sum(chillHrsSouthernHem)
      
      #print(endCompleteLoop - startTime)
      
      fileNameNH <- paste0("chillHrs_NorthernHem_", modelName.lower, "_", k, "_", yearSpan, ".tif")
      fileNameSH <- paste0("chillHrs_SouthernHem_", modelName.lower, "_", k, "_", yearSpan, ".tif")
      
      writeRaster(chillHrsNorthernHem, filename = paste0("data/chillingHours/", fileNameNH),  overwrite = TRUE)
      writeRaster(chillHrsSouthernHem, filename = paste0("data/chillingHours/", fileNameSH),  overwrite = TRUE)
    }
  }
}
#stopCluster(cl)

# do same calculations on observed data
tmin <- tasmin.observed
tmax <- tasmax.observed
yearSpan <- paste0(l, "_", l + yearRange)

tmp <- sds(tmin, tmax)
print(system.time(chillHrs <- lapp(tmp, chill)))

#chillHrs <- overlay(tmin, tmax, fun = f.chillhrs)
print("Done with chillHrs function")
startDate <- paste0(2001, "-01-01")
endDate <- paste0(l + yearRange, "-12-31")
indices <- seq(as.Date(startDate), as.Date(endDate), 1)
indices <- format(indices, format = "%m")
indices <- as.numeric(indices)

# do several count days in a month

# first count the number of days with temp below zero
system.time(monthZeroCount <- tapp(tmin, indices, fun = function(x, ...){sum(x <= 0)}, na.rm = TRUE))
names(monthZeroCount) <- month.abb
fileName_outZero <- paste0("belowZeroCount", "_observed_", yearSpan, ".tif")
writeRaster(monthZeroCount, filename = paste0("data/belowZero/", fileName_outZero),  overwrite = TRUE)

# # now do count above tmax limit
# f.tmaxLimit <- function(tmax, tmaxLimit, indices) {
#   tmaxSum <- tapp(tmax, indices, fun = function(x, ...){sum(x >= tmaxLimit)}) 
#   names(tmaxSum) <- month.abb
#   fileName_out <- paste0("tmaxGT_", tmaxLimit, "_observed_", yearSpan, ".tif")
#   writeRaster(tmaxSum, filename = paste0("data/tmaxMonthlySums/", fileName_out),  overwrite = TRUE)
# }
# tmaxfunctionStart <- Sys.time()
# #tmax > 31
# f.tmaxLimit(tmax, tmaxLimit = 31, indices)
# print(paste("Completed tmaxlimit for 31°C"))
# #tmax > 35
# f.tmaxLimit(tmax, tmaxLimit = 35, indices)
# #tmax > 38
# f.tmaxLimit(tmax, tmaxLimit = 38, indices)
# #tmax > 45
# f.tmaxLimit(tmax, tmaxLimit = 45, indices)
# #tmax > 48
# f.tmaxLimit(tmax, tmaxLimit = 48, indices)
# print(paste("Completed tmaxlimit for 48°C"))

rm(list = c("tmax", "tmin"))
chillHrs.sumMonth <- tapp(chillHrs, indices, fun = sum, na.rm = TRUE)
chillHrs.sumMonth <- chillHrs.sumMonth/10 # to get to the monthly average over 10 years
names(chillHrs.sumMonth) <- month.abb
chillHrsNorthernHem <- dropLayer(chillHrs.sumMonth, southernHemWinter) # note dropping layers for southern hemisphere winter
chillHrsSouthernHem <- dropLayer(chillHrs.sumMonth, northernHemWinter) # note dropping layers for northern hemisphere winter
chillHrsNorthernHem <- sum(chillHrsNorthernHem)
chillHrsSouthernHem <- sum(chillHrsSouthernHem)

fileNameNH <- paste0("chillHrs_NorthernHem", "_observed_", yearSpan, ".tif")
fileNameSH <- paste0("chillHrs_SouthernHem", "_observed_", yearSpan, ".tif")

writeRaster(chillHrsNorthernHem, filename = paste0("data/chillingHours/", fileNameNH),  overwrite = TRUE)
writeRaster(chillHrsSouthernHem, filename = paste0("data/chillingHours/", fileNameSH),  overwrite = TRUE)

gc(reset = FALSE, full = TRUE)
