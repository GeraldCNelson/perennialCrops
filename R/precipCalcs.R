# analysis of precipitation data

source("R/globallyUsed.R")
library(doParallel) #Foreach Parallel Adaptor 
library(foreach) #Provides foreach looping construct
library(stringr)

sspChoices <- c("ssp585") #"ssp126", "ssp585"
modelChoices <- c("GFDL-ESM4", "MRI-ESM2-0", "MPI-ESM1-2-HR", "UKESM1-0-LL",  "IPSL-CM6A-LR") #, "MPI-ESM1-2-HR", "MRI-ESM2-0") # "GFDL-ESM4", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL", "IPSL-CM6A-LR"
#modelChoices <- c( "MRI-ESM2-0")
climateVars <- c( "pr") # "tasmax", "pr", "hurs") # "tasmin", tasmax
startYearChoices <-  c(2021, 2051, 2091) #2011, 2041, 2051, 2081) # c(2091) # c(2006) #, 2041, 2051, 2081)
locOfFiles <- locOfCMIP6ncFiles
yearRange <- 9

# commented out, now in the globallyUsed.R script
#cropCharacteristics_annual <- as.data.table(read_excel("data-raw/crops/cropCharacteristics_annual_summary_02052020.xlsx", range = "A1:S26"))
#setnames(cropCharacteristics_annual, old = names(cropCharacteristics_annual), new = make.names(names(cropCharacteristics_annual)))

cropChoices <- unique(cropCharacteristics_annual$crop)

# info on managing raster disk use - https://stackoverflow.com/questions/25426405/raster-package-taking-all-hard-drive

#Sys.setenv(PROJ_LIB = "/usr/local/Cellar/proj/6.2.1/share/proj") # use until the sf and gdal issues get sorted out. See https://github.com/r-spatial/sf/issues/1060

useCores <- detectCores() - 1 # max number of cores
useCores <- 3 # better for memory intensive activities

varList <- c("startYearChoices", "sspChoices", "modelChoices", "wrld_land",  "locOfFiles", "cropChoices")
libList <- c("raster", "ncdf4", "stringr")

#test values
i <- "IPSL-CM6A-LR"
k <- "ssp585"
l <- 2091
j = "pr"
m = 1

cl <- clusterSetup(varList, libList, useCores) # function created in globallyUsed.R

foreach(i = modelChoices) %:%
  foreach(l = startYearChoices) %:%
  foreach(j = climateVars) %:%
  foreach(k = sspChoices) %dopar% {
    require(raster)
    tmpDirName <- paste0(locOfFiles, "rasterTmp_", Sys.getpid(), "/")
    rasterOptions(tmpdir = tmpDirName)
    dir.create(tmpDirName)
    
    print(paste0("working on start year: ", l, ", variable: ", j, ", ssp choice: ", k, ", model: ", i,  ", pid: ", Sys.getpid(), ", systime: ", Sys.time()))
    modelName.lower <- tolower(i)
    startTime <-  Sys.time()
    yearSpan <- paste0(l, "_", l + yearRange)
    
    fileName_in <- paste(modelName.lower, k, j, "global_daily", yearSpan, sep = "_")
    fileName_in <- paste0(fileName_in, ".nc")
    
    temp <- paste(locOfFiles, k, "/", i, "/", fileName_in, sep = "")
    print(paste0("Working on : ", temp))
    brick.pr <- rasttemp, varname = j) # because there is no explicit projection info in the netcdf files, this is assumed - +proj=longlat +datum=WGS84"

    # the overlay function needs a user defined function on the relationship between the two rasters
    overlayfunction <- function(x,y) {
      return(x * y)
    }
    
    for (m in 1:length(cropChoices)) { #fruits variable is from the globallyUsed.R file
      cropName <- cropChoices[m]
      fileNameMask.in <- paste0("data/crops/rasterMask_", tolower(cropName), ".tif")
      print(paste0("fileNameMaskIn: ", fileNameMask.in))
      mask <- rast(fileNameMask.in)
      startTime <- Sys.time()
      brick.pr.masked <- overlay(brick.pr, mask, fun = overlayfunction)
      endTime <- Sys.time()
      temp <- endTime - startTime
      print(paste0("time to crop precip data: ", temp))
      
      
       fileNameMean.masked <- paste0("data/cmip6/chillingHours/chillHrs_", hemisphereName, "_ensembleMean_masked_", cropName, "_",  yearSpan, "_", k, ".tif")
      fileNameCV.masked <- paste0("data/cmip6/chillingHours/chillHrs_", hemisphereName, "_ensembleCV_masked_", cropName, "_",  yearSpan, "_", k, ".tif")
      print(paste("fileNameMean.masked: ", fileNameMean.masked))
      writeRaster(mean.masked, filename = fileNameMean.masked,  overwrite = TRUE)
      writeRaster(CV.masked, filename = fileNameCV.masked,  overwrite = TRUE)
    }
    
    
    
    
    
    fileNameMask.in <- paste0("data/crops/rasterMask_", tolower(speciesName), ".tif")
    
    cropName <- m
    
    fileNameMask.in <- paste0("data/crops/rasterMask_", tolower(m), ".tif")
    cropMask <- rast(fileNameMask.in)
    
    cropCalendarName <- cropCharacteristics_annual[crop %in% cropName, crop.calendar]
    cropCalFilesLoc <- paste0("data-raw/crops/cropCalendars/ALL_CROPS_netCDF_0.5deg_filled/")
    fileName_in <- paste0(cropCalendarName, ".crop.calendar.fill.nc")
    #    locNFileIn <- paste0(filesLoc, fileName_in, ".gz")
    locNFileIn <- paste0(cropCalFilesLoc, fileName_in)
    R.utils::gunzip(paste0(locNFileIn, ".gz"), remove = FALSE)
    croppingCalendar <- rast(locNFileIn)
    crs(croppingCalendar) <- crs(cropMask) # needed because cropping calendar doesn't have an explicit crs
    croppingCalendar_mask <- mask(croppingCalendar, cropMask)
    croppingCalendar_plant <- croppingCalendar_mask$plant
    croppingCalendar_harvest <- croppingCalendar_mask$harvest
    cal <- c(croppingCalendar_plant_crop, croppingCalendar_plant_crop)
    unlink(locNFileIn) # delete the .nc file when no longer needed.
