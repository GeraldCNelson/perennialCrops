require(chillR)
require(sp)
require(tidyverse)
require(raster)
require(lubridate)
require(terra)
require(pbapply)
require(rasterVis)
require(future)
require(future.apply)
require(terra)
require(Rcpp)

options(future.globals.maxSize= 12000*1024^2)

source('R/chillSpatial_functions.R')

plan(multiprocess, workers=3, gc=TRUE)

# Rcpp.package.skeleton('CP', cpp_files = '../../chillPortions/global_chill/Cpp/chill_func.cpp')
# devtools::install('CP')
# require(CP)

# models <- tolower(c("GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL"))
# scenarios <- c("historical", "ssp126", "ssp585")

# an alternative approach is some nested for loops that replace all the individaul lines of code below
source("R/ISIMIPconstants.R") # loads a bunch of constants

# # chill portions, scenarios -----
# for (k in sspChoices) {
#   message(k)
#   for (l in startYearChoices) {
#     yearSpan <- l:as.numeric(paste0(l + yearRange))
#     message(paste0(' ', l))
#     for (modelChoice in modelChoices_lower) {
#     message(paste0('  ', modelChoice))
#     try(getChillWorld(scenario = k, model = modelChoice, year_range = yearSpan))
#     # cpName <- paste0(k, "_", modelChoice)
#     # assign(cpName, cp)
#     }
#   }
# }
# 
# # chill portions, historical -----
# k <- "historical"
# l <- 1991
#   for (l in startYearChoices) {
#     yearSpan <- l : as.numeric(paste0(l + yearRange))
#     for (modelChoice in modelChoices_lower) {
#       cp <- getChillWorld(scenario=k, model=modelChoice, year_range=yearSpan)
#       cpName <- paste0(k, "_", modelChoice)
#       assign(cpName, cp)
#     }
# }
#     
# t1 <- proc.time()
# cp <- getChillWorld(scenario=k, model=modelChoice, year_range=yearSpan)
# t2 <- t1-proc.time()
# t2

models <- tolower(c("GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL"))
scenarios <- c("historical", "ssp126", "ssp585")

plan(multiprocess, workers=3, gc=TRUE)
historical_GFDL_ESM4 <- getChillWorld(scenario=scenarios[1], model=models[1], year_range=1991:2010)
future:::ClusterRegistry("stop")
gc() 

plan(multiprocess, workers=3, gc=TRUE)
historical_IPSL_CM6A_LR <- getChillWorld(scenario=scenarios[1], model=models[2], year_range=1991:2010)
future:::ClusterRegistry("stop")
gc() 

plan(multiprocess, workers=3, gc=TRUE)
historical_MPI_ESM1_2_HR <- getChillWorld(scenario=scenarios[1], model=models[3], year_range=1991:2010)
future:::ClusterRegistry("stop")
gc() 

plan(multiprocess, workers=3, gc=TRUE)
historical_MRI_ESM2_0 <- getChillWorld(scenario=scenarios[1], model=models[4], year_range=1991:2010)
future:::ClusterRegistry("stop")
gc() 

plan(multiprocess, workers=3, gc=TRUE)
historical_UKESM1_0_LL <- getChillWorld(scenario=scenarios[1], model=models[5], year_range=1991:2010)
future:::ClusterRegistry("stop")
gc() 




plan(multiprocess, workers=3, gc=TRUE)
ssp126_GFCDL_ESM4 <- getChillWorld(scenario=scenarios[2], model=models[1], year_range=2041:2060) # should be GFDL
future:::ClusterRegistry("stop")
gc() 

plan(multiprocess, workers=3, gc=TRUE)
ssp126_IPSL_CM6A_LR <- getChillWorld(scenario=scenarios[2], model=models[2], year_range=2041:2060)
future:::ClusterRegistry("stop")
gc() 

plan(multiprocess, workers=3, gc=TRUE)
ssp126_MPI_ESM1_2_HR <- getChillWorld(scenario=scenarios[2], model=models[3], year_range=2041:2060)
future:::ClusterRegistry("stop")
gc() 

plan(multiprocess, workers=3, gc=TRUE)
ssp126_MRI_ESM2_0 <- getChillWorld(scenario=scenarios[2], model=models[4], year_range=2041:2060)
future:::ClusterRegistry("stop")
gc() 

plan(multiprocess, workers=3, gc=TRUE)
ssp126_UKESM1_0_LL <- getChillWorld(scenario=scenarios[2], model=models[5], year_range=2041:2060)
future:::ClusterRegistry("stop")
gc() 



plan(multiprocess, workers=3, gc=TRUE)
ssp585_GFDL_ESM4 <- getChillWorld(scenario=scenarios[3], model=models[1], year_range=2041:2060)
future:::ClusterRegistry("stop")
gc() 


########
plan(multiprocess, workers=3, gc=TRUE)
ssp585_IPSL_CM6A_LR <- getChillWorld(scenario=scenarios[3], model=models[2], year_range=2041:2060)
future:::ClusterRegistry("stop")
gc() 
#########

plan(multiprocess, workers=3, gc=TRUE)
ssp585_MPI_ESM1_2_HR <- getChillWorld(scenario=scenarios[3], model=models[3], year_range=2041:2060)
future:::ClusterRegistry("stop")
gc() 

plan(multiprocess, workers=3, gc=TRUE)
ssp585_MRI_ESM2_0 <- getChillWorld(scenario=scenarios[3], model=models[4], year_range=2041:2060)
future:::ClusterRegistry("stop")
gc() 

plan(multiprocess, workers=3, gc=TRUE)
ssp585_UKESM1_0_LL <- getChillWorld(scenario=scenarios[3], model=models[5], year_range=2041:2060)
future:::ClusterRegistry("stop")
gc() 





plan(multiprocess, workers=3, gc=TRUE)
ssp126_GFCDL_ESM4 <- getChillWorld(scenario=scenarios[2], model=models[1], year_range=2081:2100)
future:::ClusterRegistry("stop")
gc() 

plan(multiprocess, workers=3, gc=TRUE)
ssp126_IPSL_CM6A_LR <- getChillWorld(scenario=scenarios[2], model=models[2], year_range=2081:2100)
future:::ClusterRegistry("stop")
gc() 

plan(multiprocess, workers=3, gc=TRUE)
ssp126_MPI_ESM1_2_HR <- getChillWorld(scenario=scenarios[2], model=models[3], year_range=2081:2100)
future:::ClusterRegistry("stop")
gc() 

plan(multiprocess, workers=3, gc=TRUE)
ssp126_MRI_ESM2_0 <- getChillWorld(scenario=scenarios[2], model=models[4], year_range=2081:2100)
future:::ClusterRegistry("stop")
gc() 

plan(multiprocess, workers=3, gc=TRUE)
ssp126_UKESM1_0_LL <- getChillWorld(scenario=scenarios[2], model=models[5], year_range=2081:2100)
future:::ClusterRegistry("stop")
gc() 



plan(multiprocess, workers=3, gc=TRUE)
ssp585_GFDL_ESM4 <- getChillWorld(scenario=scenarios[3], model=models[1], year_range=2081:2100)
future:::ClusterRegistry("stop")
gc() 

plan(multiprocess, workers=3, gc=TRUE)
ssp585_IPSL_CM6A_LR <- getChillWorld(scenario=scenarios[3], model=models[2], year_range=2081:2100)
future:::ClusterRegistry("stop")
gc() 

plan(multiprocess, workers=3, gc=TRUE)
ssp585_MPI_ESM1_2_HR <- getChillWorld(scenario=scenarios[3], model=models[3], year_range=2081:2100)
future:::ClusterRegistry("stop")
gc() 

plan(multiprocess, workers=3, gc=TRUE)
ssp585_MRI_ESM2_0 <- getChillWorld(scenario=scenarios[3], model=models[4], year_range=2081:2100)
future:::ClusterRegistry("stop")
gc() 

plan(multiprocess, workers=3, gc=TRUE)
ssp585_UKESM1_0_LL <- getChillWorld(scenario=scenarios[3], model=models[5], year_range=2081:2100)
future:::ClusterRegistry("stop")
gc() 
