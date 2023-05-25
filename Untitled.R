library(terra)
library(raster)
model <- tolower(c("GFDL-ESM4"))
scenario <- "historical"
JDay <- JDay.south <- 92:306
dat.dir <- paste0('climdata')
#  dat.files <- list.files(dat.dir, pattern=model, recursive=TRUE, full.names = TRUE)
dat.files <- list.files(dat.dir, pattern=model, recursive=TRUE, full.names = TRUE)
dat.files <- grep('xml', dat.files, invert = TRUE, value = TRUE)
# dat.files <- grep(year_range, dat.files, invert = TRUE, value = TRUE)
tmin.files <- grep("tasmin", dat.files, value=TRUE)
ext_SH <- ext(c(-180, 180, -60, 0))
tmin.south <- rast(tmin.files[1], win = ext_SH) 
dates <- as.Date(1:nlyr(tmin.south), origin=as.Date(paste0(years[1], "-01-01"))-1)
lat <- tmin.south[[c(JDay.south[1]-1, JDay.south, JDay.south[length(JDay.south)]+1)]]
xy <- coordinates(raster(tmin.south[[1]][[1]]))
values(lat) <- xy[,2]
daytimes <-  DL(lat, JDay)
template.ras <- tmin.south[[1]]
year_range=1991:2010
years <- year_range.south <- year_range

i = 1

chill_portions.south <- getChillSpatial(years=year_range.south, lat, JDay.south, tmin=tmin.south, tmax=tmax.south, template=template.ras)
tmin=tmin.south; tmax=tmax.south; template=template.ras

dates.year <- lubridate::year(dates)

firstDay <- which(dates.year %in% years[i])[JDay][1]
lastDay <- firstDay + length(JDay)-1
period <- firstDay:lastDay
period.dates <- vector("list", length(unique(dates.year)))

period.dates[[i]] <- dates[period]
# print("starting brick NH")
# print(Sys.time())
#  browser()
t <- tmin[[period]]
 tmin.tmp <- brick(t)
 t_1 <- tmin.tmp[[1]]
 t_1[] <- as.matrix(t_1)
tmin.tmp <- as.matrix(tmin[[period]]) # note that [] has been removed

getChillSpatial <- function(years, lat, JDay, tmin, tmax, template, writeToDisk=FALSE,...) {
  message('Annualizing temperature data..')
  print(Sys.time())
  browser()
  # add days to JDay to account for omissions on either side
  JDay <- c(JDay[1]-1, JDay, JDay[length(JDay)]+1)
  # print("Done with JDays")
  # print(Sys.time())
  # print(class(tmin))
  
  if(class(tmin)=='list') {
    
    datclass='list'
    
    # interpolate hourly temperatures  
    print("starting DL calcs")
    print(Sys.time())
    
    daytimes <-  DL(lat, JDay)
    print("end DL calcs")
    print(Sys.time())
    CP <- future_lapply(seq_along(years), function(x) 
      getCP(tmin=tmin[[x]],
            tmax=tmax[[x]],
            Day_times=daytimes,
            template=tmin[[1]][[1]],
            dates = data.frame(Year = years[x], JDay = JDay), datClass=datclass))
  } else {
    print("starting DL calcs with stack")
    print(Sys.time())
    
    datclass='stack'
    
    if(class(tmin[[1]]) %in% c("RasterBrick", "RasterStack", "RasterLayer")) {
      dates <- as.Date(1:nlayers(tmin), origin=as.Date(paste0(years[1], "-01-01"))-1)
    } else {
      if(class(tmin[[1]]) %in% c("SpatRaster")) {
        dates <- as.Date(1:nlyr(tmin), origin=as.Date(paste0(years[1], "-01-01"))-1)
      }
    }
    
    dates.year <- lubridate::year(dates)
    print("vector tmin and tmax out")
    print(Sys.time())
    
    tmin.out <- vector("list", length(unique(dates.year)))
    tmax.out <- vector("list", length(unique(dates.year)))
    period.dates <- vector("list", length(unique(dates.year)))
    
    yearRange <- years
    
    # extract layers corresponding to each year in yearRange and assign each year's worth of data to list elements
    
    pb <- txtProgressBar(max=length(yearRange), style=3)
    for (i in seq_along(yearRange)) {
      print(paste0("year: ", i))
      setTxtProgressBar(pb, i)
      firstDay <- which(dates.year %in% years[i])[JDay][1]
      lastDay <- firstDay + length(JDay)-1
      period <- firstDay:lastDay
      period.dates[[i]] <- dates[period]
      # print("starting brick NH")
      # print(Sys.time())
      #  browser()
      # tmin.tmp <- brick(tmin[[period]])
      # tmin.tmp[] <- as.matrix(tmin.tmp)
      tmin.tmp <- as.matrix(tmin[[period]]) # note that [] has been removed
      # tmax.tmp <- brick(tmax[[period]])
      # tmax.tmp[] <- as.matrix(tmax.tmp)
      tmax.tmp <- as.matrix(tmax[[period]]) # note that [] has been removed
      print("end brick tmax NH")
      print(Sys.time())
      
      # tmin.tmp <- tmin[[period]]
      # tmin.tmp[] <- as.matrix(tmin.tmp)
      # 
      # tmax.tmp <- tmax[[period]]
      # tmax.tmp[] <- as.matrix(tmax.tmp)
      
      tmin.out[[i]] <- tmin.tmp
      tmax.out[[i]] <- tmax.tmp
    }
    close(pb)
    
    # run getCP function across years range
    
    if(max(JDay)>365) {
      JDay[JDay>365] <- JDay[JDay>365]-365
    }
    
    message('  Calculating daylengths..')
    daytimes <-  DL(brick(lat), JDay)
    
    message('  Calculating chill portions..')
    tmp <- raster(tmin[[1]][[1]])
    dates.list <- lapply(seq_along(years), function(x) data.frame(Year = year(period.dates[[x]]),
                                                                  JDay = yday(period.dates[[x]])))
    rm(tmin, tmax)
    gc()
    
    CP <- future_lapply(seq_along(years), function(x) getCP(tmin=tmin.out[[x]],
                                                            tmax=tmax.out[[x]],
                                                            Day_times=daytimes,
                                                            template=tmp,
                                                            dates = dates.list[[x]],
                                                            datClass=datclass))
    
  }
  
  return(CP)
  
}

