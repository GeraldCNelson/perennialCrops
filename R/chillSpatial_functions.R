# calculate day lengths from latitudes
# derived from chillR package

DL <- function (latitude, JDay, notimes.as.na = FALSE) 
{  
  
  Gamma <- 2 * pi/365 * ((JDay) - 1)
  Delta <- 180/pi * (0.006918 - 0.399912 * cos(Gamma) + 0.070257 * 
                       sin(Gamma) - 0.006758 * cos(Gamma) + 0.000907 * sin(Gamma) - 
                       0.002697 * cos(3 * (Gamma)) + 0.00148 * sin(3 * (Gamma)))
  
  CosWo.1 <- latitude
  values(CosWo.1) <- sin(-0.8333/360 * 2 * pi)
  
  Delta.sin <- sin(Delta/360 * 2 * pi)
  Delta.cos <- cos(Delta/360 * 2 * pi)
  
  CosWo.a <- (CosWo.1 - sin(latitude/360 * 2 * pi))
  CosWo.b <- (cos(latitude/360 * 2 * pi))
  
  termA <- stack(lapply(seq_along(Delta.sin), function(x) Delta.sin[x] * CosWo.a[[x]]))
  termB <- stack(lapply(seq_along(Delta.cos), function(x) Delta.cos[x] * CosWo.b[[x]]))
  
  CosWo <- termA/termB

  normal_days <- as.vector(CosWo[] >= -1 & CosWo[] <= 1)
  
  Sunrise <- rep(-99, length(CosWo[]))
  Sunrise[normal_days[]] <- 12 - acos(CosWo[][normal_days])/(15/360 *
                                                               2 * pi)
  Sunset <- rep(-99, length(CosWo[]))
  Sunset[normal_days[]] <- 12 + acos(CosWo[][normal_days])/(15/360 * 
                                                              2 * pi)
  Daylength <- Sunset - Sunrise
  Daylength[which(CosWo[] > 1)] <- 0
  Daylength[which(CosWo[] < (-1))] <- 24
  Sunrise[which(Daylength == 24)] <- 99
  Sunset[which(Daylength == 24)] <- 99
  if (notimes.as.na) {
    Sunrise[which(Sunrise %in% c(-99, 99))] <- NA
    Sunset[which(Sunset %in% c(-99, 99))] <- NA
  }
  Sunset[which(is.na(JDay))] <- NA
  Sunrise[which(is.na(JDay))] <- NA
  Daylength[which(is.na(JDay))] <- NA
  
  Sunrise[which(Sunrise == 99)] <- 0
  Sunrise[which(Sunrise == -99)] <- 12
  Sunset[which(Sunset == 99)] <- 24
  Sunset[which(Sunset == -99)] <- 12
  
  times=rep(JDay, each=length(Sunrise)/length(JDay))
  # cell <- rep(1:ncell(lat), times=length(unique(JDay)))
  cell <- rep(1:ncell(latitude), times=length(unique(JDay)))
  
  Sunrise <- round(Sunrise,2)
  Sunset <- round(Sunset,2)
  Daylength <- round(Daylength,2)
  
  Sunrise <- split(Sunrise, f=cell, drop=TRUE)
  Sunset <- split(Sunset, f=cell, drop=TRUE)
  Daylength <- split(Daylength, f=cell, drop=TRUE)
  
  return(list(Sunrise = Sunrise, Sunset = Sunset, Daylength = Daylength, JDay=rep(JDay, each=length(Sunrise))))#/length(JDay))))
}

# calculate chill portions from daily tmin and tmax
getCP <- function (tmin=tmin, tmax=tmax, dates, Day_times=daytimes, keep_sunrise_sunset = FALSE, template=NULL, datClass, ...) 
{
  
  datcells <- which(!is.na(template[]))
  
  # message('  Loading temperature data into memory..')
  
  tmax <- as.matrix(tmax)
  tmin <- as.matrix(tmin)
  
  # message('  Processing daylengths..')
  
  dates.all <- vector('list', length(template[]))
  
  # pb <- txtProgressBar(max=length(datcells), style=3)
  for(n in datcells) {
    
    if(datClass=='list') {
      
      dates.cell <- data.frame(Year = dates$Year,
                               JDay = dates$JDay,
                               Tmax = tmax[n,][dates$JDay],
                               Tmin = tmin[n,][dates$JDay],
                               cell = n,
                               Sunrise = NA,
                               Sunset = NA,
                               Daylength = NA,
                               prev_Sunset = NA,
                               next_Sunrise = NA,
                               prev_max = NA,
                               next_min = NA,
                               prev_min = NA)    
    } else if (datClass=='stack') {
      
      dates.cell <- data.frame(Year = dates$Year,
                               JDay = dates$JDay,
                               Tmax = tmax[n,], # tmax is pre-cropped to JDay range
                               Tmin = tmin[n,], # tmin is pre-cropped to JDay range
                               cell = n,
                               Sunrise = NA,
                               Sunset = NA,
                               Daylength = NA,
                               prev_Sunset = NA,
                               next_Sunrise = NA,
                               prev_max = NA,
                               next_min = NA,
                               prev_min = NA)            
      any(is.na(dates.cell$Tmax))
      
    } else {
      stop('Invalid datclass')
    }
    
    
    Day_times.cell <- list(Sunrise=Day_times$Sunrise[[n]], Sunset=Day_times$Sunset[[n]], Daylength=Day_times$Daylength[[n]])
    
    dates.all[[n]] <- MHT(dates.cell= dates.cell, 
                          Day_times.cell = list(Sunrise=Day_times$Sunrise[[n]], Sunset=Day_times$Sunset[[n]], Daylength=Day_times$Daylength[[n]]))
    
    # message('  Calculating chill portions..')

    getChill <- function(x) {
      if(!is.null(x)) {
        chill_func(na.omit(as.vector(t(x))))
      } else {
        NA
      }
    }
    
  }

  CP <- lapply(dates.all, getChill) # 6:29 refers to columns the hourly temperature data
    
  template[] <- NA
  template[] <- unlist(CP)
    
  return(template)
  
}

# parent function for calculating chill portions
getChillSpatial <- function(years, lat, JDay, tmin, tmax, template, writeToDisk=FALSE,...) {
  message('Annualizing temperature data..')
  
  # add days to JDay to account for omissions on either side
  JDay <- c(JDay[1]-1, JDay, JDay[length(JDay)]+1)
  
  if(class(tmin)=='list') {
    
    datclass='list'
    
    # interpolate hourly temperatures  
    daytimes <-  DL(lat, JDay)
    CP <- future_lapply(seq_along(years), function(x) getCP(tmin=tmin[[x]],
                                                            tmax=tmax[[x]],
                                                            Day_times=daytimes,
                                                            template=tmin[[1]][[1]],
                                                            dates = data.frame(Year = years[x],
                                                                               JDay = JDay),
                                                            datClass=datclass))
    
  } else {
    
    datclass='stack'
    
    if(class(tmin[[1]]) %in% c("RasterBrick", "RasterStack", "RasterLayer")) {
      dates <- as.Date(1:nlayers(tmin), origin=as.Date(paste0(years[1], "-01-01"))-1)
    } else {
      if(class(tmin[[1]]) %in% c("SpatRaster")) {
        dates <- as.Date(1:nlyr(tmin), origin=as.Date(paste0(years[1], "-01-01"))-1)
      }
    }
    
    dates.year <- lubridate::year(dates)
    yearRange <- years
    
    
  # loop through years, finding layers associated with each year's worth of data
  # assign a year's worth of data to each element of a list
    
    tmin.out <- vector("list", length(unique(dates.year)))
    tmax.out <- vector("list", length(unique(dates.year)))
    period.dates <- vector("list", length(unique(dates.year)))
    
    # create template raster stack to assign daily values to
    ext_hem <- ext(tmin)
    
    pb <- txtProgressBar(max=length(yearRange), style=3)
    for(i in seq_along(yearRange)) {
      setTxtProgressBar(pb, i)
      firstDay <- which(dates.year %in% years[i])[JDay][1]
      lastDay <- firstDay + length(JDay)-1
      period <- firstDay:lastDay
      period.dates[[i]] <- dates[period]
      
      brick_empty <- brick(nrows=dim(tmin)[1], ncols=dim(tmin)[2], xmn=-180, xmx=180, ymn=ext_hem[3], ymx=ext_hem[4], nl=length(period), crs = "+proj=longlat +datum=WGS84 +no_defs" )
      
      tmin.tmp <- brick_empty
      m_tmin <- as.matrix(tmin[[period]])
      tmin.tmp[] <- m_tmin
      
      tmax.tmp <- brick_empty
      m_tmax <- as.matrix(tmax[[period]])
      tmax.tmp[] <- m_tmax

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

getChillWorld <- function(scenario, model, year_range) {
  
  # get data files from scenario and model arguments
  dat.files <- list.files(dat.dir, pattern=model, recursive=TRUE, full.names = TRUE)
  dat.files <- grep('xml', dat.files, invert = TRUE, value = TRUE)
  tmin.files <- grep("tasmin", dat.files, value=TRUE) %>% grep("aux|xml|json", ., value=TRUE, invert=TRUE)
  tmax.files <- grep("tasmax", dat.files, value=TRUE) %>% grep("aux|xml|json", ., value=TRUE, invert=TRUE)
  
  if(all(year_range==1991:2010)) {
    tmin.in <- rast(tmin.files[1])
    tmax.in <- rast(tmax.files[1])
  } else if (all(year_range==2041:2060)) {
    tmin.in <- rast(tmin.files[1])
    tmax.in <- rast(tmax.files[1])
  } else if(all(year_range==2081:2100)) {
    tmin.in <- rast(tmin.files[2])
    tmax.in <- rast(tmax.files[2])
  } else {
    stop("Unsupported year range specified.")
  }
  
  # test files
  
  testVal_min <- min(minmax(tmin.in))
  testVal_max <- max(minmax(tmin.in))
  if(testVal_min < -80 | testVal_max > 100) {
    stop(paste0("min: ", testVal_min, ", max: ", testVal_max, ", file: ", scenario, " ", model, " tmin ", year_range[1]))
  }
  
  testVal_min <- min(minmax(tmax.in))
  testVal_max <- max(minmax(tmax.in))
  if(testVal_min < -75 | testVal_max > 100) {
    stop(paste0("min: ", testVal_min, ", max: ", testVal_max, ", file: ", scenario, " ", model, " tmin ", year_range[1]))
  }
  
  # set dormancy period as a vector of julian days (for northern hemisphere dormancy, days in second calendar year should be represented as 365+n)
  
  JDay.south <- 92:306
  JDay.north <- 92:306 + 180
  
  # chill portions spatial calculations
  
  #### run calculations on northern hemisphere ####
  message('  Cropping extent to northern hemisphere..')
  tmin.north <- crop(tmin.in, ext(c(-180,180,0,90)))
  # tmin.north <- crop(tmin.in, ext(c(-100,-80,30,50)))
  gc()
  tmax.north <- crop(tmax.in, ext(c(-180,180,0,90)))
  # tmax.north <- crop(tmax.in, ext(c(-100,-80,30,50)))
  gc()

  # produce spatial layer of latitudes and use to calculate day lengths
  lat <- tmin.north[[c(JDay.north[1]-1, JDay.north, JDay.north[length(JDay.north)]+1)]]
  xy <- coordinates(raster(tmin.north[[1]][[1]]))
  values(lat) <- xy[,2]

  rm(tmin.in, tmax.in)
  gc()

  # set a template raster to define the attributes of the final output raster
  # this can be a single layer of temperature data, for example
  template.ras <- tmin.north[[1]]

  # process daily temperatures
  year_range.north <- year_range[1:(length(year_range)-1)]
  t <- proc.time()
  chill_portions.north <- getChillSpatial(years=year_range.north, lat, JDay.north, tmin=tmin.north, tmax=tmax.north, template=template.ras)
  t1 <- t-proc.time()
  t1/60

  # dir.create(paste0('data/chill_portions/', scenario, '/', model, '/', year_range[10]), recursive = TRUE)
  
  outD <- paste0('data/chill_portions/', scenario, '/', model, '/', year_range[10], "/")
  if (!dir.exists(outD)) dir.create(outD, recursive = TRUE)

  chill_portions.north.raster <- rast(stack(chill_portions.north))

  writeRaster(chill_portions.north.raster,
              paste0('data/chill_portions/',
                     scenario, '/', model, '/', scenario, '_', model, '_', year_range[10], '_', 'chill_portions_north.tif'),
              overwrite=TRUE)

  rm(tmin.north, tmax.north)
  gc()
  
  #### run calculations on southern hemisphere ####
  
  if(all(year_range==1991:2010)) {
    tmin.in <- rast(tmin.files[1])
    tmax.in <- rast(tmax.files[1])
  } else if (all(year_range==2041:2060)) {
    tmin.in <- rast(tmin.files[1])
    tmax.in <- rast(tmax.files[1])
  } else if(all(year_range==2081:2100)) {
    tmin.in <- rast(tmin.files[2])
    tmax.in <- rast(tmax.files[2])
  } else {
    stop("Unsupported year range specified.")
  }
  
  future:::ClusterRegistry("stop")
  gc()
  plan(multisession, workers=3, gc=TRUE)
  
  message('  Cropping extent to southern hemisphere (excluding Antarctica)..')
  tmin.south <- crop(tmin.in, ext(c(-180,180,-60,0))) # -60 no antarctica
  gc()
  tmax.south <- crop(tmax.in, ext(c(-180,180,-60,0)))
  
  rm(tmin.in, tmax.in)
  gc()
  
  # produce spatial layer of latitudes and use to calculate day lengths
  lat <- tmin.south[[c(JDay.south[1]-1, JDay.south, JDay.south[length(JDay.south)]+1)]]
  
  xy <- coordinates(raster(tmin.south[[1]][[1]]))
  values(lat) <- xy[,2]
  # set a template raster to define the attributes of the final output raster
  # this can be a single layer of temperature data, for example
  template.ras <- tmin.south[[1]]#[[1]]
  
  # process daily temperatures
  t <- proc.time()
  chill_portions.south <- getChillSpatial(years=year_range, lat, JDay.south, tmin=tmin.south, tmax=tmax.south, template=template.ras)
  t1 <- t-proc.time()
  t1/60
  
  chill_portions.south.raster <- rast(stack(chill_portions.south))
  writeRaster(chill_portions.south.raster,
              paste0('data/chill_portions/',
                     scenario, '/', model, '/', scenario, '_', model, '_', year_range[10], '_', 'chill_portions_south.tif'),
              overwrite=TRUE)
  rm(tmin.south, tmax.south)
  gc()
  
  # # stitch northern and southern hemisphere rasters together
  # 
  # world <- lapply(seq_along(1:length(chill_portions.north)), function(x) {
  #   north.x <- raster(chill_portions.north[[x]])
  #   south.x <- raster(chill_portions.south[[x]])
  #   world.x <- mosaic(north.x, south.x, fun=min)
  #   return(world.x)
  # })
  # 
  # writeRaster(rast(stack(world)), 
  #             paste0('data/chill_portions/', 
  #                    scenario, '/', model, '/', scenario, '_', model, '_', year_range[10], '_', 'chill_portions_world.tif'),
  #             overwrite=TRUE)
  # 
  # return(world)
  
}



# interpolate daily tmin and tmax to hourly temperatures
# code modified from chillR::make_hourly_temps()

MHT <- function (dates.cell, Day_times.cell, keep_sunrise_sunset = FALSE) 
{
  if(!is.null(dates.cell)) {
    year_file <- dates.cell
    year_file <- year_file[which(!is.na(year_file$Tmin) & !is.na(year_file$Tmax)),]
    year_file <- year_file[2:(nrow(year_file)-1),]

    preserve_columns <- colnames(year_file)
    Day_times <- Day_times.cell
    
    Day_times$Sunrise[which(Day_times$Sunrise == 99)] <- 0
    Day_times$Sunrise[which(Day_times$Sunrise == -99)] <- 12
    Day_times$Sunset[which(Day_times$Sunset == 99)] <- 24
    Day_times$Sunset[which(Day_times$Sunset == -99)] <- 12
    
    year_file$Sunrise <- Day_times$Sunrise[2:(length(Day_times$Sunrise) - 1)]
    year_file$Sunset <- Day_times$Sunset[2:(length(Day_times$Sunset) - 1)]
    year_file$Daylength <- Day_times$Daylength[2:(length(Day_times$Daylength) -  1)]
    year_file$prev_Sunset <- Day_times$Sunset[1:(length(Day_times$Sunset) -  2)]
    year_file$next_Sunrise <- Day_times$Sunrise[3:length(Day_times$Sunrise)]
    year_file$prev_max <- year_file$Tmax[c(NA, 1:(nrow(year_file) -   1))]
    year_file$next_min <- year_file$Tmin[c(2:nrow(year_file),  NA)]
    year_file$prev_min <- year_file$Tmin[c(NA, 1:(nrow(year_file) -  1))]
    year_file$Tsunset <- year_file$Tmin + (year_file$Tmax - year_file$Tmin) * 
      sin((pi * (year_file$Sunset - year_file$Sunrise)/(year_file$Daylength +  4)))
    year_file$prev_Tsunset <- year_file$prev_min + 
      (year_file$prev_max - year_file$prev_min) * sin((pi * (year_file$Daylength)/(year_file$Daylength +   4)))
    colnum <- ncol(year_file) + 1
    hourcol <- c(colnum:(colnum + 23))
    for (hour in 0:23) {
      hourcount <- hour + 1
      no_riseset <- which(year_file$Daylength %in% c(0, 24, -99))
      year_file[no_riseset, colnum + hour] <- ((year_file$Tmax + 
                                                  year_file$Tmin)/2)[no_riseset]
      c_morn <- which(hour <= year_file$Sunrise)
      if (1 %in% c_morn) 
        if (!length(c_morn) == 1) 
          c_morn <- c_morn[2:length(c_morn)]
      else c_morn <- c()
      c_day <- which(hour > year_file$Sunrise & hour <= year_file$Sunset)
      c_eve <- which(hour >= year_file$Sunset)
      if (nrow(year_file) %in% c_eve) 
        c_eve <- c_eve[1:(length(c_eve) - 1)]
      year_file[c_morn, colnum + hour] <- year_file$prev_Tsunset[c_morn] - 
        ((year_file$prev_Tsunset[c_morn] - year_file$Tmin[c_morn])/
           log(max(1,24 - (year_file$prev_Sunset[c_morn] - year_file$Sunrise[c_morn]))) * 
           log(hour + 24 - year_file$prev_Sunset[c_morn] +  1))
      year_file[c_day, colnum + hour] <- year_file$Tmin[c_day] + 
        (year_file$Tmax[c_day] - year_file$Tmin[c_day]) * 
        sin((pi * (hour - year_file$Sunrise[c_day])/(year_file$Daylength[c_day] + 
                                                       4)))
      year_file[c_eve, colnum + hour] <- year_file$Tsunset[c_eve] - 
        ((year_file$Tsunset[c_eve] - year_file$next_min[c_eve])/
           log(24 - (year_file$Sunset[c_eve] - year_file$next_Sunrise[c_eve]) + 
                 1) * log(hour - year_file$Sunset[c_eve] + 1))
    }
    
    colnames(year_file)[(ncol(year_file) - 23):(ncol(year_file))] <- c(paste("Hour_", 0:23, sep = ""))
    if (!keep_sunrise_sunset) 
      year_file <- year_file[, c(preserve_columns, paste("Hour_", 0:23, sep = ""))]
    if (keep_sunrise_sunset) 
      year_file <- year_file[, c(preserve_columns, "Sunrise", 
                                 "Sunset", "Daylength", paste("Hour_",  0:23, sep = ""))]
    year_file[1, (ncol(year_file) - 23):(ncol(year_file))][
      which(is.na(year_file[1, 
                            (ncol(year_file) - 23):(ncol(year_file))]))] <- year_file[1,"Tmin"]
    year_file[nrow(year_file), (ncol(year_file) - 23):(ncol(year_file))][
      which(is.na(year_file[nrow(year_file), 
                            (ncol(year_file) - 23):(ncol(year_file))]))] <- year_file[nrow(year_file), "Tmin"]
    matrix.out <- as.matrix(year_file[,14:37])
    return(matrix.out)
  }
}

