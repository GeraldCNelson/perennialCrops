library(terra)
JDay.south <- 92:306
test <- rast("climdata/gfdl-esm4_tasmax.south_historical_1991_2010.tif")
period <- 92:94
test[[period]]
testm <- as.matrix(test[[period]])
str(testm)
max(testm, na.rm = T)
tmin.tmp <- raster::brick(test[[period]]) 
tmin.tmp
tmin.tmp[] <- raster::as.matrix(tmin.tmp)
str(tmin.tmp[])
max(tmin.tmp[], na.rm = T)


r22 <- rast("data/chill_portions/historical/gfdl-esm4/historical_gfdl-esm4_2000_chill_portions_north.tif") # today

r21 <- rast("data/chillPortions/chill_portions/historical/gfdl-esm4/historical_gfdl-esm4_2000_chill_portions_north.tif") # 2021


library(geosphere)
library(terra)
r <- rast(res=0.5, ymin=-60, nlyr=365)
lat <- yFromRow(r, 1:nrow(r))
d <- sapply(1:365, function(i) daylength(lat, i))
values(r) <- rep(d, each=ncol(r))

DL <- function (latitude, JDay, notimes.as.na = FALSE) 
{  
  # if (missing(latitude)) 
  #   stop("'latitude' not specified")
  # if (missing(JDay)) 
  #   stop("'JDay' not specified")
  # if (!isTRUE(all(is.numeric(JDay)))) 
  #   stop("'JDay' contains non-numeric values")
  # if (length(latitude) > 1) 
  #   stop("'latitude' has more than one element")
  # if (!is.numeric(latitude)) 
  #   stop("'latitude' is not numeric")
  # if (latitude > 90 | latitude < (-90)) 
  #   warning("'latitude' is usually between -90 and 90")
  
  # add days to JDay to account for omissions on either side
  # JDay <- c(JDay[1]-1, JDay, JDay[length(JDay)]+1)
  print("Check 6.1"); print(Sys.time())
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
  
  # termA <- rast(stack(lapply(seq_along(Delta.sin), function(x) Delta.sin[x] * CosWo.a[[x]])))
  # termB <- rast(stack(lapply(seq_along(Delta.cos), function(x) Delta.cos[x] * CosWo.b[[x]])))
  # termA <- lapply(seq_along(Delta.sin), function(x) Delta.sin[x] * CosWo.a[[x]]) %>% do.call('c', .)
  # termB <- lapply(seq_along(Delta.cos), function(x) Delta.cos[x] * CosWo.b[[x]]) %>% do.call('c', .)
  
  termA <- stack(lapply(seq_along(Delta.sin), function(x) Delta.sin[x] * CosWo.a[[x]]))
  termB <- stack(lapply(seq_along(Delta.cos), function(x) Delta.cos[x] * CosWo.b[[x]]))
  
  CosWo <- termA/termB
  print("Check 6.2"); print(Sys.time())
  # CosWo <- (CosWo.1 - sin(latitude/360 * 2 * pi) * sin(Delta/360 * 2 * pi))/(cos(latitude/360 * 
  # 2 * pi) * cos(Delta/360 * 2 * pi))
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


