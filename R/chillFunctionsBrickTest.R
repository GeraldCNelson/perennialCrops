library(terra)
test <- rast("climdata/gfdl-esm4_tasmax.south_historical_1991_2010.tif")
period <- 92:94
test[[period]]
ptm <- proc.time()
testm[] <- as.matrix(test[[period]])
proc.time() - ptm
str(testm)
max(testm, na.rm = T)
ptm <- proc.time()
tmin.tmp <- raster::brick(test[[period]]) 
#tmin.tmp
tmin.tmp[] <- as.matrix(tmin.tmp)
proc.time() - ptm
str(tmin.tmp[])
max(tmin.tmp[], na.rm = T)

# test to see if they are the same

matequal <- function(x, y) {
  is.matrix(x) && is.matrix(y) && dim(x) == dim(y) && all(x == y)
  print(is.matrix(x))
  print(is.matrix(y))
  print(dim(x))
  print(dim(y))
  
}
matequal(testm, tmin.tmp[])

all.equal(testm, tmin.tmp[])

r22 <- rast("data/chill_portions/historical/gfdl-esm4/historical_gfdl-esm4_2000_chill_portions_north.tif") # today

r21 <- rast("data/chillPortions/chill_portions/historical/gfdl-esm4/historical_gfdl-esm4_2000_chill_portions_north.tif") # 2021

all(testm == tmin.tmp)


library(geosphere)
library(terra)
r <- rast(res=0.5, ymin=-60, nlyr=365)
lat <- yFromRow(r, 1:nrow(r))
d <- sapply(1:365, function(i) daylength(lat, i))
values(r) <- rep(d, each=ncol(r))

