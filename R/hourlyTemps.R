library(terra)
ext_SH <- ext(-180, 180, -60, 0)
ext_NH <- ext(-180, 180, 0, 90)
ext_globe <- ext(-180, 180, -60, 90)
JDay.SH <- (92-2):(306)
JDay.NH <- JDay.SH + 180
JDay.globe <- 1:365
path <- paste0('climdata')

# test data
year <- 1991
hem <- "SH"; jDayHem <- JDay.SH# test data
year_range <- 1991:2010
model <- tolower("GFDL-ESM4")

#hem choices are globe, NH, SH
create_lats <- function(hem, jDayHem, year) { # generate spatraster with latitude values and number of layers between j days
  jdays <- tail(jDayHem,1) - head(jDayHem,1) + 1 # not checked for validity in all cases
  
  if (hem == "NH") e <- ext_NH
  if (hem == "SH") e <- ext_SH
  if (hem == "globe") {
    e <- ext(-180, 180, -60, 90)
    if (lubridate::leap_year(year)) jdays <- 366
    jDayHem <- 1:jdays
  }
  lat_dates <- as.Date(jDayHem, origin = as.Date(paste0(year,"-01-01"))) 
  lat <- rast(extent = e) |> init("y") |> rep(jdays)
  names(lat) <- paste0("X", lat_dates)
  time(lat) <- lat_dates
  
  dl <- meteor::photoperiod(lat)
  sunrise <- 12 - dl/2
  sunset <- 12 + dl/2
  latsInfo <- sds(lat, dl, sunrise, sunset)
  return(latsInfo)
}

f_comb <- function(tmax, tmin) {
  td <- (tmax - tmin) / 2
  ts <- (tmax + tmax) / 2
  tout <- sds(td, ts)
}

f_hourly_tempWdeltas <- function(hod, td, ts) { #td and ts are from the f_comb function
  print(paste0("hr of day: ", hod))
  htemps <- (td) * sin((pi / 12) * (hod - 3)) + ts
}

dat.files <- list.files(path, pattern = model, recursive=TRUE, full.names = TRUE)
dat.files <- grep('xml', dat.files, invert = TRUE, value = TRUE)
tmin.files <- grep("tasmin", dat.files, value=TRUE)
tmax.files <- grep("tasmax", dat.files, value=TRUE)

# yrDayCt - vector of days in each year to deal with leap years
yrDayCt <- numeric()
j = 1
for (i in as.numeric(year_range)) {
  year <- year_range[i]
  dates <- seq(as.Date(paste0(i, "-01-01")), as.Date(paste0(i, "-12-31")),by="1 day")
  yrDayCt[j] <- length(dates)
  j = j+1
}

# use the yrDayCt vector to subset tmin and tmax for each year in year_range
# JDay.SH <- (92-2):(306)
# JDay.NH <- JDay.SH + 180
subsetter <- function(t, start, end, e) {
  r <- rast(t, lyrs = start:end, win = e)
}

j <- i <- start <- 1
for (i in 1:length(yrDayCt)) {
  year <- year_range[[i]]
  start_j <- start + get(paste0("JDay.", hem))[1] -1
  end <- start + yrDayCt[i] - 1
  end_j <- start_j + tail(get(paste0("JDay.", hem)), 1)
  
  tmin.in <- subsetter(tmin.files[1], start_j, end_j, e)
  tmax.in <- subsetter(tmax.files[1], start_j, end_j, e)
  # to get to relevant days, need to adjust start and end. Add jday start to start and add jday end and use this for end
  print(tmin.in)
  j = j + 1
  start <- start + yrDayCt[i]
  print(paste0("start: ", start))
  temp <- create_lats(hem, jDayHem, year)
  dl <- temp[[2]]
  sunrise <- temp[[3]]
  sunset <- temp[[4]]
  
  system.time(tcombs <- f_comb(tmax.in, tmin.in))
  temp <- NULL
  gc()
  hodRange <- 4:18
  system.time(t <- lapply(hodRange, FUN = f_hourly_tempWdeltas, tcombs[[1]], tcombs[[2]]))
  t_r <- rast(t)
  
  test <- clamp(t_r, lower = -10, upper = 2, values = FALSE)
}

# test
hourlyFromDailyTemp(tmin, tmax, doy, latitude)
system.time(t <- lapply(hodRange, FUN = hourlyFromDailyTemp, tmin.in, tmax.in, jDayHem, lat))

# notes:
# chill portions are related to the number of hours below a certain temp value in a window of j days
# t is a 36x layer raster with the temperature for one of the hours in hodRange for each rownum (days)
# in t_r so at a particular cellnum you have a column of 

