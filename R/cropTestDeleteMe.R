dat.dir <- paste0('climdata/')
model <- tolower("GFDL-ESM4")
scenario = "historical"
#  dat.files <- list.files(dat.dir, pattern=model, recursive=TRUE, full.names = TRUE)
dat.files <- list.files(dat.dir, pattern=model, recursive=TRUE, full.names = TRUE)
dat.files <- grep('xml', dat.files, invert = TRUE, value = TRUE)
# dat.files <- grep(year_range, dat.files, invert = TRUE, value = TRUE)
tmin.files <- grep("tasmin", dat.files, value=TRUE)
tmax.files <- grep("tasmax", dat.files, value=TRUE)
ext_NH <- ext(c(-180, 180, 0, 90))
ext_SH <- ext(c(-180, 180, -60, 0))
system.time(tmin.in <- rast(tmin.files[1]))


period <- 92:306
system.time(tmin.north.win <- rast(tmin.files[1], win = ext_NH))
system.time(tmin_north_period <- tmin.north.win[[period]])
tmin_north_period
m <- as.matrix(tmin_north_period)

test_nh <- brick(nrows=180, ncols=720, xmn=-180, xmx=180, ymn=0, ymx=90, nl=length(period), crs = "+proj=longlat +datum=WGS84 +no_defs" )
plot(tmin_north_period, 1)
system.time(test_nh[] <- m)
plot(tmin_north_period, 1)
tmin.north <- crop(tmin.in, ext_NH)
system.time(tmin.north.win <- tmin.north.win * 1)

system.time(tmin.north <- rast(tmin.files[1], win = ext_NH, lyrs = JDay.north))
system.time(m <- as.matrix(tmin.north))
ext_hem <- ext(c(-180, 180, 0, 90))
period <- JDay.north <- 92:306 + 180

brick_empty <- brick(nrows=180, ncols=720, xmn=-180, xmx=180, ymn=ext_hem[3], ymx=ext_hem[4], nl=length(period), crs = "+proj=longlat +datum=WGS84 +no_defs")
system.time(brick_empty[] <- m)
