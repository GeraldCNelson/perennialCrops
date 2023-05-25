# accumulate world chill data!
source("../../R/MCA_utilities.R")

require(terra)
require(rasterVis)
require(raster)

models <- tolower(c("GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL"))
scenarios <- c("historical", "ssp126", "ssp585")

datdir <- '../../chillPortions/worldchill/chill_portions/'

# load chill portions data 
for(i in seq_along(scenarios)) {
  for(j in seq_along(models)) {
    files <- list.files(paste0(datdir, scenarios[i], '/', models[j]), pattern='world', include.dirs = FALSE, full.names = TRUE)
    for(k in seq_along(files)) {
      year <- str_split(files[k], "_", simplify = TRUE)[,4]
      objname <- paste0(scenarios[i], "_", models[j], "_", year, ".dat") %>% str_replace_all('-', '.')
      assign(objname, rast(files[k]))
    }
  }
}

# aggregate to common 1deg resolution

  historical_gfdl.esm4_2000.dat <- aggregate(historical_gfdl.esm4_2000.dat, fact=2)
  historical_ipsl.cm6a.lr_2000.dat <- aggregate(historical_ipsl.cm6a.lr_2000.dat, fact=2)
  historical_mpi.esm1.2.hr_2000.dat <- aggregate(historical_mpi.esm1.2.hr_2000.dat, fact=2)
  historical_mri.esm2.0_2000.dat <- aggregate(historical_mri.esm2.0_2000.dat, fact=2)
  historical_ukesm1.0.ll_2000.dat <- aggregate(historical_ukesm1.0.ll_2000.dat, fact=2)
  
  ssp126_gfdl.esm4_2050.dat <- aggregate(ssp126_gfdl.esm4_2050.dat, fact=2)
  ssp126_ipsl.cm6a.lr_2050.dat <- aggregate(ssp126_ipsl.cm6a.lr_2050.dat, fact=2)
  ssp126_mpi.esm1.2.hr_2050.dat <- aggregate(ssp126_mpi.esm1.2.hr_2050.dat, fact=2)
  ssp126_mri.esm2.0_2050.dat <- aggregate(ssp126_mri.esm2.0_2050.dat, fact=2)
  ssp126_ukesm1.0.ll_2050.dat <- aggregate(ssp126_ukesm1.0.ll_2050.dat, fact=2)
  
  ssp585_gfdl.esm4_2050.dat <- aggregate(ssp585_gfdl.esm4_2050.dat, fact=2)
  ssp585_ipsl.cm6a.lr_2050.dat <- aggregate(ssp585_ipsl.cm6a.lr_2050.dat, fact=2)
  ssp585_mpi.esm1.2.hr_2050.dat <- aggregate(ssp585_mpi.esm1.2.hr_2050.dat, fact=2)
  ssp585_mri.esm2.0_2050.dat <- aggregate(ssp585_mri.esm2.0_2050.dat, fact=2)
  # ssp585_ukesm1.0.ll_2050.dat <- aggregate(ssp585_ukesm1.0.ll_2050.dat, fact=2)
  
  # ssp126_gfdl.esm4_2090.dat <- aggregate(ssp126_gfdl.esm4_2090.dat, fact=2)
  # ssp126_ipsl.cm6a.lr_2090.dat <- aggregate(ssp126_ipsl.cm6a.lr_2090.dat, fact=2)
  ssp126_mpi.esm1.2.hr_2090.dat <- aggregate(ssp126_mpi.esm1.2.hr_2090.dat, fact=2)
  ssp126_mri.esm2.0_2090.dat <- aggregate(ssp126_mri.esm2.0_2090.dat, fact=2)
  ssp126_ukesm1.0.ll_2090.dat <- aggregate(ssp126_ukesm1.0.ll_2090.dat, fact=2)
  
  # ssp585_gfdl.esm4_2090.dat <- aggregate(ssp585_gfdl.esm4_2090.dat, fact=2)
  # ssp585_ipsl.cm6a.lr_2090.dat <- aggregate(ssp585_ipsl.cm6a.lr_2090.dat, fact=2)
  ssp585_mpi.esm1.2.hr_2090.dat <- aggregate(ssp585_mpi.esm1.2.hr_2090.dat, fact=2)
  ssp585_mri.esm2.0_2090.dat <- aggregate(ssp585_mri.esm2.0_2090.dat, fact=2)


# get ensemble means

historical_models <- ls(pattern="historical")
historical_mean <- lapply(historical_models, get)
historical_mean <- lapply(historical_mean, mean)
historical_mean <- mean(stack(lapply(historical_mean, raster)))

ssp126_2050_models <- ls(pattern="ssp126") %>% grep("2050", ., value = TRUE)
ssp126_2050_mean <- lapply(ssp126_2050_models, get)
ssp126_2050_mean <- lapply(ssp126_2050_mean, mean)
ssp126_2050_mean <- mean(stack(lapply(ssp126_2050_mean, raster)))

ssp585_2050_models <- ls(pattern="ssp585") %>% grep("2050", ., value = TRUE)
ssp585_2050_mean <- lapply(ssp585_2050_models, get)
ssp585_2050_mean <- lapply(ssp585_2050_mean, mean)
ssp585_2050_mean <- mean(stack(lapply(ssp585_2050_mean, raster)))

ssp126_2090_models <- ls(pattern="ssp126") %>% grep("2090", ., value = TRUE)
ssp126_2090_mean <- lapply(ssp126_2090_models, get)
ssp126_2090_mean <- lapply(ssp126_2090_mean, mean)
ssp126_2090_mean <- mean(stack(lapply(ssp126_2090_mean, raster)))

ssp585_2090_models <- ls(pattern="ssp585") %>% grep("2090", ., value = TRUE)
ssp585_2090_mean <- lapply(ssp585_2090_models, get)
ssp585_2090_mean <- lapply(ssp585_2090_mean, mean)
ssp585_2090_mean <- mean(stack(lapply(ssp585_2090_mean, raster)))

levelplot(historical_mean, par.settings=plasmaTheme(), margin=FALSE, main="Historical chill portions")
levelplot(ssp126_2090_mean, par.settings=plasmaTheme(), margin=FALSE, main="ssp126 2090 chill portions")
levelplot(ssp585_2090_mean, par.settings=plasmaTheme(), margin=FALSE, main="ssp585 2090 chill portions")


# get baseline differences

getBaselineDiff <- function(scenario, future_year, models) {
  out.list <- vector('list', length(models))
  names(out.list) <- models
  for(i in seq_along(models)) {
    future <- get(paste0(scenario, "_", models[i], "_", future_year, ".dat"))
    historical <- get(paste0("historical", "_", models[i], "_", 2000, ".dat"))
    out.list[[i]] <- raster(future-historical)
  }
  return(out.list)
}

models <- str_replace_all(models, "-",".")

# get baseline differences and construct objects

Z_ssp126_2050 <- stack(getBaselineDiff(scenarios[2], c(2050,2090)[1], models))
Z_ssp126_2090 <- stack(getBaselineDiff(scenarios[2], c(2050,2090)[2], models))
Z_ssp585_2050 <- stack(getBaselineDiff(scenarios[3], c(2050,2090)[1], models))
Z_ssp585_2090 <- stack(getBaselineDiff(scenarios[3], c(2050,2090)[2], models))

# plot difference in chill portions for each GCM
levelplot(Z_ssp126_2050, par.settings = RdBuTheme(), main="ssp126 diff 2000-2050", at=seq(-100,100, len=101))
levelplot(Z_ssp585_2050, par.settings = RdBuTheme(), main="ssp585 diff 2000-2050", at=seq(-100,100, len=101))
levelplot(Z_ssp126_2090, par.settings = RdBuTheme(), main="ssp126 diff 2000-2090", at=seq(-100,100, len=101))
levelplot(Z_ssp585_2090, par.settings = RdBuTheme(), main="ssp585 diff 2000-2090", at=seq(-100,100, len=101))

# plot ensemble mean difference in chill portions
levelplot(mean(Z_ssp126_2050), par.settings = RdBuTheme(), main="ssp126 change 2000-2050", at=seq(-120,120, len=101), margin=FALSE)
levelplot(calc(Z_ssp126_2050,sd), par.settings = plasmaTheme(), main="SD ssp126 change 2000- 2050", at=seq(0,50, len=101), margin=FALSE)
levelplot(mean(Z_ssp585_2050), par.settings = RdBuTheme(), main="ssp585 change 2000-2050", at=seq(-120,120, len=101), margin=FALSE)
levelplot(calc(Z_ssp585_2050,sd), par.settings = plasmaTheme(), main="SD ssp585 change 2000-2050", at=seq(0,50, len=101), margin=FALSE)
levelplot(mean(Z_ssp126_2090), par.settings = RdBuTheme(), main="ssp126 change 2000-2090", at=seq(-120,120, len=101), margin=FALSE)
levelplot(calc(Z_ssp126_2090,sd), par.settings = plasmaTheme(), main="SD ssp126 change 2000-2090", at=seq(0,50, len=101), margin=FALSE)
levelplot(mean(Z_ssp585_2090), par.settings = RdBuTheme(), main="ssp585 change 2000-2090", at=seq(-120,120, len=101), margin=FALSE)
levelplot(calc(Z_ssp585_2090,sd), par.settings = plasmaTheme(), main="SD ssp585 change 2000-2090", at=seq(0,50, len=101), margin=FALSE)


## extract site level data ##

sites <- read.csv('../../chillPortions/data/Chill_sites.csv')
sites.sp <- vect(sites[,2:3])

historical_gfdl.esm4_2000.sites <- data.frame(t(extract(historical_gfdl.esm4_2000.dat, sites.sp)))
colnames(historical_gfdl.esm4_2000.sites) <- sites[,1]
historical_gfdl.esm4_2000.sites <- historical_gfdl.esm4_2000.sites[-1,]
historical_gfdl.esm4_2000.sites$year <- 1991:2009
historical_gfdl.esm4_2000.sites$ssp <- 'historical'
historical_gfdl.esm4_2000.sites$model <- "gfdl.esm4"

historical_ipsl.cm6a.lr_2000.sites <- data.frame(t(extract(historical_ipsl.cm6a.lr_2000.dat, sites.sp)))
colnames(historical_ipsl.cm6a.lr_2000.sites) <- sites[,1]
historical_ipsl.cm6a.lr_2000.sites <- historical_ipsl.cm6a.lr_2000.sites[-1,]
historical_ipsl.cm6a.lr_2000.sites$year <- 1991:2009
historical_ipsl.cm6a.lr_2000.sites$ssp <- 'historical'
historical_ipsl.cm6a.lr_2000.sites$model <- "ipsl.cm6a.lr"

historical_mpi.esm1.2.hr_2000.sites <- data.frame(t(extract(historical_mpi.esm1.2.hr_2000.dat, sites.sp)))
colnames(historical_mpi.esm1.2.hr_2000.sites) <- sites[,1]
historical_mpi.esm1.2.hr_2000.sites <- historical_mpi.esm1.2.hr_2000.sites[-1,]
historical_mpi.esm1.2.hr_2000.sites$year <- 1991:2009
historical_mpi.esm1.2.hr_2000.sites$ssp <- 'historical'
historical_mpi.esm1.2.hr_2000.sites$model <- "mpi.esm1.2.hr"

historical_mri.esm2.0_2000.sites <- data.frame(t(extract(historical_mri.esm2.0_2000.dat, sites.sp)))
colnames(historical_mri.esm2.0_2000.sites) <- sites[,1]
historical_mri.esm2.0_2000.sites <- historical_mri.esm2.0_2000.sites[-1,]
historical_mri.esm2.0_2000.sites$year <- 1991:2009
historical_mri.esm2.0_2000.sites$ssp <- 'historical'
historical_mri.esm2.0_2000.sites$model <- "mri.esm2.0"

historical_ukesm1.0.ll_2000.sites <- data.frame(t(extract(historical_ukesm1.0.ll_2000.dat, sites.sp)))
colnames(historical_ukesm1.0.ll_2000.sites) <- sites[,1]
historical_ukesm1.0.ll_2000.sites <- historical_ukesm1.0.ll_2000.sites[-1,]
historical_ukesm1.0.ll_2000.sites$year <- 1991:2009
historical_ukesm1.0.ll_2000.sites$ssp <- 'historical'
historical_ukesm1.0.ll_2000.sites$model <- "ukesm1.0.ll"

# ssp126 2050

ssp126_gfdl.esm4_2050.sites <- data.frame(t(extract(ssp126_gfdl.esm4_2050.dat, sites.sp)))
colnames(ssp126_gfdl.esm4_2050.sites) <- sites[,1]
ssp126_gfdl.esm4_2050.sites <- ssp126_gfdl.esm4_2050.sites[-1,]
ssp126_gfdl.esm4_2050.sites$year <- 2041:2059
ssp126_gfdl.esm4_2050.sites$ssp <- 'ssp126'
ssp126_gfdl.esm4_2050.sites$model <- "gfdl.esm4"

ssp126_ipsl.cm6a.lr_2050.sites <- data.frame(t(extract(ssp126_ipsl.cm6a.lr_2050.dat, sites.sp)))
colnames(ssp126_ipsl.cm6a.lr_2050.sites) <- sites[,1]
ssp126_ipsl.cm6a.lr_2050.sites <- ssp126_ipsl.cm6a.lr_2050.sites[-1,]
ssp126_ipsl.cm6a.lr_2050.sites$year <- 2041:2059
ssp126_ipsl.cm6a.lr_2050.sites$ssp <- 'ssp126'
ssp126_ipsl.cm6a.lr_2050.sites$model <- "ipsl.cm6a.lr"

ssp126_mpi.esm1.2.hr_2050.sites <- data.frame(t(extract(ssp126_mpi.esm1.2.hr_2050.dat, sites.sp)))
colnames(ssp126_mpi.esm1.2.hr_2050.sites) <- sites[,1]
ssp126_mpi.esm1.2.hr_2050.sites <- ssp126_mpi.esm1.2.hr_2050.sites[-1,]
ssp126_mpi.esm1.2.hr_2050.sites$year <- 2041:2059
ssp126_mpi.esm1.2.hr_2050.sites$ssp <- 'ssp126'
ssp126_mpi.esm1.2.hr_2050.sites$model <- "mpi.esm1.2.hr"

ssp126_mri.esm2.0_2050.sites <- data.frame(t(extract(ssp126_mri.esm2.0_2050.dat, sites.sp)))
colnames(ssp126_mri.esm2.0_2050.sites) <- sites[,1]
ssp126_mri.esm2.0_2050.sites <- ssp126_mri.esm2.0_2050.sites[-1,]
ssp126_mri.esm2.0_2050.sites$year <- 2041:2059
ssp126_mri.esm2.0_2050.sites$ssp <- 'ssp126'
ssp126_mri.esm2.0_2050.sites$model <- "mri.esm2.0"

ssp126_ukesm1.0.ll_2050.sites <- data.frame(t(extract(ssp126_ukesm1.0.ll_2050.dat, sites.sp)))
colnames(ssp126_ukesm1.0.ll_2050.sites) <- sites[,1]
ssp126_ukesm1.0.ll_2050.sites <- ssp126_ukesm1.0.ll_2050.sites[-1,]
ssp126_ukesm1.0.ll_2050.sites$year <- 2041:2059
ssp126_ukesm1.0.ll_2050.sites$ssp <- 'ssp126'
ssp126_ukesm1.0.ll_2050.sites$model <- "ukesm1.0.ll"

# ssp585 2050

# ssp585_gfdl.esm4_2050.sites <- data.frame(t(extract(ssp585_gfdl.esm4_2050.dat, sites.sp)))
# colnames(ssp585_gfdl.esm4_2050.sites) <- sites[,1]
# ssp585_gfdl.esm4_2050.sites <- ssp585_gfdl.esm4_2050.sites[-1,]
# ssp585_gfdl.esm4_2050.sites$year <- 2041:2059
# ssp585_gfdl.esm4_2050.sites$ssp <- 'ssp585'
# ssp585_gfdl.esm4_2050.sites$model <- "gfdl.esm4"

ssp585_ipsl.cm6a.lr_2050.sites <- data.frame(t(extract(ssp585_ipsl.cm6a.lr_2050.dat, sites.sp)))
colnames(ssp585_ipsl.cm6a.lr_2050.sites) <- sites[,1]
ssp585_ipsl.cm6a.lr_2050.sites <- ssp585_ipsl.cm6a.lr_2050.sites[-1,]
ssp585_ipsl.cm6a.lr_2050.sites$year <- 2041:2059
ssp585_ipsl.cm6a.lr_2050.sites$ssp <- 'ssp585'
ssp585_ipsl.cm6a.lr_2050.sites$model <- "ipsl.cm6a.lr"

ssp585_mpi.esm1.2.hr_2050.sites <- data.frame(t(extract(ssp585_mpi.esm1.2.hr_2050.dat, sites.sp)))
colnames(ssp585_mpi.esm1.2.hr_2050.sites) <- sites[,1]
ssp585_mpi.esm1.2.hr_2050.sites <- ssp585_mpi.esm1.2.hr_2050.sites[-1,]
ssp585_mpi.esm1.2.hr_2050.sites$year <- 2041:2059
ssp585_mpi.esm1.2.hr_2050.sites$ssp <- 'ssp585'
ssp585_mpi.esm1.2.hr_2050.sites$model <- "mpi.esm1.2.hr"

ssp585_mri.esm2.0_2050.sites <- data.frame(t(extract(ssp585_mri.esm2.0_2050.dat, sites.sp)))
colnames(ssp585_mri.esm2.0_2050.sites) <- sites[,1]
ssp585_mri.esm2.0_2050.sites <- ssp585_mri.esm2.0_2050.sites[-1,]
ssp585_mri.esm2.0_2050.sites$year <- 2041:2059
ssp585_mri.esm2.0_2050.sites$ssp <- 'ssp585'
ssp585_mri.esm2.0_2050.sites$model <- "mri.esm2.0"

ssp585_ukesm1.0.ll_2050.sites <- data.frame(t(extract(ssp585_ukesm1.0.ll_2050.dat, sites.sp)))
colnames(ssp585_ukesm1.0.ll_2050.sites) <- sites[,1]
ssp585_ukesm1.0.ll_2050.sites <- ssp585_ukesm1.0.ll_2050.sites[-1,]
ssp585_ukesm1.0.ll_2050.sites$year <- 2041:2059
ssp585_ukesm1.0.ll_2050.sites$ssp <- 'ssp585'
ssp585_ukesm1.0.ll_2050.sites$model <- "ukesm1.0.ll"

# ssp126 2090

ssp126_gfdl.esm4_2090.sites <- data.frame(t(extract(ssp126_gfdl.esm4_2090.dat, sites.sp)))
colnames(ssp126_gfdl.esm4_2090.sites) <- sites[,1]
ssp126_gfdl.esm4_2090.sites <- ssp126_gfdl.esm4_2090.sites[-1,]
ssp126_gfdl.esm4_2090.sites$year <- 2081:2099
ssp126_gfdl.esm4_2090.sites$ssp <- 'ssp126'
ssp126_gfdl.esm4_2090.sites$model <- "gfdl.esm4"

ssp126_ipsl.cm6a.lr_2090.sites <- data.frame(t(extract(ssp126_ipsl.cm6a.lr_2090.dat, sites.sp)))
colnames(ssp126_ipsl.cm6a.lr_2090.sites) <- sites[,1]
ssp126_ipsl.cm6a.lr_2090.sites <- ssp126_ipsl.cm6a.lr_2090.sites[-1,]
ssp126_ipsl.cm6a.lr_2090.sites$year <- 2081:2099
ssp126_ipsl.cm6a.lr_2090.sites$ssp <- 'ssp126'
ssp126_ipsl.cm6a.lr_2090.sites$model <- "ipsl.cm6a.lr"

ssp126_mpi.esm1.2.hr_2090.sites <- data.frame(t(extract(ssp126_mpi.esm1.2.hr_2090.dat, sites.sp)))
colnames(ssp126_mpi.esm1.2.hr_2090.sites) <- sites[,1]
ssp126_mpi.esm1.2.hr_2090.sites <- ssp126_mpi.esm1.2.hr_2090.sites[-1,]
ssp126_mpi.esm1.2.hr_2090.sites$year <- 2081:2099
ssp126_mpi.esm1.2.hr_2090.sites$ssp <- 'ssp126'
ssp126_mpi.esm1.2.hr_2090.sites$model <- "mpi.esm1.2.hr"

ssp126_mri.esm2.0_2090.sites <- data.frame(t(extract(ssp126_mri.esm2.0_2090.dat, sites.sp)))
colnames(ssp126_mri.esm2.0_2090.sites) <- sites[,1]
ssp126_mri.esm2.0_2090.sites <- ssp126_mri.esm2.0_2090.sites[-1,]
ssp126_mri.esm2.0_2090.sites$year <- 2081:2099
ssp126_mri.esm2.0_2090.sites$ssp <- 'ssp126'
ssp126_mri.esm2.0_2090.sites$model <- "mri.esm2.0"

ssp126_ukesm1.0.ll_2090.sites <- data.frame(t(extract(ssp126_ukesm1.0.ll_2090.dat, sites.sp)))
colnames(ssp126_ukesm1.0.ll_2090.sites) <- sites[,1]
ssp126_ukesm1.0.ll_2090.sites <- ssp126_ukesm1.0.ll_2090.sites[-1,]
ssp126_ukesm1.0.ll_2090.sites$year <- 2081:2099
ssp126_ukesm1.0.ll_2090.sites$ssp <- 'ssp126'
ssp126_ukesm1.0.ll_2090.sites$model <- "ukesm1.0.ll"

# ssp585 2090 

# ssp585_gfdl.esm4_2090.sites <- data.frame(t(extract(ssp585_gfdl.esm4_2090.dat, sites.sp)))
# colnames(ssp585_gfdl.esm4_2090.sites) <- sites[,1]
# ssp585_gfdl.esm4_2090.sites <- ssp585_gfdl.esm4_2090.sites[-1,]
# ssp585_gfdl.esm4_2090.sites$year <- 2081:2099
# ssp585_gfdl.esm4_2090.sites$ssp <- 'ssp585'
# ssp585_gfdl.esm4_2090.sites$model <- "gfdl.esm4"

ssp585_ipsl.cm6a.lr_2090.sites <- data.frame(t(extract(ssp585_ipsl.cm6a.lr_2090.dat, sites.sp)))
colnames(ssp585_ipsl.cm6a.lr_2090.sites) <- sites[,1]
ssp585_ipsl.cm6a.lr_2090.sites <- ssp585_ipsl.cm6a.lr_2090.sites[-1,]
ssp585_ipsl.cm6a.lr_2090.sites$year <- 2081:2099
ssp585_ipsl.cm6a.lr_2090.sites$ssp <- 'ssp585'
ssp585_ipsl.cm6a.lr_2090.sites$model <- "ipsl.cm6a.lr"

ssp585_mpi.esm1.2.hr_2090.sites <- data.frame(t(extract(ssp585_mpi.esm1.2.hr_2090.dat, sites.sp)))
colnames(ssp585_mpi.esm1.2.hr_2090.sites) <- sites[,1]
ssp585_mpi.esm1.2.hr_2090.sites <- ssp585_mpi.esm1.2.hr_2090.sites[-1,]
ssp585_mpi.esm1.2.hr_2090.sites$year <- 2081:2099
ssp585_mpi.esm1.2.hr_2090.sites$ssp <- 'ssp585'
ssp585_mpi.esm1.2.hr_2090.sites$model <- "mpi.esm1.2.hr"

ssp585_mri.esm2.0_2090.sites <- data.frame(t(extract(ssp585_mri.esm2.0_2090.dat, sites.sp)))
colnames(ssp585_mri.esm2.0_2090.sites) <- sites[,1]
ssp585_mri.esm2.0_2090.sites <- ssp585_mri.esm2.0_2090.sites[-1,]
ssp585_mri.esm2.0_2090.sites$year <- 2081:2099
ssp585_mri.esm2.0_2090.sites$ssp <- 'ssp585'
ssp585_mri.esm2.0_2090.sites$model <- "mri.esm2.0"

ssp585_ukesm1.0.ll_2090.sites <- data.frame(t(extract(ssp585_ukesm1.0.ll_2090.dat, sites.sp)))
colnames(ssp585_ukesm1.0.ll_2090.sites) <- sites[,1]
ssp585_ukesm1.0.ll_2090.sites <- ssp585_ukesm1.0.ll_2090.sites[-1,]
ssp585_ukesm1.0.ll_2090.sites$year <- 2081:2099
ssp585_ukesm1.0.ll_2090.sites$ssp <- 'ssp585'
ssp585_ukesm1.0.ll_2090.sites$model <- "ukesm1.0.ll"

sitedata.historical <- rbind(historical_gfdl.esm4_2000.sites, historical_ipsl.cm6a.lr_2000.sites, historical_mpi.esm1.2.hr_2000.sites, historical_mri.esm2.0_2000.sites, historical_ukesm1.0.ll_2000.sites)
sitedata.126.2050 <- rbind(ssp126_gfdl.esm4_2050.sites, ssp126_ipsl.cm6a.lr_2050.sites, ssp126_mpi.esm1.2.hr_2050.sites, ssp126_mri.esm2.0_2050.sites, ssp126_ukesm1.0.ll_2050.sites)
sitedata.585.2050 <- rbind(ssp585_ipsl.cm6a.lr_2050.sites, ssp585_mpi.esm1.2.hr_2050.sites, ssp585_mri.esm2.0_2050.sites, ssp585_ukesm1.0.ll_2050.sites)
sitedata.126.2090 <- rbind(ssp126_gfdl.esm4_2090.sites, ssp126_ipsl.cm6a.lr_2090.sites, ssp126_mpi.esm1.2.hr_2090.sites, ssp126_mri.esm2.0_2090.sites, ssp126_ukesm1.0.ll_2090.sites)
sitedata.585.2090 <- rbind(ssp585_ipsl.cm6a.lr_2090.sites, ssp585_mpi.esm1.2.hr_2090.sites, ssp585_mri.esm2.0_2090.sites, ssp585_ukesm1.0.ll_2090.sites)

sitedata.all <- rbind(sitedata.historical, sitedata.126.2050, sitedata.585.2050, sitedata.126.2090, sitedata.585.2090)
sitedata.all <- sitedata.all[,c('ssp','model','year','Manjimup','Ithaca')]

write_csv(sitedata.all, '../../chillPortions/review/worldchill_sitedata.csv')
