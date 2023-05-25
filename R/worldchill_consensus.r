# accumulate world chill data!
# consensus maps
source("../../R/MCA_utilities.R")

require(terra)
require(rasterVis)
require(raster)

models <- tolower(c("GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL"))
scenarios <- c("historical", "ssp126", "ssp585")

# set this to where the data is
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

# collect models by scenario / year range

historical <- lapply(ls(pattern='historical'), get)
ssp126_2050 <- lapply(ls(pattern='ssp126') %>% grep('2050', ., value=TRUE), get)
ssp126_2090 <- lapply(ls(pattern='ssp126') %>% grep('2090', ., value=TRUE), get)
ssp585_2050 <- lapply(ls(pattern='ssp585') %>% grep('2050', ., value=TRUE), get)
ssp585_2090 <- lapply(ls(pattern='ssp585') %>% grep('2090', ., value=TRUE), get)

#### binarise chill portions stacks by a CP threshold of choice ####
outdir <- '../../chillPortions/outputs/binarised_chill/almond'
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  # helper function to write binarised rasters to file. modify as necessary
  write_models <- function(x, list, models, ccp_time) {
      writeRaster(list[[x]], paste0('../../chillPortions/outputs/binarised_chill/almond/', ccp_time, '_', models[[x]], '.tif'), overwrite=TRUE)
  }

  # binarise by CP > 9 (for almonds and write raster stacks to disk
  historical_threshold <- lapply(seq_along(models), function(x) historical[[x]]>9)
  names(historical_threshold) <- models
  lapply(seq_along(models), write_models, historical_threshold, models, 'historical')
  
  ssp126_2050_threshold <- lapply(seq_along(models), function(x) ssp126_2050[[x]]>9)
  names(ssp126_2050_threshold) <- models
  lapply(seq_along(models), write_models, ssp126_2050_threshold, models, 'ssp126_2050')
  
  ssp126_2090_threshold <- lapply(seq_along(models), function(x) ssp126_2090[[x]]>9)
  names(ssp126_2090_threshold) <- models
  lapply(seq_along(models), write_models, ssp126_2090_threshold, models, 'ssp126_2090')
  
  ssp585_2050_threshold <- lapply(seq_along(models), function(x) ssp585_2050[[x]]>9)
  names(ssp585_2050_threshold) <- models
  lapply(seq_along(models), write_models, ssp585_2050_threshold, models, 'ssp585_2050')
  
  ssp585_2090_threshold <- lapply(seq_along(models), function(x) ssp585_2090[[x]]>9)
  names(ssp585_2090_threshold) <- models
  lapply(seq_along(models), write_models, ssp585_2090_threshold, models, 'ssp585_2090')

# calculate safe chill (10th percentile thresholding at >9CP, 1 layer per timeseries)

  # calculate 10th percentile across stack
  binarise_safechill <- function(x, thresh,...) {
    safe <- quantile(x, 0.1,na.rm=TRUE)
    safe[safe < thresh] <- 0
    safe[safe >= thresh] <- 1
    return(safe)
  }
  
historical_safechill <- stack(lapply(seq_along(models), function(x) raster(app(historical[[x]], binarise_safechill, thresh=9))))
ssp126_2050_safechill <- stack(lapply(seq_along(models), function(x) raster(app(ssp126_2050[[x]], binarise_safechill, thresh=9))))
ssp126_2090_safechill <- stack(lapply(seq_along(models), function(x) raster(app(ssp126_2090[[x]], binarise_safechill, thresh=9))))
ssp585_2050_safechill <- stack(lapply(seq_along(models), function(x) raster(app(ssp585_2050[[x]], binarise_safechill, thresh=9))))
ssp585_2090_safechill <- stack(lapply(seq_along(models), function(x) raster(app(ssp585_2090[[x]], binarise_safechill, thresh=9))))
