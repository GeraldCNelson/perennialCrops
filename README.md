# Perennial Crops - R and c++ code and data

This dryad site contains all the R and c++ code and data to generate the results in the Global Change Biology, "Assessing temperature-based adaptation limits to climate change of temperate perennial fruit crops". 

# Directory structure
- data-raw - The raw input data are included in this directory. It has the following sub directories
  * climdata - climate data files. Land only files from the ISIMIP project. These data also can be downloaded from 
[https://doi.org/10.48364/ISIMIP.842396.1](https://doi.org/10.48364/ISIMIP.842396.1)
Follow the directions to download for tasmin, tas, and tasmax for each of the five ESMs (GFDL-ESM4, UKESM1-0-LL, MPI-ESM1-2-HR, MRI-ESM2-0, and IPSL-CM6A-LR).
  * crops - information on crop area from [www.earthstat.org](http://www.earthstat.org), chill portions requirements
  * perennials - an excel version of the data table in the supplementary materials of the paper
  * regionInformation - files used to construct the boarders used in the graphics
- data - processed data that are used to generate the final results. Generated by `chillingCalcs.R` and `chillPortions.R`
  * crops - raster mask files for each crop created from the earthstat data.
  * growingDegreeDays - daily growing degree days by ESM and ensemble means across the models, locations where growing degree days are not limiting for one of the crops.
  * perennials - various files related to crop specific suitability
  * runs - files with the results of calculations on the number of days with a value below or above a limit, a run.
- graphics - destination for graphics results

# Order of operations

## Install libraries

Here's some code to install any of the needed packages that are not already installed

`packages <- c("terra", "viridis", "data.table", "flextable", "officer", "crayon", "magrittr", "doParallel", "foreach", "Rcpp", "readxl")`
`installed_packages <- packages %in% rownames(installed.packages())`
`if (any(installed_packages == FALSE)) {  install.packages(packages[!installed_packages]) }`

The two figures for the paper use the pdf graphics driver built in to base R. The code used to produce Figure 1 and Figure 2 also uses the pdfcrop program to trim the white spaces around the edges of the pdf. 

For Mac users, run this code from a terminal to install pdfcrop. This uses the homebrew system (https://brew.sh). A plug - homebrew is a really nice way to manage all sorts of code that the MacOS should have but doesn't. It also works for linux users.
`brew install --cask mactex`

For windows users. I have found these directions but can't verify them

"You need first to install Perl support for Windows (for example https://www.activestate.com/activeperl) and after that use Miktex package manager to install pdfcrop."

More directions here [https://pdfcrop.sourceforge.net]

## Run scripts

- chillingCalcs.R - output the chilling portions file stored in data-raw. This will take a long time and the output files are already available. Run only if new versions are needed.
- perennialCalcs.R  - directions for running are included in the file. Read through the code before sourcing the file

The graphics for the paper are generated by the following files

- GCBPerennials_Figure_1.R
- GCBPerennials_Figure_2.R
- GCBPerennials_Table_4.R
- GCBPerennials_Table_5.R
- GCBPerennials_Table_7.R

Tables 1, 2, 6, and 6 are just word tables that are not generated in R.

