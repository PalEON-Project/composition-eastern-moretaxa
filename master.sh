#!/bin/bash
# controller script for fitting the PalEON settlement era composition model
# for eastern township data with finer-grained taxonomic resolution
# note that these steps are intended for use on UNIX-like machines and will need to be modified for Windows (and possibly for Mac OS X)

# this is not intended to be run as a full script as later components depend on earlier ones having finished (and the fitting steps take many days); rather it is intended to allow one to run all of the steps of the model fitting/analysis

# this code is being run under R 4.2.0 and with package versioning controlled by packrat
# restore any packages that are not installed on the system
Rscript -e "require(packrat); packrat::restore()"

# modify the contents of the config file to reflect the data versions to be used, relevant directories, and parameters of the MCMC
# in general, it's good to create a version of config, say config_0.1-0, specific to each run and then copy that file to 'config'

# for example:
# \cp config_0.1-0 config

source config

export OMP_NUM_THREADS=1

########################################################################
# create directories    ------------------------------------------------
########################################################################

if [ ! -e $plotDir ]; then
    mkdir $plotDir
fi
if [ ! -e $outputDir ]; then
    mkdir $outputDir
fi
## if [ ! -e $dataDir ]; then
##    mkdir $dataDir
## fi
if [ ! -e $tmpDir ]; then
    mkdir $tmpDir
fi

########################################################################
# obtain eastern township data from Charlie
########################################################################

cd $projectDir

########################################################################
# preprocess eastern township data -------------------------------------
########################################################################

if [ ! -e $dataDir/data_eastern_${productVersion}.Rda ]; then
    ./intersect_towns_cells.R  >& log.intersect_towns_cells 
# creates intersection_${runID}.Rda

    ./build_eastern.R  >& log.build_eastern_${productVersion} &
# this reads intersection_eastern_${runID}.Rda and creates 'data_eastern_${runID}.Rda
fi

########################################################################
# fit Bayesian composition model to eastern data -----------------------
########################################################################

./fit_eastern.R >& log.fit_eastern_${runID} &
# this creates 'PLScomposition_eastern_${runID}_full.nc'
# note that this netCDF has y-values from N to S (contradicting the 
# dim info for the y dim)

########################################################################
# combine runs and subset final output to burned-in samples
########################################################################

# this also flips the N->S orientation of netCDF files to their correct
# S->N to match the y dimension as given in the netCDF

# eastern
numRuns=10
burnin=50000
domain=eastern
./combine_runs.R $burnin $domain $numRuns
# this creates 'PLScomposition_eastern_${runID}.nc'

cp $outputDir/PLScomposition_eastern_${runID}.nc /server/web/share/paciorek/paleon/composition_east_v${productVersion}.nc

########################################################################
# make netCDFs with summary stats (posterior mean, sd, etc.) -----------
########################################################################

## not used for this analysis
## ./summarize_posterior.R 


########################################################################
# make plots -----------------------
########################################################################

## Not used for this analysis (yet)

# download us_alb.zip, water tiff, domain tiff from Wiki
# http://144.92.235.115/dokuwiki/doku.php/public_data%3Brasters
wget $cookieArgs "https://paleon.geography.wisc.edu/lib/exe/fetch.php/data_and_products%3Bus_alb.zip" -O $dataDir/tmp_alb.zip
wget $cookieArgs "https://paleon.geography.wisc.edu/lib/exe/fetch.php/data_and_products%3Bwater_pct_alb_v0.1.tif" -O $dataDir/water.tif 
wget $cookieArgs "https://paleon.geography.wisc.edu/lib/exe/fetch.php/data_and_products%3Bdomain%3Bpaleon_full_alb_v0.1.tif" -O $dataDir/paleonDomain.tif

cd $dataDir
unzip tmp_alb.zip

# this is deprecated as we concentrate on plots of full domain below
# cd $projectDir
# ./plot_eastern.R
# ./plot_western.R

# this creates the mask (paleonMask.nc) for screening out water and non-paleon state cells
./create_mask.R >& log.create_mask &

./plot_full.R
