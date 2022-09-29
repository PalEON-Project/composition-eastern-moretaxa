#!/usr/bin/Rscript
# combines runs and subsets to post-burnin samples

# also flips N-S so that values are increasing from S to N
# interim netCDF files go from N to S becaus that is how the tif files
# are set up

source("config")

# grab from input args after sourcing config as we want to overwrite burnin
args <- commandArgs(TRUE)
# now args is a character vector containing the arguments
burnin <- as.numeric(args[1])
domain <- args[2]
numRuns <- as.numeric(args[3])

runIDplus <- paste0(domain, "_", runID)

source(file.path(codeDir, 'netCDF.R'))
source(file.path(codeDir, 'set_domain.R'))
require(ncdf4)

finalThin <- 8  # have to get down to at least 300 samples or get netCDF format violation; this goes down to 250 so dividing 200 per each run evenly

numPerRun <- floor((S-burnin)/(thin*secondThin))/finalThin

load(file.path(dataDir, paste0('data_', runIDplus, '.Rda')))

finalNcdfName <- paste0('PLScomposition_', domain, '_', easternVersionID, '-', easternVersion, '_', runID, '.nc')
makeAlbersNetCDF(name = 'proportion', units = 'unitless (proportion from 0 to 1)', longname = 'relative composition, relative to all tree taxa,', fn = finalNcdfName, dir = outputDir, x = xGrid[easternDomainX], y = yGrid[easternDomainY], taxa = taxa$taxonName, numSamples = numRuns*numPerRun)

finalNcdfPtr <- nc_open(file.path(outputDir, finalNcdfName), write = TRUE)


currentSamples <- floor(S/(thin*secondThin))
finalSamples <- floor((S-burnin)/(thin*secondThin))
wh <- seq(1, finalSamples, by = finalThin)

for(i in seq_len(numRuns)) {
    outputNcdfName <- paste0('PLScomposition_', runIDplus, '-', i, '_full.nc')
    outputNcdfPtr <- nc_open(file.path(outputDir, outputNcdfName))
    
    for(p in seq_len(nTaxa)) {
        tmp <- ncvar_get(outputNcdfPtr, varid = taxa$taxonName[p], start = c(1, 1, (currentSamples - finalSamples + 1)), count = c(-1, -1, finalSamples))
        ## m2:1 flips the y-axis so that y values go from S to N
        ncvar_put(finalNcdfPtr, taxa$taxonName[p], tmp[ , m2:1, wh], start = c(1, 1, 1+(i-1)*numPerRun), count = c(-1, -1, numPerRun))
    }
    nc_close(outputNcdfPtr)
}

nc_close(finalNcdfPtr)

if(FALSE) {
    p <- p+1
    out=ncvar_get(finalNcdfPtr, varid = taxa$taxonName[p], start = c(1, 1, 1), count = c(-1, -1, -1))
    smp1 <- sample(1:180, 20)
    smp2 <- sample(1:140, 20)
    par(mfrow=c(4,5))
    for(i in 1:20)
        ts.plot(out[smp1[i],smp2[i],], ylim=c(0,1))
    
    i=0
    load('~/research/jmac/composition-eastern-moretaxa/data/data_eastern_0.1.Rda')
    pp1=nc_open('PLScomposition_eastern_0.1-1_full.nc')
    pp2=nc_open('PLScomposition_eastern_0.1-2_full.nc')
    pp3=nc_open('PLScomposition_eastern_0.1-3_full.nc')
    pp4=nc_open('PLScomposition_eastern_0.1-4_full.nc')
    smp = sample(1:180, 36)
    smp2 = sample(1:140, 36)
    
    i=i+1
    out1 <- ncvar_get(pp1, varid = taxa[i,2], start = c(1, 1, 1), count = c(-1,-1,-1))
    out2 <- ncvar_get(pp1, varid = taxa[i,2], start = c(1, 1, 1), count = c(-1,-1,-1))
    out3 <- ncvar_get(pp1, varid = taxa[i,2], start = c(1, 1, 1), count = c(-1,-1,-1))
    out4 <- ncvar_get(pp1, varid = taxa[i,2], start = c(1, 1, 1), count = c(-1,-1,-1))
    
    
    par(mfrow = c(6,6))
    print(taxa[i,2])
    for(j in 1:36) {
        vec <- c(out1[smp[j], smp2[j],101:300],
                 out2[smp[j], smp2[j],101:300],
                 out3[smp[j], smp2[j],101:300],
                 out4[smp[j], smp2[j],101:300]
                 )
        print(effectiveSize(vec))
        ts.plot(vec)
    }
}
