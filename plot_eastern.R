#!/usr/bin/Rscript

source("config")

runID <- paste0('eastern_', runID)

source(file.path(codeDir, "plot.R"))
source(file.path(codeDir, "set_domain.R"))

require(ncdf4)
require(ggplot2)
require(maptools)
require(rgdal)
require(raster)

usShp <- readShapeLines(file.path(dataDir, 'us_alb.shp'), proj4string=CRS('+init=epsg:3175'))
usShp@data$id <- rownames(usShp@data)
usFortified <- fortify(usShp, region='id')




########################################################################################
## get raw data in form to use
########################################################################################

# rev() flips tif N->S to match netCDF S->N
region = t(as.matrix(raster(file.path(dataDir, 'paleonDomain.tif'))))[easternDomainX, rev(easternDomainY)]
water = t(as.matrix(raster(file.path(dataDir, 'water.tif'))))[easternDomainX, rev(easternDomainY)]
# t()  manipulates matrix so plots correctly W-E and N-S in R

region[region %in% c(2,3,5,6,11,12)] <- NA
water[water == 100] <- NA
mask = is.na(region)
maskWater = is.na(water)

#waterE = paleon.water[easternDomainX, sort(180-easternDomainY)]
#regE = paleon.reg[easternDomainX, sort(180- easternDomainY)]

load(file.path(dataDir, paste0('data_', runID, '.Rda')))

load(file.path(dataDir, paste0('intersection_eastern_', productVersion, '.Rda')))

coord <- expand.grid(X = xGrid[easternDomainX], Y = yGrid[yRes+1-rev(easternDomainY)])
  

finalNcdfName <- paste0('PLScomposition_', runID, '.nc')

ncdfPtr <- nc_open(file.path(outputDir, finalNcdfName))
test <- ncvar_get(ncdfPtr, "Oak", c(1, 1, 1), c(-1, -1, -1))

nSamples <- dim(test)[3]

nCells <- m1*m2

if(nCells != prod(dim(test)[1:2]))
  stop("nCells does not match first dimension of netCDF file.")


preds <- array(0, c(nCells, nTaxa, nSamples))
for(p in 1:nTaxa) {
  preds[ , p, ] <- ncvar_get(ncdfPtr, taxa$taxonName[p], c(1, 1, 1), c(-1, -1, -1))
}

attributes(preds)$dimnames[[2]] <- gsub("/", "ZZZ", taxa$taxonName) 

pm <- apply(preds, c(1, 2), 'mean')
psd <- apply(preds, c(1, 2), 'sd')

pm[mask, ] <- NA
psd[mask, ] <- NA
 

#########################################################################################
## plot the data, plots in large grid
#########################################################################################

propBreaks = c(0, 0.01, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1)

figHgt = 20
figWth = 22

pdf(file.path(outputDir, paste0('PLScomposition_', runID, '_fits.pdf')), height = figHgt, width = figWth)
make_veg_map(data = pm, breaks = propBreaks, coords = coord, legendName = 'fitted proportions', map_data = usFortified, facet = TRUE)
dev.off()

psd[psd > .25] = 0.25
psdBreaks = c(0, 0.01, 0.03, 0.05, 0.075, 0.10, 0.15, 0.2, 0.25)
#psd[psd > .3] = 0.3
#psdBreaks = c(0, 0.01, 0.03, 0.05, 0.10, 0.15, 0.2, 0.25, 0.3)
  
pdf(file.path(outputDir, paste0('PLScomposition_', runID, '_uncertainty.pdf')), height = figHgt, width = figWth)
make_veg_map(data = psd, breaks = psdBreaks, coords = coord, legendName = 'std. error', map_data = usFortified, facet = TRUE)
dev.off()


# plot raw data as colored polygons
# it would be nice to have this as a function, but there is a lot of pre-processing...

nTowns <- dim(townCellOverlap)[1]
easternDataDir <- "eastern"
ohioDataDir <- "ohio"
eastern_townships <- readOGR(file.path(dataDir, easternDataDir), paste0(easternVersionID, 'polygonsver', easternVersion))
ohio_townships <- readOGR(file.path(dataDir, ohioDataDir), paste0('OH', ohioVersionID, 'polygonsver', ohioVersion))
ohio_townships <- spTransform(ohio_townships, CRSobj=CRS('+init=epsg:3175'))  # transform to Albers

raw <- matrix(0, nrow = nTowns, ncol = nTaxa)
for(p in 1:nTaxa) {
  tbl <- table(data$town[data$taxon == p])
  raw[as.numeric(names(tbl)) , p] <- tbl
}
raw <- raw / rowSums(raw)
attributes(raw)$dimnames[[2]] <- gsub("/", "ZZZ", taxa$taxonName) 

for(p in seq_len(nTaxa))
  raw[ , p] <- cut(raw[ , p], propBreaks, include.lowest = TRUE, labels = FALSE)

east_fort <- fortify(eastern_townships)
ohio_fort <- fortify(ohio_townships)

east_idMap <- data.frame(id = sort(eastern_townships@data$ID), town = seq_along(eastern_townships))
ohio_idMap <- data.frame(id = sort(ohio_townships@data$ID), town = length(eastern_townships) + seq_along(ohio_townships))

east_fort <- merge(east_fort, east_idMap, by.x = 'id', by.y = 'id', all.x = TRUE, all.y = FALSE)
east_fort <- east_fort[order(as.numeric(east_fort$id), east_fort$order), ]

ohio_fort <- merge(ohio_fort, ohio_idMap, by.x = 'id', by.y = 'id', all.x = TRUE, all.y = FALSE)
ohio_fort <- ohio_fort[order(as.numeric(ohio_fort$id), ohio_fort$order), ]

fort <- rbind(east_fort, ohio_fort)
fort <- cbind(fort, raw[fort$town,])
names(fort)[1:2] <- c('X', 'Y')

#ggplot(fort, aes(long, lat, group = id)) + geom_polygon(aes(fill = Oak))

names(fort) <- gsub("ZZZ", "/", names(fort))
taxon_dat_long <- melt(fort, c('X', 'Y', 'order', 'hole', 'piece', 'group', 'id' ))

breaks = propBreaks
breaklabels <- apply(cbind(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)]), 1, 
                       function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })

legendName = "raw proportions"

d <- ggplot(taxon_dat_long, aes(X, Y, group = town)) +
  geom_polygon(aes(fill = as.factor(value))) + 
  scale_fill_manual(values = rev(terrain.colors(length(breaks))), labels = breaklabels, name = legendName) +
  theme(strip.text.x = element_text(size = 16), legend.key.size = unit(1.5, "cm"), legend.text = element_text(size = 16), legend.title = element_text(size = 16)) +
  coord_fixed() +
  facet_wrap( ~ variable, ncol=5)
d <- add_map_albers(plot_obj = d, map_data = usFortified, dat = taxon_dat_long)
d <- theme_clean(d)

pdf(file.path(outputDir, paste0('PLScomposition_', runID, '_rawData.pdf')), height = figHgt, width = figWth)
print(d)
dev.off()


