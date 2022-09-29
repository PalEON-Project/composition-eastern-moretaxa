source("config")

source(file.path(codeDir, "set_domain.R"))

require(ncdf4)

require(maptools)
require(rgdal)
require(raster)
require(ggplot2)

paleonStates <- c("Minnesota","Wisconsin","Michigan", "Illinois", "Indiana", "Ohio", "Pennsylvania", "New Jersey", "New York", "Massachusetts", "Rhode Island", "Connecticut", "New Hampshire", "Vermont", "Maine")

usShp <- readShapeLines(file.path(dataDir, 'us_alb.shp'), proj4string=CRS('+init=epsg:3175'))
usShp@data$id <- rownames(usShp@data)
usFortified <- fortify(usShp, region='id')

paleonShp <- usShp[usShp$STATE_NAME %in% paleonStates, ]
paleonShp@data$id <- rownames(paleonShp@data)
paleonFortified <- fortify(paleonShp, region='id')

########################################################################################
## get raw data in form to use
########################################################################################

# rev() flips from tif N->S orientation to netCDF S->N
region = t(as.matrix(raster(file.path(dataDir, 'paleonDomain.tif'))))[ , yRes:1]
water = t(as.matrix(raster(file.path(dataDir, 'water.tif'))))[ , yRes:1]
# t()  manipulates matrix so plots correctly W-E and N-S in R

# region[region %in% c(2,3,5,6,11,12)] <- NA
water[water == 100] <- NA
mask = is.na(region)
maskWater = is.na(water)

# western data

load(file.path(dataDir, paste0('data_western_', runID, '.Rda')))
taxa_west <- taxa

raw_west <- matrix(0, nrow = nCells, ncol = nTaxa)
for(p in 1:nTaxa) {
  tbl <- table(data$cell[data$taxon == p])
  raw_west[as.numeric(names(tbl)) , p] <- tbl
}
total <- rowSums(raw_west)
raw_west[total == 0] = NA
total[total == 0] <- 1
raw_west <- raw_west / total
dimnames(raw_west)[[2]] <- taxa$taxonName  



# eastern data
load(file.path(dataDir, paste0('data_eastern_', runID, '.Rda')))
taxa_east <- taxa
load(file.path(dataDir, paste0('intersection_eastern_', productVersion, '.Rda')))

nTowns <- dim(townCellOverlap)[1]
easternDataDir <- 'eastern'
ohioDataDir <- 'ohio'
eastern_townships <- readOGR(file.path(dataDir, easternDataDir), paste0(easternVersionID, 'polygons_v', easternVersion))
ohio_townships <- readOGR(file.path(dataDir, ohioDataDir), paste0('OH', ohioVersionID, 'polygons_v', ohioVersion))
ohio_townships <- spTransform(ohio_townships, CRSobj=CRS('+init=epsg:3175'))  # transform to Albers

raw_east <- matrix(0, nrow = nTowns, ncol = nTaxa)
for(p in 1:nTaxa) {
  tbl <- table(data$town[data$taxon == p])
  raw_east[as.numeric(names(tbl)) , p] <- tbl
}
raw_east <- raw_east / rowSums(raw_east)
attributes(raw_east)$dimnames[[2]] <- taxa$taxonName


taxaNames <- unique(taxa_east$taxonName, taxa_west$taxonName)


# put NA for masked areas for taxa not modeled in west
taxaMissingInWest <- !(taxaNames %in% dimnames(raw_west)[[2]])
old_raw_west <- raw_west
raw_west <- matrix(as.numeric(NA), nrow(old_raw_west), ncol(raw_east))
dimnames(raw_west)[[2]] <- taxaNames
raw_west[ , dimnames(old_raw_west)[[2]]] <- old_raw_west


# to fill in non-NA for cells lacking any trees but not water or out of PalEON states
west_presence <- read.csv(file.path(dataDir, 'western.csv'))
west_presence <- west_presence[west_presence$x <= easternLimitOfWesternDomain, ]
dataPresent <- as.numeric(rowSums(west_presence[ , c(3:ncol(west_presence))]) > 0)
tmp <- matrix(dataPresent, length(westernDomainY), length(westernDomainX))
dataPresent <- c(t(tmp[rev(westernDomainY),]))
rm(west_presence)
for(p in 1:ncol(raw_west)) {
    raw_west[is.na(raw_west[,p]) & dataPresent == 1, p] <- 0
}



raw_west_poly <- rasterToPolygons(raster(ncols = length(westernDomainX), nrows = length(westernDomainY),
                                       xmn = xRange[1], xmx = easternLimitOfWesternDomain + gridSpacing/2,
                         ymn = yRange[1], ymx = yRange[2]),
                      fun=NULL, n=4, na.rm=TRUE, digits=12, dissolve=FALSE)
west_fort <- fortify(raw_west_poly)
names(west_fort)[1:2] <- c('X', 'Y')


west_fort$id <- as.numeric(west_fort$id)
west_fort <- cbind(west_fort, raw_west[west_fort$id, ])

northeast_fort <- fortify(eastern_townships)
northeast_fort$id <- as.numeric(northeast_fort$id)
ohio_fort <- fortify(ohio_townships)
ohio_fort$id <- as.numeric(ohio_fort$id)

northeast_idMap <- data.frame(id = sort(eastern_townships@data$ID), town = seq_along(eastern_townships))
ohio_idMap <- data.frame(id = sort(ohio_townships@data$ID), town = length(eastern_townships) + seq_along(ohio_townships))

northeast_fort <- merge(northeast_fort, northeast_idMap, by.x = 'id', by.y = 'id', all.x = TRUE, all.y = FALSE)
northeast_fort <- northeast_fort[order(as.numeric(northeast_fort$id), northeast_fort$order), ]

ohio_fort <- merge(ohio_fort, ohio_idMap, by.x = 'id', by.y = 'id', all.x = TRUE, all.y = FALSE)
ohio_fort <- ohio_fort[order(as.numeric(ohio_fort$id), ohio_fort$order), ]

east_fort <- rbind(northeast_fort, ohio_fort)
east_fort <- cbind(east_fort, raw_east[east_fort$town,])
east_fort$id <- east_fort$town + max(west_fort$id)
names(east_fort)[2:3] <- c('X', 'Y')
east_fort$town <- NULL
east_fort <- east_fort[ , names(west_fort)]

fort <- rbind(west_fort, east_fort)


#names(fort) <- gsub("ZZZ", "/", names(fort))
taxon_dat_long <- melt(fort, c('X', 'Y', 'order', 'hole', 'piece', 'group', 'id' ))

# strip out chestnut/atl wh cedar as not modeled in west
# max_west_id <- max(west_fort$id)
# taxon_dat_long <- taxon_dat_long[!(taxon_dat_long$id <= max_west_id & taxon_dat_long$variable %in% taxaNames[taxaMissingInWest]), ]

finalNcdfName <- paste0('composition_v', productVersion, '.nc')

ncdfPtr <- nc_open(file.path(outputDir, finalNcdfName))
test <- ncvar_get(ncdfPtr, "Oak", c(1, 1, 1), c(-1, -1, -1))

nCells <- xRes*yRes
nTaxa <- length(taxaNames)
nSamples <- dim(test)[3]


preds <- array(0, c(nCells, nTaxa, nSamples))
#dimnames(preds)[[2]] <- taxa$taxonName
for(p in 1:nTaxa) {
  preds[ , p, ] <- ncvar_get(ncdfPtr, taxa$taxonName[p], c(1, 1, 1), c(-1, -1, -1))
}

attributes(preds)$dimnames[[2]] <- gsub("/", "ZZZ", taxaNames)

pm <- apply(preds, c(1, 2), 'mean')
psd <- apply(preds, c(1, 2), 'sd')


pm[mask, ] <- NA
psd[mask, ] <- NA

coord <- expand.grid(X = xGrid, Y = yGrid)

