#!/usr/bin/Rscript

# code to intersect Charlie's townships with Simon's 8 km grid
# result is 180x296x(2082) array of proportion intersection
# where 180 is y dimension and 296 is x dimension
# open issue: intersection only uses discrete approximation with 100 points

require(rgdal)
require(raster)

source("config")

####################################################################
# read in shape file info and create raster for PalEON Albers grid
####################################################################

P4S.latlon <- "+proj=longlat +datum=WGS84"
eastern_townships <- readOGR(file.path(dataDir), paste0(easternVersionID, '-', easternVersion), p4s = P4S.latlon)

eastern_townships <- spTransform(eastern_townships, CRSobj=CRS('+init=epsg:3175'))  # transform to Albers

nTowns <- length(eastern_townships) 

source(file.path(codeDir, "set_domain.R"))

rast <- raster(crs = CRS('+init=epsg:3175'),
               xmn = xRange[1], xmx = xRange[2],
               ymn = yRange[1], ymx = yRange[2],
               ncols = xRes, nrows = yRes)

# this raster has rows as y and columns as x and starts with 1,1 in the NW corner, so 180,1 is the SW corner

# NOTE: do not subset as townships@polygons[[index]] as the sorting of this is different than townships[index, ]


####################################################################
# intersect grid with townships ------------------------------------
####################################################################

# intersect with eastern townships

ord <- order(as.numeric(eastern_townships$TARGET_FID))

for(i in seq_len(nrow(eastern_townships))) {
  aa <- rasterize(x = eastern_townships[ord[i], ], y = rast, getCover = TRUE)  
  if(i == 1){
    poly.stack <- stack((aa)) 
  } else {
    poly.stack <- addLayer(poly.stack, aa)
  }
  if(i%%100 == 0) print(i)
}

# result is proportion of each grid cell covered by the polygon
interTmp <- as.array(poly.stack) * 64  # 64 converts to km2 (not really necessary since normalize to proportion below)

# check
if(FALSE) {
  area <- unlist(sapply(eastern_townships@polygons, function(x) x@area/1000000))
  plot(area, apply(interTmp, 3, sum))
}

inter <- array(0, c(xRes, yRes, nTowns))
for(i in 1:nTowns)
  inter[ , , i] <- t(interTmp[ , , i]/sum(interTmp[ , , i]))


# inter goes NW to NE and proceeds in lat bands southward. Last cell is the SE corner
# for plotting, I'll need to reverse the columns

nCells <- xRes*yRes
ids <- 1:nCells
usedIds <- c(matrix(ids, xRes, yRes)[easternDomainX, easternDomainY])

inter <- inter[easternDomainX, easternDomainY, ]

save(inter, usedIds, nTowns, file = file.path(dataDir, paste0('intersection_eastern_', productVersion, '.Rda')))


