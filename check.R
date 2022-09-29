require(rgdal)
require(raster)

eastern_townships <- readOGR('.', '2082-42')

eastern_townships <- spTransform(eastern_townships, CRSobj=CRS('+init=epsg:3175'))  # transform to Albers

counts <- eastern_townships@data