# auxiliary code for creating a netCDF file and returning a pointer to it
require(ncdf4)

makeAlbersNetCDF <- function(name = NULL, units = '', longname = '', fn = NULL, dir = '', x, y, taxa, numSamples) {
  if(is.null(fn)) stop("makeAlbersNetCDF: requires filename as argument.")

  x_dim <-  ncdim_def("x", "meters_east", x, longname = 'x coordinate of grid cell centroid in Albers projection (Great Lakes St Lawrence Albers [Proj4 +init=epsg:3175])')
  y_dim <-  ncdim_def("y", "meters_north", y, longname = 'y coordinate of grid cell centroid in Albers projection (Great Lakes St Lawrence Albers [Proj4 +init=epsg:3175])')
  draw_dim <- ncdim_def("sample", "number", 1:numSamples, longname = "MCMC sample")
  vars <- list()
  length(vars) == length(taxa)
  for(k in seq_along(taxa)) {
    vars[[k]] <- ncvar_def(name = taxa[k], dim=list(x_dim, y_dim, draw_dim), units = units, longname = paste0(longname, " for taxon ", taxa[k]), prec="double")
  }
  ncdfPtr <- nc_create(file.path(dir, fn), vars)
  nc_close(ncdfPtr)
  invisible(NULL)
}

makeAlbersNetCDFsummary <- function(name = NULL, units = '', longname = '', fn = NULL, dir = '', dims, taxa) {
  if(is.null(fn)) stop("makeAlbersNetCDFsummary: requires filename as argument.")
  vars <- list()
  length(vars) == length(taxa)
  for(k in seq_along(taxa)) {
    vars[[k]] <- ncvar_def(name = taxa[k], dim = dims, units = units, longname = paste0(longname, " for taxon ", taxa[k]), prec="double")
  }
  ncdfPtr <- nc_create(file.path(dir, fn), vars)
  nc_close(ncdfPtr)
  invisible(NULL)
}
