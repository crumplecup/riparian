

#' BISECT EXTENT
#'
#' Given an extent `ext`, returns the extents of the bisections in a list.
#'
#' @param ext is an object of class `extent`
#' @param logical indictor of vertical or horizontal bisection, default TRUE
#' @return a list of the extents of the bisections
#' @export

bisecter <- function(ext, vertical=TRUE)	{
  if (class(ext) != 'Extent')	{
    ext <- extent(ext)
  }
  x_dif <- ext[2] - ext[1]
  y_dif <- ext[4] - ext[3]
  if (vertical)	{
    aframe <- c(ext[1], ext[1] + x_dif * .5, ext[3:4])
    bframe <- c(ext[1] + x_dif * .5, ext[2:4])
  }
  if (!vertical)	{
    aframe <- c(ext[1:2], ext[3], ext[3] + y_dif * .5)
    bframe <- c(ext[1:2], ext[3] + y_dif * .5, ext[4])
  }
  return(list(aframe,bframe))
}



#' Build Cover Change Table
#'
#'
#'
#' @param csv is a csv file with cols(id, year, type, obs(1:50))
#' @param year1 is a number matching a value in csv$year
#' @param year2 is a number matching a value in csv$year
#' @param polys is a spatial polygons object (sampling boxes)
#' @param lots is a spatial polygons object (riparian taxlots)
#' @param permits is a csv file with cols(maptaxlot, record_no)
#' @return prints a csv file with cols(id, change, maptaxlot, permits) to the working directory
#' @export


build_change_table <- function(csv,
                               year1,
                               year2,
                               polys,
                               lots,
                               permits)  {
  dt <- data.table::setDT(csv)
  dt$mn <- apply(dt[,4:53], 1, function(x) sum(x)/100)
  yr1 <- dt[year == year1]
  yr2 <- dt[year == year2]
  chng <- rep(NA,nrow(yr2))
  for (i in 1:nrow(yr2))  {
    try(chng[i] <- yr2$mn[i] - yr1$mn[yr1$id == yr2$id[i]], silent=T)
  }
  ids <- 1:length(polys)
  change <- array(0,length(polys))
  for (i in seq_along(change))  {
    change[i] <- chng[yr2$id %in% as.character(ids[i])]
  }
  mtls <- lookup_lots(polys, lots)
  perms <- lookup_permits(polys, lots, permits)
  dt <- data.table::setDT(data.frame(id = ids, change = change, maptaxlot = mtls[ids], permits = perms[ids]))
  write.csv(dt[order(change)], file = paste0('change_table_', year1, '-', year2,'.csv'))
}



#' COORDINATE LIST FROM SPATIALPOLYGONS
#'
#'
#'
#' @param poly is a `SpatialPolygons` object
#' @return a list coordinates for vertices of polygons in `poly`
#' @export

coord_lis <- function(poly)  {
  coords <- lapply(slot(poly,'polygons'),
                   function(x) lapply(slot(x,'Polygons'),
                                      function(y) slot(y, 'coords')))
  coords[[1]]
}



#' MERGE RASTERS TO FILL EXTENT
#'
#' Given a `SpatialPolygons` object and a directory filepath character string,
#' returns merged rasters from `dir` over extent of `poly`.
#'
#' @param poly is a `SpatialPolygons` object
#' @param dir is a directory filepath character string
#' @param band is the band to extract from multiband imagery
#' @return merged rasters from `dir` over extent of `poly`
#' @export



fill_extent <- function(poly, dir, band)  {
  og_dir <- getwd()
  setwd(dir)
  files <- get_rasters(dir)
  poly <- match_crs(poly, raster::raster(files[1]))
  ext <- raster::extent(poly)
  hit <- 0
  
  for(i in seq_along(files))  {
    ras <- raster::raster(files[i])
    if (!is.null(unlist(ras[poly,])))  {
      hit <- c(hit,i)
    }
  }
  hit <- hit[-1]
  ras <- raster::raster(files[hit[1]], band=band)
  if (length(hit)>1)  {
    for (i in 2:length(hit)) {
      ras <- raster::mosaic(ras,raster::raster(files[hit[i]], band=band),fun=max, tolerance = 1)
    }
  }
  setwd(og_dir)
  raster::crop(ras, ext)
}


#' box from extent
#' 
#' produce polygon representation of extent
#' 
#' @param ras is a spatial object (raster)
#' @return a Spatial Polygon drawn around the extent of `ras`
ext_box <- function(ras, crs_ref)  {
  ext <- raster::extent(ras)
  crs_ref <- raster::crs(ras)
  box <- matrix(ncol = 2, nrow = 5)
  box[1,] <- ext[c(1,3)]
  box[2,] <- ext[c(1,4)]
  box[3,] <- ext[c(2,4)]
  box[4,] <- ext[c(2,3)]
  box[5,] <- ext[c(1,3)]
  spatialize(box, crs_ref)
}
  
fill_extent1 <- function(poly, dir, band)  {
  og_dir <- getwd()
  setwd(dir)
  files <- get_rasters(dir)
  poly <- match_crs(poly, raster::raster(files[1]))
  ext <- raster::extent(poly)
  hit <- 0
  
  for(i in seq_along(files))  {
    if (in_extent(ext,raster::raster(files[i])))  {
      hit <- c(hit,i)
    }
  }
  hit <- hit[-1]
  ras <- raster::raster(files[hit[1]], band=band)
  if (length(hit)>1)  {
    for (i in 2:length(hit)) {
      ras <- raster::mosaic(ras,raster::raster(files[hit[i]], band=band),fun=max, tolerance = .1)
    }
  }
  setwd(og_dir)
  raster::crop(ras, ext)
}





#' Get colors from 4 band raster
#'
#' Returns data.table of mean ndvi and rgb values for a sampling box
#'
#' @param poly is the sampling box
#' @param dir is the path to the rasters
#' @return a dt with cols (ndvi, red, grn, blu)
#' @export


get_cols_4band <- function(poly, dir)  {
  red <- read_sample(poly, dir, 1)
  grn <- read_sample(poly, dir, 2)
  blu <- read_sample(poly, dir, 3)
  nir <- read_sample(poly, dir, 4)
  ndvi <- (nir - red) / (nir + red)
  setDT(data.frame(ndvi = ndvi, red = red, grn = grn, blu = blu))
}



#' Get CRS
#'
#' Get a reference CRS from rasters in a directory
#'
#' @param dir is a path to a directory of rasters
#' @return the CRS character string
#' @export


get_crs <- function(dir)  {
  og_dir <- getwd()
  setwd(dir)
  files <- get_rasters(dir)
  crs <- raster::crs(raster::raster(files[1]))
  setwd(og_dir)
  crs
}



#' get raster names
#'
#' scan a path directory for raster file names, includes .tif, .hdr
#'
#' @param path is a character string specifying a file path to a directory of rasters
#' @return a character vector of raster file names
#' @export



get_rasters <- function(path = getwd())  {
  og_dir <- getwd()
  setwd(path)
  files <- dir()
  files <- files[!grepl('.xml', files)]
  files <- files[!grepl('.ovr', files)]
  files <- files[!grepl('.ini', files)]
  tifs <- files[grepl('.tif', files)]
  hdrs <- files[grepl('.hdr', files)]
  setwd(og_dir)
  c(tifs,hdrs)
}





#' IN EXTENT LOGICAL TEST
#'
#' Tests whether the point (or extent) \code{pt} falls entirely within the extent of a raster \code{ras},
#' returning TRUE or FALSE.
#'
#' @param pt is a coordinate vector or extent class object
#' @param so is a spatial object
#' @return TRUE if pt is within extent of ras, else FALSE
#' @export

in_extent <- function(pt, so)	{
  is_in=FALSE
  ext <- raster::extent(so)
  
  if (class(pt) == 'numeric')	{
    if (pt[1] >= ext[1] &
        pt[1] <= ext[2] &
        pt[2] >= ext[3] &
        pt[2] <= ext[4])	{ is_in <- TRUE }
  }
  
  if (class(pt) == 'Extent')	{
    if (pt[1] >= ext[1] &
        pt[2] <= ext[2] &
        pt[3] >= ext[3] &
        pt[4] <= ext[4])	{ is_in <- TRUE }
  }
  
  
  return(is_in)
}



#' Look up MapTaxlot`
#'
#' Given a list of polygons, returns character vector of all maptaxlot numbers
#' within each polygon.
#'
#' @param so is a spatial object (sampling boxes)
#' @param lots is an spdf of polygons (bc taxlots)
#' @return the MapTaxlot numbers of each input polygon in a character vector
#' @export


lookup_lots <- function(so, lots)  {
  mtls <- vector(length = length(so), mode = 'character')
  polys <- sp::spTransform(lots, raster::crs(so))
  for (i in seq_along(so)) {
    lot <- raster::crop(polys, raster::extent(so[i,]))
    mtls[i] <- squash(levels(factor(lot$MapTaxlot)))
  }
  mtls
}



#' Look up permits
#'
#'
#'
#' @param so is a spatial polygons object (sampling boxes)
#' @param lots is a spatial polygons object (taxlots)
#' @param permits is a csv file with cols(maptaxlot, record_no)
#' @return a character vector listing permit record numbers for taxlots in polys


lookup_permits <- function(so, lots, permits){
  names(permits[,1]) <- c('record_no')
  names(permits[,5]) <- c('maptaxlot')
  perms <- vector(length = length(so), mode = 'character')
  polys <- sp::spTransform(lots, raster::crs(so))
  for (i in seq_along(so)) {
    if (class(so) %in% c('SpatialPolygons', 'SpatialPolygonsDataFrame')) {
      lot <- raster::crop(polys, raster::extent(so[i,]))
    }
    if (class(so) %in% c('RasterLayer', 'RasterStack', 'RasterBrick'))  {
      lot <- raster::crop(polys, so[i,], mask = T)
    }
    map_nos <- levels(factor(lot$MapTaxlot))
    perms[i] <- squash(permits$record_no[permits$maptaxlot %in% map_nos])
  }
  perms
}




#' Convert coordinates to SpatialLines
#'
#' Wrapper for sp workflow of nesting spatial objects in lists like russian dolls.
#'
#' @param coords is a matrix with cols(x,y)
#' @param crs_ref is a crs character string
#' @return a SpatialLines object
#' @export

make_line <- function(coords,crs_ref)  {
  sp::SpatialLines(
    list(
      sp::Lines(
        list(
          sp::Line(coords)
        ), ID = 'first'
      )
    ), crs_ref
  )
  # coords %>% (sp::Line) %>% list %>% (sp::Lines(ID='first')) %>% list %>% (sp::SpatialLines(crs_ref))
}




#' MATCH CRS
#'
#' Make sure your extents overlap by converting polygons
#' into a common reference crs.
#'
#' @param polys is a spatial polygons objects
#' @param crs_ref is a character string representing the crs
#' @return a list of spatial polygons projected in the common crs
#' @export



match_crs <- function(polys, crs_ref)  {
  if (class(polys) %in% c('RasterLayer', 'RasterStack'))  {
    polys <- raster::projectRaster(polys, crs_ref)
  }
  if (class(polys) %in% c('SpatialPolygons', 'SpatialPolygonsDataFrame'))  {
    polys <- sp::spTransform(polys, raster::crs(crs_ref))
  }
  polys
}





#' QUAD CUT RASTER
#'
#' Bisects raster into four equal parts, returns the parts in a list.
#'
#' @param ras is a raster
#' @return four equal parts of \code{ras} in a list
#' @export

quadcut <- function(ras)	{
  start <- bisecter(ras)
  left <- bisecter(start[[1]], vertical=F)
  right <- bisecter(start[[2]], vertical=F)
  return(list(left[[1]], left[[2]], right[[1]], right[[2]]))
}



#' Plot Cover Change
#'
#' Prints a graph of cover change.
#'
#' @param csv is a csv with cols (id, year, type, results(1:50))
#' @param year1 is an integer matching a value in csv$year
#' @param year2 is an integer matching a value in csv$year
#' @param title is the character vector to name the output file (png)
#' @param heading is the character vector used as the main title of the plot
#' @param leg_pos is the legend position argument (character)
#' @param y_lab is the label for the y axis (character)
#' @return prints a graph of mean cover change into the working directory
#' @import data.table
#' @export

plot_change <- function(csv,
                        year1,
                        year2,
                        title='cover_change.png',
                        heading=' ',
                        leg_pos='topleft',
                        y_lab = 'No. of Sample Sites with Cover Change of X')  {
  dt <- as.data.table(csv)
  dt$mn <- apply(dt[,4:53],1,function(x) sum(x)/100)
  yr1 <- dt[year == year1]
  yr2 <- dt[year == year2]
  chng <- rep(NA,nrow(yr2))
  for (i in 1:nrow(yr2))  {
    try(chng[i] <- yr2$mn[i] - yr1$mn[yr1$id == yr2$id[i]], silent=T)
  }
  vars <- matrix(0, ncol=5)
  colnames(vars) <- c('mean', 'sd', 'n', 'upr', 'lwr')
  vars[1,1] <- mean(chng[!is.na(chng)])
  vars[1,2] <- sd(chng[!is.na(chng)])
  vars[1,3] <- length(chng[!is.na(chng)])
  vars[1,4] <- vars[1,1] + 1.96 * (vars[1,2] / sqrt(vars[1,3]))
  vars[1,5] <- vars[1,1] - 1.96 * (vars[1,2] / sqrt(vars[1,3]))
  pmf_chng <- density(chng[!is.na(chng)])
  
  png(title, width=9, height=6, units='in', res=300)
  plot(pmf_chng, lwd=2, col='forestgreen',
       main=heading, ylab=y_lab, axes=F)
  axis(1, at=seq(-1,1,.1),
       labels=c('-100%','-90%','-80%', '-70%', '-60%', '-50%', '-40%', '-30%', '-20%', '-10%', '0%',
                '10%', '20%', '30%', '40%','50%','60%','70%','80%','90%','100%'))
  axis(2)
  abline(v=vars[1,4],col='steelblue',lty=2)
  abline(v=vars[1,5],col='steelblue',lty=2)
  abline(v=vars[1,1],col='slategrey')
  legend(leg_pos,legend=c('No. of Samples at X% Change','Mean Cover Change',
                          'Upper & Lower 95% CI on Mean Change'),
         fill = c('forestgreen','slategrey','steelblue'))
  dev.off()
  print(vars)
}


#' Plot Cover Extent
#'
#' Prints a graph of cover extent.
#'
#' @param csv is a csv with cols (id, year, type, results(1:50))
#' @param title is the character vector to name the output file (png)
#' @param heading is the character vector used as the main title of the plot
#' @param leg_pos is the legend position argument (character)
#' @param y_lab is the label for the y axis (character)
#' @return prints a graph of mean cover extent into the working directory
#' @import data.table
#' @export
plot_cover <- function(csv,
                       title='cover_change.png',
                       heading=' ',
                       leg_pos='topleft',
                       y_lab = 'Proportion of Samples at or below X% Cover')  {
  dt <- as.data.table(csv)
  dt$mn <- apply(dt[,4:53],1,function(x) sum(x)/100)
  vars <- matrix(0, ncol=5)
  colnames(vars) <- c('mean', 'sd', 'n', 'upr', 'lwr')
  vars[1,1] <- mean(dt$mn[!is.na(dt$mn)])
  vars[1,2] <- sd(dt$mn[!is.na(dt$mn)])
  vars[1,3] <- length(dt$mn[!is.na(dt$mn)])
  vars[1,4] <- vars[1,1] + 1.96 * (vars[1,2] / sqrt(vars[1,3]))
  vars[1,5] <- vars[1,1] - 1.96 * (vars[1,2] / sqrt(vars[1,3]))
  
  png(title, width=9, height=6, units='in', res=300)
  plot(staTools::cdf(dt$mn), lwd=2, col='forestgreen', type = 'l',
       main=heading, xlab='Cover Extent', ylab=y_lab, axes=F)
  axis(1, at=seq(0,1,.1),
       labels=c('0%', '10%', '20%', '30%', '40%','50%','60%','70%','80%','90%','100%'))
  axis(2)
  abline(v=vars[1,4],col='steelblue',lty=2)
  abline(v=vars[1,5],col='steelblue',lty=2)
  abline(v=vars[1,1],col='slategrey')
  legend(leg_pos,legend=c('No. of Samples at X% Change','Mean Cover Extent',
                          'Upper & Lower 95% CI on Mean Extent'),
         fill = c('forestgreen','slategrey','steelblue'))
  dev.off()
  print(vars)
}



#' predict change
#' 
#' Given a directories of matching rasters from different years,
#' calculate change by subtracting year 1 from year 2.
#' 
#' @param year1_path is a character string indicating path to year 1 rasters
#' @param year2_path is a character string indicating path to year 2 rasters
#' @param out_path is a character string indicating path to output change rasters
#' @return prints .tif files to `out_path` of year 2 minus year 1
#' @export
pred_change <- function(year1_path, year2_path, out_path)  {
  og_dir <- getwd()
  year1_files <- get_rasters(year1_path)
  year2_files <- get_rasters(year2_path)
  for (i in seq_along(year1_files))  {
    setwd(year1_path)
    year1 <- raster::raster(year1_files[i])
    setwd(year2_path)
    year2 <- raster::raster(year2_files[i])
    year1 <- raster::projectRaster(year1, year2)
    chng <- year2 - year1
    setwd(out_path)
    raster::writeRaster(chng, year1_files[i], 'GTiff')
  }
  setwd(og_dir)
}




#' predicted change report
#' 
#' Given a path to change rasters produced by `pred_change()`,
#' outputs a csv to the working directory of change per taxlot,
#' including permits for each taxlot.
#' 
#' @param chng_path is a path for directory of change rasters output by `pred_change()`
#' @param year1_path is a character string path for directory of year 1 rasters
#' @param year2_path is a character string path for directory of year 2 rasters
#' @param permits is a list of permits issued by the county
#' @param lots is a spatial polygons object (maptaxlots)
#' @param title is the name given to the output csv file
#' @return a csv to the working directory showing change per taxlot, and listing permits for each taxlot
#' @import data.table
#' @export
pred_change_report <- function(chng_path, 
                               year1_path,
                               year2_path,
                               permits = permits_13to18,
                               lots = prc_lots, 
                               title = 'pred_chng_table.csv')  {
  og_dir <- getwd()
  files <- get_rasters(year2_path)
  fil_len <- length(files)
  ids <- 1:fil_len
  
  cov_yr1 <- array(0, fil_len)
  area_yr1 <- array(0, fil_len)
  cov_yr2 <- array(0, fil_len)
  area_yr2 <- array(0, fil_len)
  
  chng_adj <- array(0, fil_len)
  
  
  chng <- array(0, fil_len)
  area <- vector(length = fil_len, mode = 'numeric')
  rat <- vector(length = fil_len, mode = 'numeric')
  mtls <- vector(length = fil_len, mode = 'character')
  perms <- vector(length = fil_len, mode = 'character')
  cover <- vector(length = fil_len, mode = 'numeric')
  for (i in seq_along(files))  {
    print(paste0('Reading ', i, ' of ', fil_len, ' files...'))
    ras <- raster::raster(paste0(chng_path, files[i]))
    vals <- raster::values(ras)
    chng[i] <- sum(vals[!is.na(vals)]) / 2
    area[i] <- length(vals[!is.na(vals)])
    cover[i] <- (sum(vals[!is.na(vals)])/2) / area[i]
    lot <- sp::spTransform(lots, raster::crs(ras))
    lot <- raster::crop(lot, ras, mask = T)
    map_nos <- levels(factor(lot$MapTaxlot))
    mtls[i] <- squash(map_nos)
    perms[i] <- squash(permits$record_no[permits$maptaxlot %in% map_nos])

    ras <- raster::raster(paste0(year1_path, files[i]))
    vals <- raster::values(ras)
    vals[vals<0] <- 0
    vals[vals>2] <- 2
    area_yr1[i] <- length(vals[!is.na(vals)])
    cov_yr1[i] <- sum(vals[!is.na(vals)]) / 2
    
    ras <- raster::raster(paste0(year2_path, files[i]))
    vals <- raster::values(ras)
    vals[vals<0] <- 0
    vals[vals>2] <- 2
    area_yr2[i] <- length(vals[!is.na(vals)])
    cov_yr2[i] <- sum(vals[!is.na(vals)]) / 2
    
  }
  rat[area > 0] <- chng[area > 0] / area[area > 0]
  cov_yr1[area_yr1 > 0] <- cov_yr1[area_yr1 > 0] / area_yr1[area_yr1 > 0]
  cov_yr2[area_yr2 > 0] <- cov_yr2[area_yr2 > 0] / area_yr2[area_yr2 > 0]
  chng_adj <- cov_yr2 - cov_yr1
  
  tot_yr1 <- sum(cov_yr1 * area_yr1) / sum(area_yr1)
  tot_yr2 <- sum(cov_yr2 * area_yr2) / sum(area_yr2)
  tot_chng <- sum(chng_adj * area_yr1) / sum(area_yr1)
  dt <- setDT(data.frame(ids = ids, 
                         year1 = cov_yr1, 
                         year2 = cov_yr2,
                         change = chng_adj,
                         acres = area / 43560,
                         ratio = rat,
                         lots = mtls, 
                         permits = perms))
  setkey(dt, ratio)
  setwd(og_dir)
  write.csv(dt, title)
  print(paste0('Predicted cover extent year 1:  ', round(tot_yr1, 2)))
  print(paste0('Predicted cover extent year 2:  ', round(tot_yr2, 2)))
  print(paste0('Predicted cover change:  ', round(tot_chng, 2)))
  
  
}




#' before and after plots
#' 
#' for a set of rasters, print images for year 1 and year 2
#' 
#' @param chng_path is a path for directory of change rasters output by `pred_change()`
#' @param year1_path is a character string path for directory of year 1 rasters
#' @param year2_path is a character string path for directory of year 2 rasters
#' @param img_path is a character string path of directory to output results
#' @param buff is a spatial polygons object (stream buffer layer)
#' @param lots is a spatial polygons object (taxlots)
#' @return plots to the working directory showing before and after
#' @export

before_and_after <- function(chng_path, 
                             year1_path,
                             year2_path,
                             img_path,
                             buff = prc_buff,
                             lots = prc_lots)  {
  og_dir <- getwd()
  setwd(chng_path)
  files <- get_rasters()
  for (i in seq_along(files))  {
    setwd(chng_path)
    ras <- raster::raster(files[i])
    crs1 <- get_crs(year1_path)
    ras <- raster::projectRaster(ras, crs = crs1)
    box <- ext_box(ras)
    r <- fill_extent(box, year1_path, 1)
    g <- fill_extent(box, year1_path, 2)
    b <- fill_extent(box, year1_path, 3)
    year1 <- raster::stack(r, g, b)
    buff1 <- sp::spTransform(buff, crs1)
    buff1 <- raster::crop(buff1, box)
    lot1 <- sp::spTransform(lots, crs1)
    lot1 <- raster::crop(lot1, box)
    crs2 <- get_crs(year2_path)
    ras <- raster::projectRaster(ras, crs = crs2)
    box <- ext_box(ras)
    r <- fill_extent(box, year2_path, 1)
    g <- fill_extent(box, year2_path, 2)
    b <- fill_extent(box, year2_path, 3)
    year2 <- raster::stack(r, g, b)
    buff2 <- sp::spTransform(buff, crs2)
    buff2 <- raster::crop(buff2, box)
    lot2 <- sp::spTransform(lots, crs2)
    lot2 <- raster::crop(lot2, box)
    setwd(img_path)
    png(paste0('site_', i, '_before', '.png'))
    raster::plotRGB(year1)
    lines(lot1, col = 'coral')
    lines(buff1, lwd = 2, col = 'slateblue')
    dev.off()
    png(paste0('site_', i, '_after', '.png'))
    raster::plotRGB(year2)
    lines(lot2, col = 'coral')
    lines(buff2, lwd = 2, col = 'slateblue')
    dev.off()
  }
  setwd(og_dir)
}





#' plot predicted cover prime
#'
#' Predicts cover extent within a polygon and returns the results as a raster plot
#'
#' @param poly is a spatial polygons object (taxlot or sampling box)
#' @param in_path is a character string indicating the path to a directory to rasters
#' @param out_path is a character string indicating the path to the output directory
#' @param rgb_path is the path to 3-band rgb data
#' @param cir_path is the path to 3-band cir data
#' @param buff is the 50-ft riparian buffer polygon for `poly`
#' @return a raster plot of predicted cover masked by `poly`
#' @export

pred_cover1 <- function(poly,
                       in_dir = NULL, 
                       out_dir = NULL, 
                       rgb_dir = NULL, 
                       cir_dir = NULL, 
                       buffer = prc_buff)  {
    if(!is.null(rgb_dir)) {
      red <- fill_extent(poly, rgb_dir, 1)
      grn <- fill_extent(poly, rgb_dir, 2)
      blu <- fill_extent(poly, rgb_dir, 3)
    }
    if(!is.null(cir_dir))  {
      nir <- fill_extent(poly, cir_dir, 1)
    } else {
      nir <- fill_extent(poly, in_dir, 4)
      red <- fill_extent(poly, in_dir, 1)
      grn <- fill_extent(poly, in_dir, 2)
      blu <- fill_extent(poly, in_dir, 3)
    }
    
    ndvi <- (nir - red) / (nir + red)
    dt <- setDT(
      data.frame(ndvi = raster::values(ndvi),
                 red = raster::values(red),
                 grn = raster::values(grn),
                 blu = raster::values(blu)))
    pred <- predict(mod18, newdata = dt)
    ras <- ndvi
    raster::values(ras) <- pred
    poly <- match_crs(poly, ras)
    ras <- raster::mask(ras, poly)
    buffer <- match_crs(buffer, ras)
    ras <- raster::mask(ras, buffer)
    ras
}


#' plot predicted cover
#'
#' Predicts cover extent within a polygon and returns the results as a raster plot
#'
#' @param poly is a spatial polygons object (taxlot or sampling box)
#' @param in_path is a character string indicating the path to a directory to rasters
#' @param out_path is a character string indicating the path to the output directory
#' @param rgb_path is the path to 3-band rgb data
#' @param cir_path is the path to 3-band cir data
#' @param buff is the 50-ft riparian buffer polygon for `poly`
#' @return a raster plot of predicted cover masked by `poly`
#' @export

pred_cover <- function(polys,
                       in_path = NULL, 
                       out_path = NULL, 
                       rgb_path = NULL, 
                       cir_path = NULL, 
                       buff = prc_buff)  {
  for (i in seq_along(polys))  {
    ras <- pred_cover1(poly = polys[i,], 
                       in_dir = in_path, 
                       out_dir = out_path,
                       rgb_dir = rgb_path,
                       cir_dir = cir_path,
                       buffer = buff) 
    
    raster::writeRaster(ras, paste0(out_path,'site_',i,'.tif'), 'GTiff', overwrite = T)
  }
}

#' predicted cover report
#' 
#' produces cover extent table from predicted cover rasters
#' 
#' @param in_path is the character string path directory to rasters
#' @param taxlots is a polygon spatial object of maptaxlots
#' @return prints a predicted cover table to the working directory
#' @export

pred_cover_report <- function(in_path, taxlots = prc_lots)  {
  files <- get_rasters(in_path)
  cover <- array(0, length(files))
  area <- array(0, length(files))
  mtls <- vector(length(files), mode = 'character')
  for (i in seq_along(files))  {
    ras <- raster::raster(paste0(in_path, files[i]))
    vals <- raster::values(ras)
    vals[vals<0] <- 0
    vals[vals>2] <- 2
    cover[i] <- sum(vals[!is.na(vals)])
    area[i] <- length(vals[!is.na(vals)])
    ras <- match_crs(ras, taxlots)
    mtls[i] <- squash(lookup_lots(ras, taxlots))
  }
  pct_cov <- mapply(function(a, b) (sum(a)/2) / b, a = cover, b = area)
  tot_cov <- (sum(cover)/2) * (area/(sum(area)))
  dt <- data.table::setDT(data.frame(lot = mtls, cover = pct_cov, area = area))
  write.csv(dt, 'pred_cov.csv')
  print(paste0('Predicted Cover Extent:  ', tot_cov, '%'))
  
}

#' assign raster
#' 
#' wrapper for assigning vals to raster
#' 
#' @param ras is a raster
#' @param vals are vals for raster
#' @return raster `ras` with values `vals`
#' @export
assign_raster <- function(ras, vals)  {
  new <- ras
  values(new) <- vals
  new
}


#' Print Sample Plots
#'
#' Given a shapefile of sampling boxes, and a path to a directory of 4-band ortho-imagery
#' prints rgb and/or ndvi plots of the sampling area with numbered boxes superimposed.
#' While the function prints both rgb and ndvi by default, the user can turn either off
#' by setting either the \code{rgb} or \code{ndvi} arguments to \code{FALSE}.
#'
#' @param in_path is a character vector containing the path to the orthoimagery
#' @param out_path is a character vector specifying the output directory for the plots
#' @param polys is a SpatialPolygonsDataFrame object (the sample boxes shapefile)
#' @param rgb is a logical indicator, prints rgb plots if true
#' @param ndvi is a logical indicator, prints ndvi plots if true
#' @return prints an rgb and/or ndvi plot into the \code{out_path} directory
#' @export


plot_samples <- function(in_path, out_path, polys = samples, 
                         print_rgb = TRUE, print_ndvi = TRUE)  {
  og_wd <- getwd()
  setwd(in_path)
  files <- get_rasters(in_path)
  crs_ref <- raster::crs(raster::raster(files[1]))
  polys <- sp::spTransform(polys, crs_ref)
  
  for (i in seq_along(polys))	{
    setwd(in_path)
    frame <- raster::extent(polys[i,])
    hit <- 0
    
    for (j in seq_along(files))	{
      if (in_extent(frame, raster::raster(files[j])))	{
        hit <- c(hit,j)
      }
    }
    hit <- hit[-1]
    
    if (length(hit) == 0)  {
      print(paste0('Failed to find raster in extent for sample ',i))
    }
    
    if (length(hit) > 0)  {
      if (print_rgb) rgb <- raster::brick(files[hit[1]])
      if (print_ndvi) {
        nir <- raster::raster(files[hit[1]], band = 4)
        red <- raster::raster(files[hit[1]], band = 1)
      }
    }
    
    if (length(hit) > 1)	{
      for (j in 2:length(hit))	{
        if (print_rgb) rgb <- raster::mosaic(rgb, raster::brick(files[hit[j]]), fun = max)
        if (print_ndvi)  {
          nir <- raster::mosaic(nir, raster::raster(files[hit[j]], band = 4), fun = max)
          red <- raster::mosaic(red, raster::raster(files[hit[j]], band = 1), fun = max)
        }
      }
    }
    
    if (print_rgb) rgb <- raster::crop(rgb, frame)
    if (print_ndvi)  {
      nir <- raster::crop(nir, frame)
      red <- raster::crop(red, frame)
      ndvi <- (nir - red) / (nir + red)
    }
    
    setwd(out_path)
    
    if (print_rgb)  {
      png(paste0('rgb_',i,'.png'))
      raster::plotRGB(rgb, main = paste0('Sample Site ',i))
      area <- lapply(methods::slot(polys[i,], 'polygons'),
                     function(x) lapply(methods::slot(x, 'Polygons'),
                                        function(y) methods::slot(y, 'coords')))
      area <- area[[1]]
      for (j in 1:50)  {
        lines(area[[j]])
      }
      lines(area[[1]], col = 'pink', lwd = 3)
      lines(area[[50]], col = 'red', lwd = 3)
      labcords <- matrix(0, nrow = 50, ncol = 2)
      for (j in 1:50)	{
        labcords[j,] <- methods::slot(
          rgeos::gCentroid(
            spatialize(area[[j]], crs_ref)
          ), 'coords')
      }
      text(labcords, as.character(1:50))
      dev.off()
    }
    
    if (print_ndvi)  {
      png(paste0('ndvi_',i,'.png'))
      plot(ndvi, main = paste0('Sample Site ',i))
      area <- lapply(methods::slot(polys[i,], 'polygons'),
                     function(x) lapply(methods::slot(x, 'Polygons'),
                                        function(y) methods::slot(y, 'coords')))
      area <- area[[1]]
      for (j in 1:50)  {
        lines(area[[j]])
      }
      lines(area[[1]], col = 'pink', lwd = 3)
      lines(area[[50]], col = 'red', lwd = 3)
      labcords <- matrix(0, nrow = 50, ncol = 2)
      for (j in 1:50)	{
        labcords[j,] <- methods::slot(
          rgeos::gCentroid(
            spatialize(area[[j]], crs_ref)
          ), 'coords')
      }
      text(labcords, as.character(1:50))
      dev.off()
    }
  }
  setwd(og_wd)
}




#' Convert SpatialPolygon to coordinate matrix
#'
#' Wrapper to convert polygons to coordinate matrices.
#'
#' @param poly is a SpatialPolygons object
#' @return a matrix of coordinates representing the spatial object
#' @export

pull_coords <- function(poly)  {
  coords <- lapply(methods::slot(poly, 'polygons'),
                   function(x) lapply(methods::slot(x, 'Polygons'),
                                      function(y) methods::slot(y, 'coords')))
  coords[[1]][[1]]
}



#' Read band
#'
#' Returns the mean value of a given raster band across a polygon.
#'
#' @param poly is a spatial polygon object
#' @param dir is a directory of raster files
#' @param band is a number specifying the desired raster band
#' @return the mean values of the raster band masked by the polygon
#' @export

read_band <- function(poly, dir, band)  {
  og_dir <- getwd()
  setwd(dir)
  files <- get_rasters(dir)
  poly <- match_crs(poly, raster::raster(files[1]))
  ras <- fill_extent(poly, dir, band)
  setwd(og_dir)
  vals <- raster::values(raster::mask(ras, poly))
  mean(vals[!is.na(vals)])
}



#' Read sample
#'
#' Get mean values of a raster band from each slice of a sampling box
#'
#' @param poly is a spatial polygons object
#' @param dir is a path to a directory of raster files
#' @param band is the desired band of the raster to pull
#' @return the mean values of the raster band across each slice of a sampling box
#' @export

read_sample <- function(poly, dir, band)  {
  crs_ref <- get_crs(dir)
  coords <- coord_lis(poly)
  unlist(
    lapply(coords, function(a, b, c) read_band(spatialize(a, crs_ref), b, c),
           b = dir, c = band))
}



#' Squash Characters Vectors
#'
#' Takes a character vector and returns it as a string.
#'
#' @param vec is a character vector
#' @return vector contents as a string (length 1 char vector)


squash <- function(vec)  {
  char <- vec[1]
  if(length(vec) > 1) {
    for (i in 2:length(vec))  {
      char <- paste(char, vec[i], sep = ' ')
    }
  }
  char
}






#' CONVERT XY COORDS TO SPATIALPOLYGONS
#'
#' Given a matrix of coords with cols (x,y), and a character string representing a coordinate
#' projection, returns a SpatialPolygon defined by coords of \code{mat} in projection \code{crs_ref}.
#'
#' @param mat is a matrix of coords with cols (x,y)
#' @param crs_ref is a character string representing a coordinate projection for proj4string
#' @return a SpatialPolygon defined by coords of mat in projection crs_ref
#' @export


spatialize <- function(mat, crs_ref)	{
  sp::SpatialPolygons(
    list(
      sp::Polygons(
        list(
          sp::Polygon(mat)
        ), ID=1
      )
    ), proj4string = crs_ref
  )
}



#' stack_extent
#'
#' Given a `SpatialPolygons` object and a directory filepath character string,
#' returns merged rasters from `dir` over extent of `poly`.
#' Returns stack of all layers of og rasters
#'
#' @param poly is a `SpatialPolygons` object
#' @param dir is a directory filepath character string
#' @return merged rasters from `dir` over extent of `poly`
#' @export



stack_extent <- function(poly, dir)  {
  files <- get_rasters(dir)
  poly <- match_crs(poly, raster::raster(file.path(dir, files[1])))
  ext <- raster::extent(poly)
  hit <- 0
  
  for(i in seq_along(files))  {
    ras <- raster::raster(file.path(dir, files[i]))
    if (!is.null(unlist(ras[poly,])))  {
      hit <- c(hit,i)
    }
  }
  hit <- hit[-1]
  ras <- raster::stack(file.path(dir, files[hit[1]]))
  if (length(hit)>1)  {
    for (i in 2:length(hit)) {
      ras <- raster::mosaic(ras,raster::stack(file.path(dir, files[hit[i]])), 
                            fun=max, tolerance = 1)
    }
  }
  raster::crop(ras, ext)
}




















