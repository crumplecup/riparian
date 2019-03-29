

#' BISECT EXTENT
#'
#' Given an extent \code{ext}, returns the extents of the bisections in a list.
#'
#' @param ext is an object of class \code{extent}
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
#' @param csv is a csv file with cols(sample_id, year, type, obs(1:50))
#' @param year1 is a number matching a value in csv$year
#' @param year2 is a number matching a value in csv$year
#' @param polys is a spatial polygons object (sampling boxes)
#' @param lots is a spatial polygons object (riparian taxlots)
#' @param permits is a csv file with cols(maptaxlot, record_no)
#' @return prints a csv file with cols(id, change, maptaxlot, permits) to the working directory
#' @import data.table
#' @export


build_change_table <- function(csv,
                               year1,
                               year2,
                               polys,
                               lots,
                               permits)  {
  dt <- setDT(csv)
  yr1 <- csv[year == year1]
  yr2 <- csv[year == year2]
  chng <- rep(NA,nrow(yr2))
  for (i in 1:nrow(yr2))  {
    try(chng[i] <- yr2$mn[i] - yr1$mn[yr1$sample_id == yr2$sample_id[i]], silent=T)
  }
  ids <- 1:length(polys)
  change <- array(0,length(polys))
  for (i in seq_along(change))  {
    change[i] <- chng[yr2$sample_id %in% as.character(ids[i])]
  }
  mtls <- lookup_lots(polys, lots)
  perms <- lookup_permits(polys, lots, permits)
  dt <- setDT(data.frame(id = ids, change = change, maptaxlot = mtls[ids], permits = perms[ids]))
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
    if (!is.null(ras[poly,]))  {
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

fill_extent1 <- function(poly, dir, band)  {
  og_dir <- getwd()
  setwd(dir)
  files <- dir()
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
  files <- dir()
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



get_rasters <- function(path)  {
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
#' @param polys is a list of polygons (sampling boxes)
#' @param lots is an spdf of polygons (bc taxlots)
#' @return the MapTaxlot numbers of each input polygon in a character vector
#' @export


lookup_lots <- function(polys, lots)  {
  mtls <- vector(length = length(polys), mode = 'character')
  box <- match_crs(samples, lots)
  for (i in seq_along(samples)) {
    lot <- crop(lots, box[i,])
    mtls[i] <- squash(levels(factor(lot$MapTaxlot)))
  }
  mtls
}



#' Look up permits
#'
#'
#'
#' @param polys is a spatial polygons object (sampling boxes)
#' @param lots is a spatial polygons object (taxlots)
#' @param permits is a csv file with cols(maptaxlot, record_no)
#' @return a character vector listing permit record numbers for taxlots in polys


lookup_permits <- function(polys, lots, permits){
  perms <- vector(length = length(polys), mode = 'character')
  box <- match_crs(samples, raster::crs(lots))
  for (i in seq_along(samples)) {
    lot <- crop(lots, box[i,])
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
#' @import magrittr
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
  polys <- sp::spTransform(polys, raster::crs(crs_ref))
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
#' @param csv is a csv with cols (sample_id, year, type, results(1:50))
#' @param year1 is a number matching csv$year
#' @param year2 is a number matching csv$year
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
                        heading='% Cover Change',
                        leg_pos='topleft',
                        y_lab = 'No. of Sample Sites with Cover Change of X')  {
  dt <- setDT(csv)
  year1 <- csv[year == year1]
  year2 <- csv[year == year2]
  chng <- rep(NA,nrow(year2))
  for (i in 1:nrow(year2))  {
    try(chng[i] <- year2$mn[i] - year1$mn[year1$sample_id == year2$sample_id[i]], silent=T)
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
#' @param csv is a csv with cols (sample_id, year, type, results(1:50))
#' @param title is the character vector to name the output file (png)
#' @param heading is the character vector used as the main title of the plot
#' @param leg_pos is the legend position argument (character)
#' @param y_lab is the label for the y axis (character)
#' @return prints a graph of mean cover extent into the working directory
#' @import data.table
#' @export
plot_cover <- function(csv,
                       title='cover_change.png',
                       heading='% Cover Change',
                       leg_pos='topleft',
                       y_lab = 'Proportion of Samples at or below X% Cover')  {
  dt <- setDT(csv)
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



#' plot predicted cover
#' 
#' Predicts cover extent within a polygon and returns the results as a raster plot
#' 
#' @param poly is a spatial polygons object (taxlot or sampling box)
#' @param path is a character string indicating the path to a directory to rasters
#' @param title is the name of the png plot printed to the working directory
#' @return a raster plot of predicted cover masked by `poly`
#' @export

plot_pred_cover <- function(poly, path, title = 'pred_change.png')  {
  files <- get_rasters(path)
  og_dir <- getwd()
  setwd(path)
  nir <- fill_extent(poly, path, 4)
  red <- fill_extent(poly, path, 1)
  grn <- fill_extent(poly, path, 2)
  blu <- fill_extent(poly, path, 3)
  ndvi <- (nir - red) / (nir + red)
  dt <- setDT(
    data.frame(ndvi = raster::values(ndvi),
               red = raster::values(red),
               grn = raster::values(grn),
               blue = raster::values(blu)))
  pred <- predict(mod, newdata = dt)
  ras <- ndvi
  raster::values(ras) <- pred
  ras <- raster::mask(ras, poly)
  png(title, width=9, height=6, units='in', res=300)
  plot(ras)
  dev.off()
  
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
  files <- dir()
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



#' SNAP POINT TO STREAM FLOWLINE
#'
#'
#'
#' @param pt is a coordinate pair (numeric vector)
#' @param mat is a numeric matrix representing the stream path
#' @return the point on the stream path closest to \code{pt}



streamsnap <- function(pt,mat)	{
  ptdif <- raster::pointDistance(pt,mat,lonlat=F)
  nebr <- mat[ptdif==min(ptdif),]
  nebrid <- ((ptdif==min(ptdif)) %>% as.numeric * 1:nrow(mat)) %>% sum
  nebrdif <- raster::pointDistance(nebr,mat,lonlat=F)
  above <- FALSE
  if (nebrid < nrow(mat))	{
    if (ptdif[nebrid+1]<nebrdif[nebrid+1]) above <- TRUE
  }
  if(above) new <- mat[1:nebrid,] %>% rbind(pt %>% rbind(mat[(nebrid+1):nrow(mat),]))
  if(!above) new <- mat[1:(nebrid-1),] %>% rbind(pt %>% rbind(mat[nebrid:nrow(mat),]))

  return(new)
}



#triangle length

tri_length <- function(pt,bear,dist,north=TRUE)	{
  #given hypotenuse length of right triangle
  #and coords of start point
  #returns coords of end point along hypotenuse

  if(bear>=90 & bear<270)	{
    C <- abs(bear-180)
    north <- FALSE
  }

  if(bear<90) C <- bear
  if(bear>=270) C <- 360-bear
  a <- (dist*sin((90-C)*pi/180))/sin(90*pi/180)
  if(north) y <- pt[2] + a
  if(!north) y <- pt[2] - a

  c <- (dist*sin(C*pi/180))/sin(90*pi/180)
  if(bear<=180) x <- pt[1] + c
  if(bear>180) x <- pt[1] - c

  pt <- c(x,y)
  return(pt)
}

get_bear <- function(coords)	{
  bear <- atan2(coords[2,1]-coords[1,1],coords[2,2]-coords[1,2])
  if(bear<0) bear <- 2*pi + bear
  bear <- bear * (180/pi)
  if(bear>360) bear <- bear - 360
  return(bear)
}

set_flag <- function(pt,mat,dist=75%>%FtoM,us=F,k=0)	{
  #given a matrix of coords from start point past end point
  #returns end point along path at distance dist

  mat <- pt %>% streamsnap(mat)
  difs <- pt %>% pointDistance(mat,lonlat=F)
  for (i in 1:nrow(mat))	{
    if (difs[i]==min(difs)) ptid <- i
  }

  if (us) mat <- mat[ptid:nrow(mat),]
  if (!us) mat <- mat[ptid:1,]

  so_far <- 0
  while (so_far <= dist)	{
    k <- k + 1
    rem <- dist - so_far
    so_far <- so_far + norm_vec(mat[k+1,]-mat[k,])
  }
  if (nrow(mat)>=(k+1))	{
    bear <- mat[k:(k+1),] %>% get_bear
  }
  if (nrow(mat)==k)	{
    bear <- mat[(k-1):k,] %>% get_bear
  }

  pt <- mat[k,] %>% tri_length(bear,rem)
  return(pt)
}



logindex <- function(log,ids=0)	{
  #given vector of logical operators
  #returns vector of index ids of true among false
  if (log %>% dim %>% length < 2)	{
    for (i in 1:length(log))	{
      if(log[i]) ids <- c(ids,i)
    }
  }
  if (log %>% dim %>% length > 1)	{
    for ( i in 1:nrow(log))	{
      if(log[i,1] & log[i,2]) ids <- c(ids,i)
    }
  }

  ids <- ids[-1]
  return(ids)
}

set_pt <- function(pt,mat,dist,k=0)	{
  #given a matrix of coords from start point past end point
  #returns end point along path at distance dist

  mat <- pt %>% snap_pt(mat)
  difs <- pt %>% pointDistance(mat,lonlat=F)
  for (i in 1:nrow(mat))	{
    if (difs[i]==min(difs)) ptid <- i
  }

  mat <- mat[ptid:nrow(mat),]

  so_far <- 0
  while (so_far <= dist)	{
    k <- k + 1
    rem <- dist - so_far
    so_far <- so_far + norm_vec(mat[k+1,]-mat[k,])
  }
  if (nrow(mat)>=(k+1))	{
    bear <- mat[k:(k+1),] %>% get_bear
  }
  if (nrow(mat)==k)	{
    bear <- mat[(k-1):k,] %>% get_bear
  }

  pt <- mat[k,] %>% tri_length(bear,rem)
  return(pt)
}


snap_pt <- function(pt,mat)	{
  ptdif <- pointDistance(pt,mat,lonlat=F)
  nebr <- mat[ptdif==min(ptdif),]
  nebrid <- (ptdif==min(ptdif)) %>% logindex
  nebrid <- nebrid[1]
  nebrdif <- pointDistance(nebr,mat,lonlat=F)
  above <- FALSE
  if (nebrid < nrow(mat))	{
    if (ptdif[nebrid+1]<nebrdif[nebrid+1]) above <- TRUE
  }
  if(above) new <- mat[1:nebrid,] %>% rbind(pt %>% rbind(mat[(nebrid+1):nrow(mat),]))
  if(!above) new <- mat[1:(nebrid-1),] %>% rbind(pt %>% rbind(mat[nebrid:nrow(mat),]))

  return(new)
}







#' Print Sample Plots
#'
#' Given a shapefile of sampling boxes, and a path to a directory of 4-band ortho-imagery
#' prints rgb and/or ndvi plots of the sampling area with numbered boxes superimposed.
#' While the function prints both rgb and ndvi by default, the user can turn either off
#' by setting either the \code{rgb} or \code{ndvi} arguments to \code{FALSE}.
#'
#' @param samples is a SpatialPolygonsDataFrame object (the sample boxes shapefile)
#' @param in_path is a character vector containing the path to the orthoimagery
#' @param out_path is a character vector specifying the output directory for the plots
#' @param rgb is a logical indicator, prints rgb plots if true
#' @param ndvi is a logical indicator, prints ndvi plots if true
#' @return prints an rgb and/or ndvi plot into the \code{out_path} directory
#' @export


plot_samples <- function(samples, in_path, out_path, rgb_test=TRUE, ndvi_test=TRUE)  {
  og_wd <- getwd()
  setwd(in_path)
  files <- dir()
  crs_ref <- raster::crs(raster::raster(files[1]))
  samples <- match_crs(list(samples), crs_ref)
  samples <- samples[[1]]

  for (i in seq_along(samples))	{
    setwd(in_path)
    frame <- raster::extent(samples[i,])
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
      if (rgb_test) rgb <- raster::brick(files[hit[1]])
      if (ndvi_test) {
        nir <- raster::raster(files[hit[1]], band = 4)
        red <- raster::raster(files[hit[1]], band = 1)
      }
    }

    if (length(hit) > 1)	{
      for (j in 2:length(hit))	{
        if (rgb_test) rgb <- raster::mosaic(rgb, raster::brick(files[hit[j]]), fun = max)
        if (ndvi_test)  {
          nir <- raster::mosaic(nir, raster::raster(files[hit[j]], band = 4), fun = max)
          red <- raster::mosaic(red, raster::raster(files[hit[j]], band = 1), fun = max)
        }
      }
    }

    if (rgb_test) rgb <- raster::crop(rgb, frame)
    if (ndvi_test)  {
      nir <- raster::crop(nir, frame)
      red <- raster::crop(red, frame)
      ndvi <- (nir - red) / (nir + red)
    }

    setwd(out_path)

    if (rgb_test)  {
      png(paste0('rgb_',i,'.png'))
      raster::plotRGB(rgb, main = paste0('Sample Site ',i))
      area <- lapply(methods::slot(samples[i,], 'polygons'),
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

    if (ndvi_test)  {
      png(paste0('ndvi_',i,'.png'))
      plot(ndvi, main = paste0('Sample Site ',i))
      area <- lapply(methods::slot(samples[i,], 'polygons'),
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













