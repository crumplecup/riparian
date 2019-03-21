

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



#' COORDINATE LIST FROM SPATIALPOLYGONS
#'
#'
#'
#' @param poly is a \code{SpatialPolygons} object
#' @return a list coordinates for vertices of polygons in \code{poly}
#' @export

coord_lis <- function(poly)  {
  coords <- lapply(slot(poly,'polygons'),
                   function(x) lapply(slot(x,'Polygons'),
                          function(y) slot(y, 'coords')))
  coords[[1]]
}



#' Cover Change Report
#'
#' Prints a graph of cover change.
#'
#' @param year1 is a matrix with cols (sample_id,year,type,results(1:50))
#' @param year2 is a matrix same format as year1
#' @return prints a graph of mean cover change into the working directory
#' @import data.table
#' @export

change_report <- function(csv,
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
  mn_chng <- mean(chng[!is.na(chng)])
  sd_chng <- sd(chng[!is.na(chng)])
  n_chng <- length(chng[!is.na(chng)])
  upr_chng <- mn_chng + 1.96 * (sd_chng / sqrt(n_chng))
  lwr_chng <- mn_chng - 1.96 * (sd_chng / sqrt(n_chng))
  pmf_chng <- density(chng[!is.na(chng)])

  png(title, width=9, height=6, units='in', res=300)
  pmf_chng %>% plot(lwd=2,col='forestgreen',
                    main=heading, ylab=y_lab, axes=F)
  axis(1, at=seq(-1,1,.1),
       labels=c('-100%','-90%','-80%', '-70%', '-60%', '-50%', '-40%', '-30%', '-20%', '-10%', '0%',
                '10%', '20%', '30%', '40%','50%','60%','70%','80%','90%','100%'))
  axis(2)
  abline(v=upr_chng,col='steelblue',lty=2)
  abline(v=lwr_chng,col='steelblue',lty=2)
  abline(v=mn_chng,col='slategrey')
  legend(leg_pos,legend=c('No. of Samples at X% Change','Mean Cover Change',
                            'Upper & Lower 95% CI on Mean Change'),
         fill = c('forestgreen','slategrey','steelblue'))
  dev.off()

}



#' MERGE RASTERS TO FILL EXTENT
#'
#' Given a \code{SpatialPolygons} object and a directory filepath character string,
#' returns merged rasters from \code{dir} over extent of \code{poly}.
#'
#' @param poly is a \code{SpatialPolygons} object
#' @param dir is a directory filepath character string
#' @param band is the band to extract from multiband imagery
#' @return merged rasters from \code{dir} over extent of \code{poly}
#' @export



fill_extent <- function(poly,dir,band)  {
  ext <- extent(poly)
  og_dir <- getwd()
  setwd(dir)
  files <- dir()
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
      ras <- raster::mosaic(ras,raster::raster(files[hit[i]], band=band),fun=max)
    }
  }
  setwd(og_dir)
  crop(ras,ext)
}




#' IN EXTENT LOGICAL TEST
#'
#' Tests whether the point (or extent) \code{pt} falls entirely within the extent of a raster \code{ras},
#' returning TRUE or FALSE.
#'
#' @param pt is a coordinate vector or extent class object
#' @param so is a spatial object
#' @param is_in is a logical binary set to FALSE
#' @return TRUE if pt is within extent of ras, else FALSE
#' @export

in_extent <- function(pt, so, is_in=FALSE)	{
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



#' Look up MapTaxLot`
#'
#' Given a polygon, returns the MapTaxlot number of the centroid.
#'
#' @param poly is a polygon (sampling box)
#' @param lots is an spdf of polygons (bc taxlots)
#' @return the MapTaxlot of the centroid of the input polygon
#' @export


lookup_lot <- function(poly, lots)  {
  spot <- rgeos::gCentroid(poly)
  for (i in seq_along(lots))  {
    if(in_extent(spot, lots[i])) {
      lots$MapTaxlot
    }
  }
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
#' Make sure your extents overlap by converting all polygons in a list
#' into a common reference crs.
#'
#' @param poly_list is a list of spatial polygon objects
#' @param crs_ref is a character string representing the crs
#' @return a list of spatial polygons projected in the common crs
#' @export



match_crs <- function(poly_list, crs_ref)  {
  lapply(poly_list, function(x) sp::spTransform(x, raster::crs(crs_ref)))
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






#' CONVERT XY COORDS TO SPATIALPOLYGONS
#'
#' Given a matrix of coords with cols (x,y), and a character string representing a coordinate
#' projection, returns a SpatialPolygon defined by coords of \code{mat} in projection \code{crs_ref}.
#'
#' @param mat is a matrix of coords with cols (x,y)
#' @param crs_ref is a character string representing a coordinate projection for proj4string
#' @return a SpatialPolygon defined by coords of mat in projection crs_ref
#' @export


spatialize <- function(mat,crs_ref)	{
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













