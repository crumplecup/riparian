





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
  files <- get_ras(chng_path)
  for (i in seq_along(files))  {
    setwd(chng_path)
    ras <- raster::raster(file.path(chng_path, files[i]))
    crs1 <- get_crs(year1_path)
    ras <- raster::projectRaster(ras, crs = crs1)
    box <- ext_box(ras)
    r <- fill_extent(box, year1_path, 1)
    g <- fill_extent(box, year1_path, 2)
    b <- fill_extent(box, year1_path, 3)
    year1 <- raster::stack(r, g, b)
    buff1 <- sf::st_transform(buff, crs1)
    buff1 <- raster::crop(buff1, box)
    lot1 <- sf::st_transform(lots, crs1)
    lot1 <- raster::crop(lot1, box)
    crs2 <- get_crs(year2_path)
    ras <- raster::projectRaster(ras, crs = crs2)
    box <- ext_box(ras)
    r <- fill_extent(box, year2_path, 1)
    g <- fill_extent(box, year2_path, 2)
    b <- fill_extent(box, year2_path, 3)
    year2 <- raster::stack(r, g, b)
    buff2 <- sf::st_transform(buff, crs2)
    buff2 <- raster::crop(buff2, box)
    lot2 <- sf::st_transform(lots, crs2)
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


#' Look up MapTaxlot`
#'
#' Given a list of polygons, returns character vector of all maptaxlot numbers
#' within each polygon.
#'
#' @param so is a spatial object (sampling boxes)
#' @param lots is an spdf of polygons (bc taxlots)
#' @return the MapTaxlot numbers of each input polygon in a character vector
#' @seealso build_change_table pred_cover_report
#' @export
lookup_lots <- function(so, lots)  {
  lot <- suppressWarnings(sf::st_crop(lots, so))
  mtls <- squash(levels(factor(lot$MapTaxlot)))
  mtls
}


#' Look_up_permits
#'
#'
#'
#' @param so is a spatial polygons object (sampling boxes)
#' @param lots is a spatial polygons object (taxlots)
#' @param permits is a csv file with cols(maptaxlot, record_no)
#' @return a character vector listing permit record numbers for taxlots in polys
#' @seealso build_change_table
lookup_permits <- function(so, lots, permits){
  names(permits[,1]) <- c('record_no')
  names(permits[,5]) <- c('maptaxlot')
  perms <- vector(length = length(so), mode = 'character')
  crs_ref <- sf::st_crs(so)
  polys <- sf::st_transform(lots, crs_ref)
  for (i in seq_along(so)) {
    if (class(so)[1] %in% c('sf')) {
      lot <- suppressWarnings(sf::st_crop(polys, so[i, ]))
      
    }
    if (class(so)[1] %in% c('SpatialPolygons', 'SpatialPolygonsDataFrame')) {
      lot <- raster::crop(polys, matrix(crs_ref, nrow = 2))
    }
    if (class(so)[1] %in% c('RasterLayer', 'RasterStack', 'RasterBrick'))  {
      lot <- raster::crop(polys, so[i,], mask = T)
    }
    map_nos <- levels(factor(lot$MapTaxlot))
    perms[i] <- squash(permits$record_no[permits$maptaxlot %in% map_nos])
  }
  perms
}



#' match_crs
#'
#' Make sure your extents overlap by converting polygons
#' into a common reference crs.
#'
#' @param so is a spatial polygons objects
#' @param target is a spatial object with the target crs projection
#' @return spatial object projected in the target crs
match_crs <- function(so, target)  {
  crs_ref <- sf::st_crs(target)
  if (class(so)[1] %in% c('sf', 'sfc', 'sfg'))  {
    so <- sf::st_transform(so, crs_ref)
  }
  if (class(so)[1] %in% c('RasterLayer', 'RasterStack'))  {
    so <- raster::projectRaster(so, crs_ref[[1]])
  }
  if (class(so)[1] %in% c('SpatialPolygons', 'SpatialPolygonsDataFrame'))  {
    so <- sp::spTransform(so, crs_ref[[1]])
  }
  so
}


#' sample box
#' 
#' Given points for the centerline, left and right bank, assembles the points
#' into 50 slices, and the slices into a spatial Polgyons object.
#' 
#' @param mid is the centerline matrix output by centerline
#' @param left is the left bank matrix output by draw_bank
#' @param right is the right bank matrix output by draw_bank
#' @param id is a unique integer identifier for the sampling box
#' @return a Polygons object (sampling box) containing 50 Polygon objects (slices)
#' @seealso centerline
#' @seealso draw_bank
#' @export
sample_box <- function(mid, left, right, id)  {
  lbox <- list()
  rbox <- list()
  box <- list()
  k <- 1
  for (i in 1:(nrow(mid)-1))	{
    rbox[[i]] <-  sp::Polygon(rbind(right[i:(i+1),], rbind(mid[(i+1):i,], right[i,])))
    lbox[[i]] <- sp::Polygon(rbind(left[i:(i+1),], rbind(mid[(i+1):i,], left[i,])))
    box[[k]] <- rbox[[i]]
    k <- k + 1
    box[[k]] <- lbox[[i]]
    k <- k + 1
  }
  sp::Polygons(box, ID=id)
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
  files <- get_rasters(dir)
  poly <- match_crs(poly, raster::raster(file.path(dir, files[1])))
  ras <- fill_extent(poly, dir, band)
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
  coords <- sf::st_coordinates(sf::st_geometry(poly))
  unlist(
    lapply(coords, function(a, b, c) read_band(sf::st_multipolygon(a), b, c),
           b = dir, c = band))
}



#' set point
#' 
#' Given a starting point `pt` and a line `mat`,
#' returns a point distance `dist` from point `pt` along line `mat`
#' 
#' @param pt is a coordinate pair (starting point)
#' @param mat is a matrix representing a line path with cols (x,y)
#' @param dist is a distance in meters
#' @return the point distance `dist` along line `mat` from starting point `pt`
#' @export

set_pt <- function(pt,mat,dist,k=0)	{
  mat <- snap_pt(pt, mat)
  difs <- raster::pointDistance(pt, mat, lonlat=F)
  for (i in 1:nrow(mat))	{
    if (difs[i] == min(difs)) ptid <- i
  }
  
  mat <- mat[ptid:nrow(mat), ]
  
  so_far <- 0
  while (so_far <= dist)	{
    k <- k + 1
    rem <- dist - so_far
    so_far <- so_far + norm(matrix(mat[k+1,] - mat[k,], ncol = 2), type = '2')
  }
  if (nrow(mat)>=(k+1))	{
    bear <- get_bear(mat[k:(k+1),])
  }
  if (nrow(mat)==k)	{
    bear <- get_bear(mat[(k-1):k,])
  }
  
  pt <- tri_length(mat[k,], bear, rem)
  return(pt)
}


#' snap point
#' 
#' given a point along a stream path, and a matrix of points representing the stream path,
#' adds the point to the matrix of points and returns a new matrix
#' 
#' @param pt is a coordinate pair
#' @param mat is a matrix with cols (x,y) and rows of points
#' @return matrix `mat` with `pt` inserted by nearest neighbor
#' @importFrom magrittr %>%
#' @export

snap_pt <- function(pt, mat)	{
  ptdif <- raster::pointDistance(pt, mat, lonlat=F)
  nebr <- mat[ptdif == min(ptdif), ]
  nebrid <- (ptdif == min(ptdif)) %>% logindex
  nebrid <- nebrid[1]
  nebrdif <- raster::pointDistance(nebr,mat,lonlat=F)
  above <- FALSE
  if (nebrid < nrow(mat))	{
    if (ptdif[nebrid+1] < nebrdif[nebrid+1]) above <- TRUE
  }
  if(above) new <- mat[1:nebrid,] %>% rbind(pt %>% rbind(mat[(nebrid+1):nrow(mat),]))
  if(!above) new <- mat[1:(nebrid-1),] %>% rbind(pt %>% rbind(mat[nebrid:nrow(mat),]))
  
  return(new)
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




#' draw bank
#' 
#' Takes output from get_buffer and draws the left and right bank points for
#' sampling box.
#' 
#' @param mat is a buffer matrix output by get_buffer
#' @param box is a vector of corner row ids output by get_buffer
#' @param bank is 'L' for left bank and 'R' for right bank
#' @return right or left bank segment of sampling box points as matrix
#' @seealso get_buffer
#' @export
draw_bank <- function(mat, box, bank='L', crs)  {
  mat <- rbind(mat, mat)
  above <- TRUE
  done <- FALSE
  if (bank=='L')  {
    start <- box[2]
    end <- box[1]
    off <- box[4]
  }
  if (bank=='R')  {
    start <- box[4]
    end <- box[3]
    off <- box[2]
  }
  
  k <- start
  while(above & !done)	{
    k <- k + 1
    if (mat[end,1]==mat[k,1] & mat[end,2]==mat[k,2])	{
      seg <- mat[start:k,]
      done <- TRUE
    }
    if (mat[off,1]==mat[k,1] & mat[off,2]==mat[k,2]) above <- FALSE
  }
  k <- start + nrow(mat)*0.5
  while(!above & !done)	{
    k <- k - 1
    if (mat[end,1]==mat[k,1] & mat[end,2]==mat[k,2])	{
      seg <- mat[(start+nrow(mat)*0.5):k,]
      done <- TRUE
    }
  }
  seglen <- sp::Line(seg) %>% list %>% sp::Lines(ID='first') %>% list %>% sp::SpatialLines(crs)
  seglen <- rgeos::gLength(seglen) / 25
  
  side <- mat[start,] %>% as.vector %>% as.matrix %>% t
  for (j in 1:24)	{
    side <- side %>% rbind(side[nrow(side), ] %>% set_pt(mat=seg, dist=seglen))
  }
  side <- side %>% rbind(mat[end,])
  side
}

#   
#   #find length of 'right' side
#   above <- TRUE
#   done <- FALSE
#   k <- ruid
#   while(above)	{
#     k <- k + 1
#     if (right_down[1] == buff[k, 1] & right_down[2] == buff[k, 2])	{
#       seg <- buff[ruid:k, ]
#       above <- FALSE
#       done <- TRUE
#     }
#     if (left_up[1] == buff[k,1] & left_up[2] == buff[k,2]) above <- FALSE
#   }
#   k <- ruid + nrow(buff)*0.5
#   while(!above & !done)	{
#     k <- k - 1
#     if (right_down[1] == buff[k,1] & right_down[2] == buff[k,2])	{
#       seg <- buff[(ruid + nrow(buff) * 0.5):k, ]
#       above <- TRUE
#     }
#   }
#   seglen <- seg %>% Line %>% list %>% Lines(ID='first') %>% list %>% SpatialLines(crs_ref) %>% lineLength / 25
#   
#   #record 'right' side points
#   right_side <- right_up %>% as.vector %>% as.matrix %>% t
#   for (j in 1:24)	{
#     right_side <- right_side %>% rbind(right_side[nrow(right_side), ] %>% set_pt(mat=seg, dist=seglen))
#   }
#   right_side <- right_side %>% rbind(right_down)
#   rights[[i]] <- right_side
# }


# 
# buffs %>% length
# lefts %>% length
# rights %>% length
# 
# shortseg <- list()
# shortbox <- list()
# shortbuff <- list()
# shortleft <- list()
# shortright <- list()
# 
# k <- 1
# for (i in 1:length(buffs))	{
#   if(length(lefts[[i]])>1)	{
#     shortseg[[k]] <- segmat[[i]]
#     shortbox[[k]] <- boxes[[i]]
#     shortbuff[[k]] <- buffs[[i]]
#     shortleft[[k]] <- lefts[[i]]
#     shortright[[k]] <- rights[[i]]
#     k <- k+1
#   }
# }
# shortseg %>% length
# 
# segmat <- shortseg
# boxes <- shortbox
# buffs <- shortbuff
# lefts <- shortleft
# rights <- shortright
# 
# 


# polys <- list()
# right_ids <- paste0('box_',c(replicate(5,0),replicate(5,1),replicate(5,2),replicate(5,3),
#                              replicate(5,4)), replicate(5,1:5))
# left_ids <- paste0('box_',c(replicate(5,9),replicate(5,8),replicate(5,7),replicate(5,6),
#                             replicate(5,5)), replicate(5,5:1))
# 
# for (i in 1:length(boxes))	{
#   if (boxes[[i]] %>% length > 1)	{
#     mid <- boxes[[i]]
#     left <- lefts[[i]]
#     right <- rights[[i]]



#' get buffer
#' 
#' Given the matrix output from the centerline function, generates the riparian
#' buffer, adds the four corners of the sampling box explicitly to the buffer matrix,
#' and returns the buffer and a vector of corner row ids in a list.
#' 
#' @param mat is a centerline matrix
#' @return a list containing the buffer matrix and vector of box corner row ids
#' @seealso centerline
#' @export
get_buffer <- function(mat, crs) {
  bear <- get_bear(mat[1:2,]) + 90
  if (bear > 360) bear <- bear - 360
  left_up <-  tri_length(mat[1,], bear, dist=50*0.3048)
  bear <- bear + 180
  if (bear > 360) bear <- bear - 360
  right_up <- tri_length(mat[1,], bear, dist=50*0.3048)
  bear <- get_bear(mat[(nrow(mat) - 1):nrow(mat), ]) + 90
  if (bear > 360) bear <- bear - 360
  left_down <- tri_length(mat[nrow(mat), ], bear, dist=50*0.3048)
  bear <- bear + 180
  if (bear > 360) bear <- bear - 360
  right_down <- tri_length(mat[nrow(mat), ], bear, dist=50*0.3048)
  
  #create buffer polygon and pull matrix of coordinates
  buff <- sp::Line(mat) %>% list %>% sp::Lines(ID='first') %>% list %>% sp::SpatialLines(crs)
  buff <- rgeos::gBuffer(buff, width=50*0.3048)
  buff <- coord_lis(buff)[[1]]
  
  #add corner points
  # streamsnap finds nearest point on buffer to sampling box corner
  buff <- snap_pt(left_up, buff)
  buff <- snap_pt(right_up, buff)
  buff <- snap_pt(left_down, buff)
  buff <- snap_pt(right_down, buff)
  
  #save corner indexes
  for (i in 1:nrow(buff))	{
    if (left_up[1]==buff[i,1] & left_up[2]==buff[i,2]) luid <- i
    if (right_up[1]==buff[i,1] & right_up[2]==buff[i,2]) ruid <- i
    if (left_down[1]==buff[i,1] & left_down[2]==buff[i,2]) ldid <- i
    if (right_down[1]==buff[i,1] & right_down[2]==buff[i,2]) rdid <- i
  }
  corners <- c(ldid, luid, rdid, ruid)
  names(corners) <- c('ldid', 'luid', 'rdid', 'ruid')
  return(list(buff, corners))
}


#' pt to line
#' 
#' returns nearest point on line to `pt` among vertices in the line path
#' 
#' @param pt is a coordinate pair
#' @param line is SpatialLines object
#' @return nearest point on line to `pt` among vertices in the line path
#' @export
pt_to_line <- function(pt, line)  {
  difs <- apply(line, 1, function(a) raster::pointDistance(a, pt, lonlat = F))
  nebr <- logindex(difs == min(difs))
  line[nebr,]
}




#' box from extent 1
#' 
#' produce polygon representation of extent
#' 
#' @param ras is a spatial object (raster)
#' @return a Spatial Polygon drawn around the extent of `ras`
# ext_box <- function(ras, crs_ref)  {
#   ext <- raster::extent(ras)
#   crs_ref <- raster::crs(ras)
#   box <- matrix(ncol = 2, nrow = 5)
#   box[1,] <- ext[c(1,3)]
#   box[2,] <- ext[c(1,4)]
#   box[3,] <- ext[c(2,4)]
#   box[4,] <- ext[c(2,3)]
#   box[5,] <- ext[c(1,3)]
#   spatialize(box, crs_ref)
# }

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


# will_fail <- function() {
#   stop("Fail and bail!")
# }
# 
# flag <- tryCatch(will_fail(), error=function(cond) {
#   message(cond)
#   return(NA)
# } )

#' poly_to_poly
#' 
#' extracts list of polygons from each spatial object
#' add polygons to master list
#' turn master list into new SpatialPolygonDataFrame object
#' 
#' @param poly_list is a list of SpatialPolygonDataFrame objects
#' @return a SpatialPolygonDataFrame with elements from `poly_list`
#' @export

poly_to_poly <- function(poly_list) {
  crs_ref <- raster::crs(poly_list[[1]])
  ss <- list()
  id <- 1
  
  for (i in seq_along(poly_list))  {
    so <- poly_list[[i]]
    if (i > 1) so <- sp::spTransform(so, crs_ref)
    polys <- methods::slot(so, 'polygons')
    
    
    for (j in seq_along(polys)) {
      methods::slot(polys[[j]], 'ID') <- as.character(id)
      ss[[id]] <- polys[[j]]
      id <- id + 1
    }
  }
  
  ids <- data.frame(ID = as.character(1:length(ss)))
  rownames(ids) <- as.character(1:length(ss))
  
  so <- sp::SpatialPolygons(ss, proj4string = crs_ref)
  so <- sp::SpatialPolygonsDataFrame(so, data = ids)
  so
}


#' thumbnails 1
#' 
#' gets list of raster names from `in_path`
#' searches for rasters in extent of `polys`
#' mosaics rasters in extent together and crops to extent
#' saves raster in `out_path`
#' 
#' @param in_path is a character vector containing the path to the orthoimagery
#' @param out_path is a character vector specifying the output directory for the plots
#' @param polys is a SpatialPolygonsDataFrame object (the sample boxes shapefile)
#' @return saves rasters of extent `polys` to `out_path`
#' @importFrom magrittr %>%
#' @export

# 
# thumbnails <- function(in_path, out_path, polys = samples)  {
#   files <- get_ras(in_path)
#   crs_ref <- raster::crs(raster::raster(file.path(in_path, files[1])))
#   polys <- sp::spTransform(polys, crs_ref)
# 
#   for (i in seq_along(polys))	{
#     frame <- raster::extent(polys[i,])
#     hit <- 0
# 
#     for (j in seq_along(files))	{
#       if (riparian::in_extent(frame, raster::raster(file.path(in_path, files[j]))))	{
#         hit <- c(hit,j)
#       }
#     }
#     hit <- hit[-1]
# 
#     if (length(hit) == 0)  {
#       print(paste0('Failed to find raster in extent for sample ',i))
#     }
# 
#     if (length(hit) > 0)  {
#       ras <- raster::stack(file.path(in_path, files[hit[1]]))
#       ras <- raster::crop(ras, frame)
#     }
# 
#     if (length(hit) > 1)	{
#       for (k in 2:length(hit))	{
#         r <- raster::stack(file.path(in_path, files[hit[k]]))
#         r <- raster::crop(r, frame)
#         ras <- raster::mosaic(ras, r, fun = max)
#       }
#     }
# 
#     ras <- raster::crop(ras, frame)
#     raster::writeRaster(ras, filename = file.path(out_path, paste0('sample_', i)),
#                         format = 'GTiff')
#   }
# }



#' Print Sample Plots
#'
#' Superimposes sampling boxes over orthographic data.
#' Predicts the level of cover using specified `method`.
#' `method` options are c('lm', 'binom', 'lm3').
#' Default is 'lm', use 'lm3' for 3-band data.
#'
#' @param in_path is a character vector containing the path to the orthoimagery
#' @param out_path is a character vector specifying the output directory for the plots
#' @param polys is a SpatialPolygonsDataFrame object (the sample boxes shapefile)
#' @param year is an integer specifying the year of orthographic survey
#' @param type is a character vector length 1 specifying sample type
#' @param method is a character vector specifying model type
#' @return prints rgb plots of `polys` to `out_path` & a csv of predicted results
#' @importFrom magrittr %>%
#' @export


plot_samples <- function(in_path, 
                         out_path, 
                         polys = samples,
                         year,
                         type = 'random',
                         method = 'lm')  {
  files <- get_ras(in_path)
  crs_ref <- raster::crs(raster::raster(file.path(in_path, files[1])))
  polys <- sp::spTransform(polys, crs_ref)
  if (method == 'ml') model <- keras::load_model_hdf5(mod_path)
  obs <- array(0, c(length(polys), 50))
  ob <- 0
  
  for (i in seq_along(polys))	{
    frame <- raster::extent(polys[i,])
    hit <- 0
    
    for (j in seq_along(files))	{
      if (riparian::in_extent(frame, raster::raster(file.path(in_path, files[j]))))	{
        hit <- c(hit,j)
      }
    }
    hit <- hit[-1]
    
    if (length(hit) == 0)  {
      print(paste0('Failed to find raster in extent for sample ',i))
    }
    
    if (length(hit) > 0)  {
      sam <- raster::stack(file.path(in_path, files[hit[1]]))
      sam <- raster::crop(sam, frame)
    }
    
    if (length(hit) > 1)	{
      for (k in 2:length(hit))	{
        r <- raster::stack(file.path(in_path, files[hit[k]]))
        r <- raster::crop(r, frame)
        sam <- raster::mosaic(sam, r, fun = max)
      }
    }
    
    sam <- raster::crop(sam, frame)
    
    
    png(file.path(out_path, paste0('sample_',i,'.png')))
    raster::plotRGB(sam, r = 1, g = 2, b = 3, main = paste0('Sample Site ',i))
    area <- lapply(methods::slot(polys[i,], 'polygons'),
                   function(x) lapply(methods::slot(x, 'Polygons'),
                                      function(y) methods::slot(y, 'coords')))
    area <- area[[1]]
    
    for (j in 1:50)  {
      box <- riparian::spatialize(area[j], crs_ref)
      frm <- raster::extent(box)
      slc <- raster::mask(sam, box)
      slc <- raster::crop(slc, frm)
      
      if (method == 'lm') {
        car <- color_array(slc)
        prd <- predict(rip_lmod, newdata = car)
        cov <- 0
        if (prd > .25) cov <- 1
        if (prd > .75) cov <- 2
        obs[i, j] <- cov
      }
      
      if (method == 'lm3') {
        car <- color_array_3band(slc)
        prd <- predict(rip_lmod3, newdata = car)
        cov <- 0
        if (prd > .25) cov <- 1
        if (prd > .75) cov <- 2
        obs[i, j] <- cov
      }
      
      
      if (method == 'binom') {
        car <- color_array(slc)
        cov <- predict(rip_bin1, newdata = car)
        if (cov > 0) cov <- cov + predict(rip_bin2, newdata = car)
        obs[i, j] <- cov
      }
      
      
      
      
      
      pal <- get_palette(c('crimson', 'gold', 'forest'), .7)
      lines(area[[j]], col = get_palette('slate', .5), lwd = 3)
      lines(area[[j]][1:2, ], col = pal[cov+1], lwd = 7)
    }
    lines(area[[1]], lty = 2, lwd = 2)
    labcords <- matrix(0, nrow = 50, ncol = 2)
    for (j in 1:50)	{
      labcords[j,] <- methods::slot(
        rgeos::gCentroid(
          riparian::spatialize(area[[j]], crs_ref)
        ), 'coords')
    }
    text(labcords, as.character(1:50))
    dev.off()
    
    
  }
  
  df <- as.data.frame(obs)
  colnames(df) <- 1:50
  df$id <- 1:nrow(df)
  df$year <- year
  df$type <- type
  df <- df[ , c(51:53,1:50)]
  write.csv(df, file = file.path(out_path, paste0('samples_', year, '.csv')))
  df
}



#' set flag 1
#' 
#' 
#' 
#' @param pt is a coordinate pair (starting point)
#' @param mat is a matrix of coords representing a line with cols (x,y)
#' @param dist is a distance in meters
#' @param us is a logical boolean indicating whether distance is upstream or downstream
#' @return the point distance `dist` from starting point `pt` along line `mat`
#' @importFrom magrittr %>%
#' @export
# 
# set_flag <- function(pt,mat,dist=75*0.3048,us=F)	{
#   mat <- streamsnap(pt, mat)
#   difs <- raster::pointDistance(pt, mat, lonlat=F)
#   for (i in 1:nrow(mat))	{
#     if (difs[i] == min(difs)) ptid <- i
#   }
#   
#   if (us) mat <- mat[ptid:nrow(mat), ]
#   if (!us) mat <- mat[ptid:1, ]
#   
#   so_far <- 0
#   k <- 0
#   while (so_far <= dist)	{
#     k <- k + 1
#     rem <- dist - so_far
#     so_far <- so_far + norm(matrix(mat[k+1,] - mat[k,], ncol = 2), type = '2')
#   }
#   if (nrow(mat) >= (k+1))	{
#     bear <- get_bear(mat[k:(k+1),])
#   }
#   if (nrow(mat) == k)	{
#     bear <- get_bear(mat[(k-1):k,])
#   }
#   
#   pt <- tri_length(mat[k,], bear, rem)
#   return(pt)
# }




#' snap point to stream 1
#'
#' Given a coordinate pair `pt` and a matrix of coords `mat`, 
#' returns the point on `mat` nearest `pt`.
#'
#' @param pt is a coordinate pair (numeric vector)
#' @param mat is a matrix of coords representing the stream path
#' @return the point on the stream path closest to `pt`
#' @export
# 
# streamsnap <- function(pt,mat)	{
#   ptdif <- raster::pointDistance(pt,mat,lonlat=F)
#   nebr <- mat[ptdif==min(ptdif),]
#   nebrid <- ((ptdif==min(ptdif)) %>% as.numeric * 1:nrow(mat)) %>% sum
#   nebrdif <- raster::pointDistance(nebr,mat,lonlat=F)
#   above <- FALSE
#   if (nebrid < nrow(mat))	{
#     if (ptdif[nebrid+1]<nebrdif[nebrid+1]) above <- TRUE
#   }
#   if(above) new <- mat[1:nebrid,] %>% rbind(pt %>% rbind(mat[(nebrid+1):nrow(mat),]))
#   if(!above) new <- mat[1:(nebrid-1),] %>% rbind(pt %>% rbind(mat[nebrid:nrow(mat),]))
#   return(new)
# }




#' find line 1
#' 
#' finds the nearest line in a sp object, returns coordinate matrix
#' 
#' @param pt is a coordinate pair (x,y)
#' @param lines is a SpatialLinesDataFrame
#' @return coordinate matrix of nearest line to `pt` in `lines`
#' @importFrom raster pointDistance
#' @importFrom magrittr %>%
#' @export

# find_line <- function (pt, lines = prc_strms)  {
#   lines <- sp::coordinates(lines)
#   lines <- lines[[1]]
#   difs <- lapply(lines, function(a) min(apply(a, 1, function(b) raster::pointDistance(b, pt, lonlat = F))))
#   difs <- unlist(difs)
#   id <- logindex(difs == min(difs))
#   lines[[id[1]]]
# }


#' sample streams
#' 
#' Generates `n` samples along `prc`, finds the closest points on `strms`
#' to each point along `prc`, returns `n` samples along `strms`.
#' 
#' @param n is the number of samples to take
#' @param lots is a polygon object representing sample extent
#' @param prc is the line path for perennial and fish-bearing streams in RR (NHD)
#' @param strms is the hydro-enforced drainage layer in RR (WSI)
#' @return `n` points along `strms`
#' @export
sample_streams <- function(n = 100, 
                           lots = prc_mtls, 
                           prc = prc_per, 
                           strms = prc_strms,
                           buff = per_buff)  {
  crs_ref <- raster::crs(strms)
  lots <- sp::spTransform(lots, crs_ref)
  prc <- sp::spTransform(prc, crs_ref)
  buff <- sp::spTransform(buff, crs_ref)
  poly <- raster::crop(buff, lots, mask=T)
  prc <- prc[poly,]
  pts <- sp::spsample(prc, 1000, 'random')
  pts <- sample(pts, n)
  pt_cords <- sp::coordinates(pts)
  finds <- apply(pt_cords, 1, function(a) find_line(a, strms))
  if(class(finds) == 'matrix')  {
    finds <- list(matrix(finds, ncol = 2))
  }
  pts <- matrix(0, nrow = length(finds), ncol = 2)
  for (i in seq_along(finds)) {
    pts[i,] <- pt_to_line(pt_cords[i,], finds[[i]])
  }
  polys <- list()
  boxes <- centerline(pts, finds)
  for(i in 1:length(boxes))  {
    print(i)
    buff <- get_buffer(boxes[[i]], crs_ref)
    corners <- buff[[2]]
    buff <- buff[[1]]
    left_bank <- try(draw_bank(buff, corners, 'L', crs_ref), silent = TRUE)
    right_bank <- try(draw_bank(buff, corners, 'R', crs_ref), silent = TRUE)
    polys[[i]] <- try(sample_box(boxes[[i]], left_bank, right_bank, i), silent = TRUE) 
  }
  polys <- sp::SpatialPolygons(polys, proj4string=crs_ref)
  pid <- lapply(slot(polys, 'polygons'), function(x) slot(x, 'ID'))
  p.df <- data.frame( ID=1:length(polys), row.names = pid)
  sp::SpatialPolygonsDataFrame(polys, p.df)
} 




#' centerline 1
#' 
#' Given sampling points, find the centerline points of the sampling box.
#' 
#' @param pts is a spatial points object (stream samples)
#' @param strms is the stream layer `pts` lie along
#' @return a list of matrices with coordinates for 25 segments along stream path
#' @export
# centerline <- function(pts, strms) {
#   boxes <- list()
#   for(i in 1:nrow(pts))	{
#     box <- matrix(0, ncol=2)
#     try(	box <- box %>% rbind(pts[i,] %>% set_flag(strms[[i]], us=T)))
#     for (j in 1:25)	{
#       try(box <- box %>% rbind(box[nrow(box),] %>% set_flag(strms[[i]],dist= 6*0.3048)))
#     }
#     box <- box[-1,]
#     boxes[[i]] <- box
#   }
#   
#   short <- list()
#   k <- 1
#   for (i in 1:length(boxes))	{
#     if(length(boxes[[i]])>1)	{
#       short[[k]] <- boxes[[i]]
#       k <- k+1
#     }
#   }
#   short
# }
# 



#' Look up permits 1
#'
#'
#'
#' @param so is a spatial polygons object (sampling boxes)
#' @param lots is a spatial polygons object (taxlots)
#' @param permits is a csv file with cols(maptaxlot, record_no)
#' @return a character vector listing permit record numbers for taxlots in polys

# 
# lookup_permits <- function(so, lots, permits){
#   names(permits[,1]) <- c('record_no')
#   names(permits[,5]) <- c('maptaxlot')
#   perms <- vector(length = length(so), mode = 'character')
#   polys <- sp::spTransform(lots, raster::crs(so))
#   for (i in seq_along(so)) {
#     if (class(so) %in% c('SpatialPolygons', 'SpatialPolygonsDataFrame')) {
#       lot <- raster::crop(polys, raster::extent(so[i,]))
#     }
#     if (class(so) %in% c('RasterLayer', 'RasterStack', 'RasterBrick'))  {
#       lot <- raster::crop(polys, so[i,], mask = T)
#     }
#     map_nos <- levels(factor(lot$MapTaxlot))
#     perms[i] <- squash(permits$record_no[permits$maptaxlot %in% map_nos])
#   }
#   perms
# }



#' MATCH CRS
#'
#' Make sure your extents overlap by converting polygons
#' into a common reference crs.
#'
#' @param polys is a spatial polygons objects
#' @param crs_ref is a character string representing the crs
#' @return a list of spatial polygons projected in the common crs
#' @export



# match_crs <- function(polys, crs_ref)  {
#   if (class(polys) %in% c('RasterLayer', 'RasterStack'))  {
#     polys <- raster::projectRaster(polys, crs_ref)
#   }
#   if (class(polys) %in% c('SpatialPolygons', 'SpatialPolygonsDataFrame'))  {
#     polys <- sp::spTransform(polys, raster::crs(crs_ref))
#   }
#   polys
# }




#' predicted change report 1
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
# pred_change_report <- function(chng_path, 
#                                year1_path,
#                                year2_path,
#                                permits = permits_13to18,
#                                lots = prc_lots, 
#                                title = 'pred_chng_table.csv')  {
#   og_dir <- getwd()
#   files <- get_rasters(year2_path)
#   fil_len <- length(files)
#   ids <- 1:fil_len
#   
#   cov_yr1 <- array(0, fil_len)
#   area_yr1 <- array(0, fil_len)
#   cov_yr2 <- array(0, fil_len)
#   area_yr2 <- array(0, fil_len)
#   
#   chng_adj <- array(0, fil_len)
#   
#   
#   chng <- array(0, fil_len)
#   area <- vector(length = fil_len, mode = 'numeric')
#   rat <- vector(length = fil_len, mode = 'numeric')
#   mtls <- vector(length = fil_len, mode = 'character')
#   perms <- vector(length = fil_len, mode = 'character')
#   cover <- vector(length = fil_len, mode = 'numeric')
#   for (i in seq_along(files))  {
#     print(paste0('Reading ', i, ' of ', fil_len, ' files...'))
#     ras <- raster::raster(paste0(chng_path, files[i]))
#     vals <- raster::values(ras)
#     chng[i] <- sum(vals[!is.na(vals)]) / 2
#     area[i] <- length(vals[!is.na(vals)])
#     cover[i] <- (sum(vals[!is.na(vals)])/2) / area[i]
#     lot <- sp::spTransform(lots, raster::crs(ras))
#     lot <- raster::crop(lot, ras, mask = T)
#     map_nos <- levels(factor(lot$MapTaxlot))
#     mtls[i] <- squash(map_nos)
#     perms[i] <- squash(permits$record_no[permits$maptaxlot %in% map_nos])
# 
#     ras <- raster::raster(paste0(year1_path, files[i]))
#     vals <- raster::values(ras)
#     vals[vals<0] <- 0
#     vals[vals>2] <- 2
#     area_yr1[i] <- length(vals[!is.na(vals)])
#     cov_yr1[i] <- sum(vals[!is.na(vals)]) / 2
#     
#     ras <- raster::raster(paste0(year2_path, files[i]))
#     vals <- raster::values(ras)
#     vals[vals<0] <- 0
#     vals[vals>2] <- 2
#     area_yr2[i] <- length(vals[!is.na(vals)])
#     cov_yr2[i] <- sum(vals[!is.na(vals)]) / 2
#     
#   }
#   rat[area > 0] <- chng[area > 0] / area[area > 0]
#   cov_yr1[area_yr1 > 0] <- cov_yr1[area_yr1 > 0] / area_yr1[area_yr1 > 0]
#   cov_yr2[area_yr2 > 0] <- cov_yr2[area_yr2 > 0] / area_yr2[area_yr2 > 0]
#   chng_adj <- cov_yr2 - cov_yr1
#   
#   tot_yr1 <- sum(cov_yr1 * area_yr1) / sum(area_yr1)
#   tot_yr2 <- sum(cov_yr2 * area_yr2) / sum(area_yr2)
#   tot_chng <- sum(chng_adj * area_yr1) / sum(area_yr1)
#   dt <- setDT(data.frame(ids = ids, 
#                          year1 = cov_yr1, 
#                          year2 = cov_yr2,
#                          change = chng_adj,
#                          acres = area / 43560,
#                          ratio = rat,
#                          lots = mtls, 
#                          permits = perms))
#   setkey(dt, ratio)
#   setwd(og_dir)
#   write.csv(dt, title)
#   print(paste0('Predicted cover extent year 1:  ', round(tot_yr1, 2)))
#   print(paste0('Predicted cover extent year 2:  ', round(tot_yr2, 2)))
#   print(paste0('Predicted cover change:  ', round(tot_chng, 2)))
#   
#   
# }


#' Look up MapTaxlot` 1
#'
#' Given a list of polygons, returns character vector of all maptaxlot numbers
#' within each polygon.
#'
#' @param so is a spatial object (sampling boxes)
#' @param lots is an spdf of polygons (bc taxlots)
#' @return the MapTaxlot numbers of each input polygon in a character vector
#' @export


# lookup_lots <- function(so, lots)  {
#   mtls <- vector(length = length(so), mode = 'character')
#   polys <- sp::spTransform(lots, raster::crs(so))
#   for (i in seq_along(so)) {
#     lot <- raster::crop(polys, raster::extent(so[i,]))
#     mtls[i] <- squash(levels(factor(lot$MapTaxlot)))
#   }
#   mtls
# }



#' Get CRS 1
#'
#' Get a reference CRS from rasters in a directory
#'
#' @param dir is a path to a directory of rasters
#' @return the CRS character string
#' @export


# get_crs <- function(dir)  {
#   og_dir <- getwd()
#   setwd(dir)
#   files <- get_rasters(dir)
#   crs <- raster::crs(raster::raster(files[1]))
#   setwd(og_dir)
#   crs
# }
