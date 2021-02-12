
#' bearing
#' 
#' Take bearing.
#' 
#' @param a is a pt
#' @param b is a pt
#' @return bearing in degrees
bearing <- function(a, b)	{
  bear <- atan2(b[1]-a[1],b[2]-a[2])
  if(bear<0) bear <- 2*pi + bear
  bear <- bear * (180/pi)
  if(bear>360) bear <- bear - 360
  return(bear)
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
  crs_ref <- sf::st_crs(lots)
  polys <- sf::st_transform(polys, crs_ref)
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
  mtls <- 0
  perms <- 0
  for (i in seq_along(change))  {
    change[i] <- chng[yr2$id %in% as.character(ids[i])]
    lot <- suppressWarnings(sf::st_crop(lots, polys[i, ]))
    mtls[i] <- squash(levels(factor(lot$MapTaxlot)))
    map_nos <- levels(factor(lot$MapTaxlot))
    perms[i] <- squash(permits$record_no[permits$maptaxlot %in% map_nos])
  }
  dt <- data.table::setDT(data.frame(id = ids, change = change, maptaxlot = mtls[ids], permits = perms[ids]))
  write.csv(dt[order(change)], file = paste0('change_table_', year1, '-', year2,'.csv'))
}




#' centerline
#' 
#' Given sampling points, find the centerline points of the sampling box.
#' 
#' @param pts is a spatial points object (stream samples)
#' @param strms is the stream layer `pts` lie along
#' @return a list of matrices with coordinates for 25 segments along stream path
centerline <- function(find, span, inc, strms, thresh) {
  print('calling centerline')
  if (!('sfc' %in% class(strms))) strms <- sf::st_as_sf(strms)
  strm_geom <- sf::st_geometry(strms)
  strm_cords <- sf::st_coordinates(strms)
  mat <- orient(find, strm_cords[strm_cords[,3] == find[3], ], 0)
  mat <- sf::st_segmentize(sf::st_sfc(sf::st_linestring(mat)), 2)[[1]]
  nxt_pt <- tryCatch(traverse(mat, span, strm_cords, thresh),
                     error=function(cond) {
                       message(cond)
                       return(NA)
                     })
  
  if (class(nxt_pt) == 'list') {
    ln <- nxt_pt[[1]]  # start line list
    mat <- orient(nxt_pt[[1]], 
                  strm_cords[strm_cords[,3] == nxt_pt[[1]][3], ], 
                  uturn(nxt_pt[[3]]))
    for (j in 1:25) {
      nxt_pt <- tryCatch(
        traverse(mat, inc, strm_cords, thresh),
        error=function(cond) {
          message(cond)
          return(NA)
        })
      if (is.na(nxt_pt[1])) {
        ln <- NA
        break
      } else {
        mat <- sf::st_segmentize(
          sf::st_sfc(
            sf::st_linestring(
              strm_cords[strm_cords[,3] == nxt_pt[[1]][3], ]
            )
          )
          , 2)[[1]]
        mat <- orient(nxt_pt[[1]], mat, nxt_pt[[3]])
        ln <- rbind(ln, nxt_pt[[1]])
      }
    }
  }
  ln
}


#' color_array
#'
#' extracts color statistics from raster image
#'
#' @param img is a raster stack of 4 bands rgbn
#' @return an array of summary color statistics
#' @seealso plot_samples
color_array <- function(img) {
  ar <- raster::values(img)
  ar <- ar[!is.na(ar[ ,1]), ]
  car <- array(0, c(1, 15))
  car[1, 1] <- mean(ar[ , 1])
  car[1, 2] <- sd(ar[ , 1])
  car[1, 3] <- sum(ar[ , 1])
  car[1, 4] <- mean(ar[ , 2])
  car[1, 5] <- sd(ar[ , 2])
  car[1, 6] <- sum(ar[ , 2])
  car[1, 7] <- mean(ar[ , 3])
  car[1, 8] <- sd(ar[ , 3])
  car[1, 9] <- sum(ar[ , 3])
  car[1, 10] <- mean(ar[ , 4])
  car[1, 11] <- sd(ar[ , 4])
  car[1, 12] <- sum(ar[ , 4])
  ndvi <- (ar[ , 4] - ar[ , 1]) / (ar[ , 4] + ar[ , 1])
  car[1, 13] <- mean(ndvi)
  car[1, 14] <- sd(ndvi)
  car[1, 15] <- sum(ndvi)
  df <- as.data.frame(car)
  colnames(df) <- c('red_mn', 'red_sd', 'red_sm',
                    'grn_mn', 'grn_sd', 'grn_sm',
                    'blu_mn', 'blu_sd', 'blu_sm',
                    'nir_mn', 'nir_sd', 'nir_sm',
                    'ndv_mn', 'ndv_sd', 'ndv_sm')
  df
}


#' color_array_3band
#'
#' extracts color statistics from raster image
#'
#' @param img is a raster stack of 3 bands rgb
#' @return an array of summary color statistics
#' @seealso plot_samples
color_array_3band <- function(img) {
  ar <- raster::values(img)
  ar <- ar[!is.na(ar[ ,1]), ]
  car <- array(0, c(1, 9))
  car[1, 1] <- mean(ar[ , 1])
  car[1, 2] <- sd(ar[ , 1])
  car[1, 3] <- sum(ar[ , 1])
  car[1, 4] <- mean(ar[ , 2])
  car[1, 5] <- sd(ar[ , 2])
  car[1, 6] <- sum(ar[ , 2])
  car[1, 7] <- mean(ar[ , 3])
  car[1, 8] <- sd(ar[ , 3])
  car[1, 9] <- sum(ar[ , 3])
  df <- as.data.frame(car)
  colnames(df) <- c('red_mn', 'red_sd', 'red_sm',
                    'grn_mn', 'grn_sd', 'grn_sm',
                    'blu_mn', 'blu_sd', 'blu_sm')
  df
}


#' color_array_set
#' 
#' color array set from polys for model calibration.
#' 
#' @param polys is an sf object
#' @param in_path is the directory of raster
#' @return color array data.table
#' @export
color_array_set <- function(polys, in_path) {
  for (i in seq_along(polys)) {
    poly <- polys[[i]]
    for (j in seq_along(poly)) {
      row <- get_cols_4band(poly[[j]], in_path)
      if (i == 1 & j == 1) {
        rec <- row
      } else {
        rec <- rbind(rec, row)
      }
    }
  }
  rec
}


#' comp_angle
#' 
#' Compare angles of two vectors.
#' 
#' @param bear1 is a bearing
#' @param bear2 is a bearing
#' @return difference in bearing
#' @seealso orient
comp_angle <- function(bear1, bear2) {
  if (abs(bear1 - bear2) > 180) {
    if (bear1 < bear2) {
      bear1 <- bear1 + 360
    } else {
      bear2 <- bear2 + 360
    }
  }
  return(abs(bear1 - bear2))
}


#' corners
#' 
#' Find the four corners of a sampling box.
#' 
#' @param mat is a matrix with cols c(x,y)
#' @param dist is the size of riparian corridor
#' @return a list with elements [buffer, corners]
#' @seealso sample_streams
corners <- function(mat, dist) {
  print('finding corners')
  bear <- get_bear(mat[1:2,]) + 90
  if (bear > 360) bear <- bear - 360
  left_up <-  tri_length(mat[1,], bear, dist)
  bear <- bear + 180
  if (bear > 360) bear <- bear - 360
  right_up <- tri_length(mat[1,], bear, dist)
  bear <- get_bear(mat[(nrow(mat) - 1):nrow(mat), ]) + 90
  if (bear > 360) bear <- bear - 360
  left_down <- tri_length(mat[nrow(mat), ], bear, dist)
  bear <- bear + 180
  if (bear > 360) bear <- bear - 360
  right_down <- tri_length(mat[nrow(mat), ], bear, dist)
  
  buff <- sf::st_segmentize(sf::st_buffer(sf::st_zm(sf::st_linestring(mat)), dist), .5)
  buff <- sf::st_coordinates(buff)
  corners <- rbind(left_up, rbind(left_down, rbind(right_up, right_down)))
  for (i in 1:4)  {
    difs <- raster::pointDistance(corners[i, 1:2], buff[, 1:2], lonlat = F)
    id <- sum((1:length(difs)) * (difs == min(difs)))
    corners[i, ] <- buff[id, 1:2]
  }
  list(buff, corners)
}



#' extent box
#' 
#' produce polygon representation of extent
#' 
#' @param ras is a spatial object (raster)
#' @return a Spatial Polygon drawn around the extent of `ras`
#' @seealso before_and_after
ext_box <- function(ras)  {
  ext <- raster::extent(ras)
  crs_ref <- sf::st_crs(ras)
  box <- matrix(ncol = 2, nrow = 5)
  box[1,] <- ext[c(1,3)]
  box[2,] <- ext[c(1,4)]
  box[3,] <- ext[c(2,4)]
  box[4,] <- ext[c(2,3)]
  box[5,] <- ext[c(1,3)]
  sf::st_sfc(sf::st_polygon(list(box)), crs = crs_ref)
}



#' fill_extent
#'
#' Given a `sf` object and a directory filepath character string,
#' returns merged rasters from `in_path` over extent of `poly`.
#'
#' @param poly is a `sf` object
#' @param in_path is a directory filepath character string
#' @param pad is the padding factor
#' @return merged rasters from `in_path` over extent of `poly`
#' @seealso before_and_after get_cols_4band plot_samples pred_cover1 sample_streams thumbnails
#' @export
fill_extent <- function(poly, in_path, pad = .1)  {
  files <- get_ras(in_path)
  crs_ref <- sf::st_crs(raster::raster(file.path(in_path, files[1])))
  poly <- sf::st_transform(poly, crs_ref)
  frame <- pad(poly, pad = pad)
  hit <- 0
  print('searching rasters')
  
  for (i in seq_along(files))	{
    if (
      touches(
        frame, sf::st_bbox(
          raster::raster(file.path(in_path, files[i]))
        )
      )
    )	{
      hit <- c(hit,i)
    }
  }
  hit <- hit[-1]
  
  
  if (length(hit) == 0)  {
    m <- message('Failed to find raster in extent for sample.')
    return(m)
  }
  print('raster found')
  
  if (length(hit) > 0)  {
    r <- raster::stack(file.path(in_path, files[hit[1]]))
    r <- tryCatch(raster::crop(r, matrix(frame, ncol = 2)), 
                  error = function(x) { 
                    message(x)
                    m <- message('Unable to crop raster to sample.')
                    return(m)
                  })
    sam <- r
  }
  
  if (length(hit) > 1)	{
    for (k in 2:length(hit))	{
      r <- raster::stack(file.path(in_path, files[hit[k]]))
      r <- tryCatch(raster::crop(r, matrix(frame, ncol = 2)), error = function(x) { 
        message(g)
        m <- message('Unable to crop raster to sample.')
        return(m)
      })
      sam <- raster::mosaic(sam, r, fun = max)
      print('raster mosaic created')
    }
  }
  raster::crop(sam, matrix(frame, nrow = 2))
}



#' find line
#' 
#' finds the nearest line in a sp object, returns coordinate matrix
#' 
#' @param pt is a coordinate pair (x,y)
#' @param lines is a SpatialLinesDataFrame
#' @return coordinate matrix of nearest line to `pt` in `lines`
find_line <- function (pt, lines = prc_strms)  {
  print('calling find line')
  if (!('sf' %in% class(lines))) lines <- sf::st_as_sf(lines)
  lines <- sf::st_coordinates(lines)
  difs <- apply(lines, 1, function(x) raster::pointDistance(x[1:2], pt[1:2], lonlat = F))
  id <- lines[logindex(difs == min(difs)), ]
  # list(id, st_linestring(lines[lines[,3] == id[3], ]))
  id
}


#' get bearing
#' 
#' Given a matrix of two coordinate pairs, returns the bearing of the angle
#' between the two points in degrees.
#' 
#' @param coords is a matrix with cols (x,y) and rows(1:2)
#' @return the angle of the bearing between the points in `coords` in degrees
#' @export
get_bear <- function(coords)	{
  bear <- atan2(coords[2,1]-coords[1,1],coords[2,2]-coords[1,2])
  if(bear<0) bear <- 2*pi + bear
  bear <- bear * (180/pi)
  if(bear>360) bear <- bear - 360
  return(bear)
}


#' Get colors from 4 band raster
#'
#' Returns data.table of mean ndvi and rgb values for a sampling box
#'
#' @param poly is the sampling box
#' @param dir is the path to the rasters
#' @return a dt with cols (ndvi, red, grn, blu)
get_cols_4band <- function(poly, dir)  {
  files <- get_ras(dir)
  crs_ref <- sf::st_crs(raster::raster(file.path(dir, files[1])))
  poly <- sf::st_sfc(sf::st_polygon(poly), crs = crs_ref)
  
  r <- fill_extent(poly, dir)
  r <- raster::mask(r, as(poly, 'Spatial'))
  red <- raster::values(raster::raster(r, layer = 1))
  red <- red[!is.na(red)]
  grn <- raster::values(raster::raster(r, layer = 2))
  blu <- raster::values(raster::raster(r, layer = 3))
  nir <- raster::values(raster::raster(r, layer = 4))
  ndvi <- (nir - red) / (nir + red)
  red <- red[!is.na(red)]
  grn <- grn[!is.na(grn)]
  blu <- blu[!is.na(blu)]
  nir <- nir[!is.na(nir)]
  ndvi <- ndvi[!is.na(ndvi)]
  
  red <- fitdistrplus::descdist(red, graph = F)
  red_mn <- red$mean
  red_sd <- red$sd
  red_sk <- red$skewness
  red_kr <- red$kurtosis
  
  grn <- fitdistrplus::descdist(grn, graph = F)
  grn_mn <- grn$mean
  grn_sd <- grn$sd
  grn_sk <- grn$skewness
  grn_kr <- grn$kurtosis
  
  blu <- fitdistrplus::descdist(blu, graph = F)
  blu_mn <- blu$mean
  blu_sd <- blu$sd
  blu_sk <- blu$skewness
  blu_kr <- blu$kurtosis
  
  nir <- fitdistrplus::descdist(nir, graph = F)
  nir_mn <- nir$mean
  nir_sd <- nir$sd
  nir_sk <- nir$skewness
  nir_kr <- nir$kurtosis
  
  ndvi <- fitdistrplus::descdist(ndvi, graph = F)
  ndvi_mn <- ndvi$mean
  ndvi_sd <- ndvi$sd
  ndvi_sk <- ndvi$skewness
  ndvi_kr <- ndvi$kurtosis
  
  rec <- setDT(data.frame(
    red_mn = red_mn,
    red_sd = red_sd,
    red_sk = red_sk,
    red_kr = red_kr,
    blu_mn = blu_mn,
    blu_sd = blu_sd,
    blu_sk = blu_sk,
    blu_kr = blu_kr,
    grn_mn = grn_mn,
    grn_sd = grn_sd,
    grn_sk = grn_sk,
    grn_kr = grn_kr,
    nir_mn = nir_mn,
    nir_sd = nir_sd,
    nir_sk = nir_sk,
    nir_kr = nir_kr,
    ndvi_mn = ndvi_mn,
    ndvi_sd = ndvi_sd,
    ndvi_sk = ndvi_sk,
    ndvi_kr = ndvi_kr))
  rec
}



#' Get CRS
#'
#' Get a reference CRS from rasters in a directory
#'
#' @param dir is a path to a directory of rasters
#' @return the CRS character string
get_crs <- function(dir)  {
  files <- get_ras(dir)
  crs <- sf::st_crs(raster::raster(file.path(dir, files[1])))
  crs
}



#' get_palette
#'
#' get a palette of colors at chosen transparency
#' colors are coral, hardwood, gold, forest, sky, ocean
#' violet, rose, crimson
#' choose a palette with a vector of names
#' get a random palette length `x` if `x` is a number
#'
#' @param x is a vector of color names or number
#' @param a is the alpha transparency (0,1)
#' @return a palette of rgb values
#' @export
get_palette <- function(x, a = .33) {
  col_names <- c('coral', 'hardwood', 'gold', 'forest', 'leaf', 'sky', 'ocean',
                 'violet', 'rose', 'crimson', 'white', 'slate', 'charcoal', 'black')
  cols <- c(rgb(.921, .251, .203, a, 'coral'),
            rgb(.321, .196, .129, a, 'hardwood'),
            rgb(.812, .675, .0, a, 'gold'),
            rgb(.0, .361, .024, a, 'forest'),
            rgb(.561, .82, .459, a, 'leaf'),
            rgb(.0, .753, .78, a, 'sky'),
            rgb(.02, .0, .612, a, 'ocean'),
            rgb(.608, .0, .89, a, 'violet'),
            rgb(.839, .369, .471, a, 'rose'),
            rgb(.788, .0, .0, a, 'crimson'),
            rgb(1, 1, 1, a, 'white'),
            rgb(.666, .666, .666, a, 'slate'),
            rgb(.333, .333, .333, a, 'charcoal'),
            rgb(0, 0, 0, a, 'black'))
  df <- data.frame(col_names = col_names, cols = cols,
                   stringsAsFactors = F)
  palette <- 0
  if (inherits(x, 'character')) {
    for (i in seq_along(x)) {
      palette[i] <- df$cols[df$col_names == x[i]]
    }
  }
  if (inherits(x, 'numeric')) {
    palette <- df$cols[sample.int(length(cols), x, replace = TRUE)]
  }
  return(palette)
}



#' get_ras
#'
#' get filenames of rasters in `path`
#'
#' @param path is a character vector specifying a directory of files
#' @return a character vector of filenames of rasters in `path`
#' @export
get_ras <- function(path = getwd())  {
  files <- list.files(path)
  files <- files[!grepl('.xml', files)]
  files <- files[!grepl('.ovr', files)]
  files <- files[!grepl('.ini', files)]
  tifs <- files[grepl('.tif', files)]
  hdrs <- files[grepl('.hdr', files)]
  c(tifs,hdrs)
}


#' insert_pt
#' 
#' Add point to line matrix.
#' 
#' @param pt is the point to add
#' @param mat is the line matrix
#' @return list of new matrix and id of pt on mat
insert_pt <- function(pt, mat) {
  # print("calling insert point")
  difs <- apply(mat, 1, function(x) raster::pointDistance(pt[1:2], x[1:2], lonlat = F))
  if (length(pt) < ncol(mat)) pt <- c(pt, mat[1, -c(1:length(pt))])
  if (length(pt) > ncol(mat)) pt <- pt[1:ncol(mat)]
  idx <- sum((1:length(difs)) * (difs == min(difs)))  # id of closest pt on mat
  if (idx < nrow(mat)) {  # idx is not end of mat
    # print('idx not head of mat')
    idx_nxt <- raster::pointDistance(mat[idx, 1:2], mat[idx+1, 1:2], lonlat = F)
    pt_nxt <- raster::pointDistance(pt[1:2], mat[idx+1, 1:2], lonlat = F)
    if (pt_nxt > idx_nxt) { # idx is ahead of pt
      pt_id <- idx
      if (idx == 1) {  # if first, rbind
        mat <- rbind(pt, mat)
      }
      if (idx > 1) {  # if not first, insert with two rbind
        mat <- rbind(
          mat[1:(idx-1), ], rbind(
            pt, mat[idx:nrow(mat), ]
          )
        )
      }
      
    } else { # idx is behind pt
      pt_id <- idx + 1
      mat <- rbind(
        mat[1:idx, ], rbind(
          pt, mat[(idx+1):nrow(mat), ]
        )
      )
    }
  } else {  # idx is end of mat
    # print('idx is end of mat')
    idx_nxt <- raster::pointDistance(mat[idx, 1:2], mat[idx-1, 1:2], lonlat = F)
    pt_nxt <- raster::pointDistance(pt[1:2], mat[idx-1, 1:2], lonlat = F)
    # print('idx_nxt')
    # print(idx_nxt)
    # print('pt_nxt')
    # print(pt_nxt)
    if (pt_nxt > idx_nxt) { # idx is ahead of pt
      pt_id <- idx + 1
      mat <- rbind(mat, pt)
    } else {
      pt_id <- idx
      mat <- rbind(
        mat[1:(idx-1), ], rbind(
          pt, mat[idx, ]
        )
      )
    }
  }
  list(mat, pt_id)
}



#' is_in
#' 
#' pt in box
#' 
#' @param pt is a point
#' @param box is sf bbox
#' @return boolean touches
is_in <- function(pt, box)	{
  is_in=FALSE
  if (class(box) == 'bbox')	{
    if (pt[1] >= box[1] &
        pt[1] <= box[3] &
        pt[2] >= box[2] &
        pt[2] <= box[4])	{ is_in <- TRUE }
  }
  return(is_in)
}



#' junction
#' 
#' Finds adjacent stream lines.
#' 
#' @param id is the id number of the stream
#' @param pt is the current position
#' @param strm_cords is a matrix of stream geometry
#' @param thresh is the connection threshhold
#' @return matrix of adjacent stream coordinates
junction <- function (id, pt, strm_cords, thresh) {
  # print("calling junction")
  nbrs <- strm_cords[strm_cords[,3] != id, ]
  strm_cols <- length(strm_cords[1, ])
  strm_dif <- raster::pointDistance(pt[1:2], nbrs[,1:2], lonlat = F)
  # if stream(s) adjoin
  if (min(strm_dif) < thresh) {
    min_nbr <- nbrs[strm_dif == min(strm_dif), ]
    # if more than one adjoining stream
    if (length(min_nbr) > strm_cols) {
      nbr_lns <- apply(min_nbr, 1, function(x) nrow(strm_cords[strm_cords[,3] == x[3], ]))
      min_nbr <- min_nbr[nbr_lns == max(nbr_lns), ]  # select longest stream
    }
    nbr <- strm_cords[strm_cords[,3] == min_nbr[3], ]
  } else { # if no streams adjoin break
    stop('Encountered stream edge.')
  }
  nbr
}


#' logical index
#' 
#' Given a vector of logical booleans, returns the indices of true elements
#' 
#' @param bool is a vector of logical operators
#' @return the indices of true elements along the vector `log`
#' @importFrom magrittr %>%
logindex <- function(bool,ids=0)	{
  if (length(dim(bool)) < 2)	{
    for (i in 1:length(bool))	{
      if(bool[i]) ids <- c(ids,i)
    }
  }
  if (length(dim(bool)) > 1)	{
    for ( i in 1:nrow(bool))	{
      if(bool[i,1] & bool[i,2]) ids <- c(ids,i)
    }
  }
  
  ids <- ids[-1]
  return(ids)
}




#' match_crs
#'
#' Make sure your extents overlap by converting polygons
#' into a common reference crs.
#'
#' @param so is a spatial polygons objects
#' @param target is a spatial object with the target crs projection
#' @return spatial object projected in the target crs
#' @seealso pred_cover1
match_crs <- function(so, target)  {
  crs_ref <- sf::st_crs(target)
  if (class(so)[1] %in% c('sf', 'sfc', 'sfg'))  {
    so <- sf::st_transform(so, crs_ref)
  }
  if (class(so)[1] %in% c('RasterLayer', 'RasterStack'))  {
    so <- raster::projectRaster(so, as(target, 'Spatial'))
  }
  if (class(so)[1] %in% c('SpatialPolygons', 'SpatialPolygonsDataFrame'))  {
    so <- sp::spTransform(so, crs_ref[[1]])
  }
  so
}


#' orient
#' 
#' Orient a matrix for forward traverse.
#' 
#' @param pt is the current position
#' @param mat is the line matrix to traverse
#' @param bear is the bearing of travel
#' @return matrix leading from pt in direction of travel
orient <- function(pt, mat, bear) {
  # print('calling orient')
  difs <- apply(mat, 1, function(x) raster::pointDistance(pt[1:2], x[1:2], lonlat = F))
  if (sum(pt %in% mat) == length(pt)) {
    # print('point is in mat')
    # print('dif length')
    # print(length(difs))
    idx <- sum((1:length(difs)) * (difs == min(difs)))  # id of closest pt on mat
  } else {
    # print('point is not in mat')
    mat_pt <- insert_pt(pt, mat)
    mat <- mat_pt[[1]]
    idx <- mat_pt[[2]]
  }
  if (idx > 1 & idx < length(difs)) { # if pt is in middle of mat
    # print('pt is in inner row')
    br <- bearing(mat[idx-1, 1:2], mat[idx, 1:2])  # taking bearing heading down rows
    junct_angl <- comp_angle(br, bear)
    if (junct_angl < 90) {  # if junction angle is acute, continue heading down
      # print('continue heading down')
      mat <- mat[idx:nrow(mat), ]
      # print('mat length')
      # print(nrow(mat))
    } else {  # if obtuse, head up
      # print('head up')
      mat <- mat[idx:1, ]
    }
  }
  if (idx == length(difs)) {
    # print('idx is end of mat')
    mat <- mat[idx:1, ]
  }
  mat
}


#' pad
#' 
#' pad sf object and return extent
#' 
#' @param so is an sf object
#' @param pad is the factor of padding
#' @return an sf boundary box
#' @seealso fill_extent
pad <- function(so, pad = .1) {
  box <- sf::st_bbox(so)
  xdif <- (box$xmax - box$xmin) * pad
  buff <- sf::st_buffer(so, xdif)
  sf::st_bbox(buff)
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
  pal <- get_palette(c('crimson', 'gold', 'forest'), .7)
  files <- get_ras(in_path)
  crs_ref <- sf::st_crs(raster::raster(file.path(in_path, files[1])))[[1]]
  polys <- sf::st_transform(polys, crs_ref)
  if (method == 'ml') model <- keras::load_model_hdf5(mod_path)
  obs <- array(0, c(length(polys), 50))
  ob <- 0
  
  for (i in seq_along(polys))	{
    sam <- fill_extent(polys[i, ], in_path)
    if (!is.null(sam)) {
      labcords <- matrix(0, nrow = 50, ncol = 2)
      print('opening png device')
      png(file.path(out_path, paste0('sample_', i, '_', year, '.png')),
          width = 36, height = 36, units = 'cm', res = 300)
      raster::plotRGB(sam, r = 1, g = 2, b = 3, main = paste0('Sample Site ',i))
      
      for (j in 1:50)  {
        box <- polys[[i]][[j]]
        box <- sf::st_sfc(sf::st_polygon(box), crs = crs_ref)
        area <- sf::st_coordinates(box)
        leg1 <- raster::pointDistance(area[1,1:2], area[2,1:2], lonlat = F)
        leg2 <- raster::pointDistance(area[2,1:2], area[3,1:2], lonlat = F)
        frm <- sf::st_bbox(box)
        slc <- raster::mask(sam, as(box, 'Spatial'))
        slc <- raster::crop(slc, raster::extent(as(box, 'Spatial')))
        
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
        
        
        lines(area, col = get_palette('slate', .5), lwd = 3)
        if (leg1 > leg2) {
          lines(area[2:3, ], col = pal[cov+1], lwd = 7)
          lines(area[4:5, ], col = pal[cov+1], lwd = 7)
        } else {
          lines(area[1:2, ], col = pal[cov+1], lwd = 7)
          lines(area[3:4, ], col = pal[cov+1], lwd = 7)
        }
        labcords[j,] <- sf::st_coordinates(sf::st_centroid(box))
      }
      
      text(labcords, as.character(1:50), cex = 1.75, col = 'white')
      dev.off()
      print('png made')
    } else {
      message(paste0('Raster not found for ', i, '.  Trying next sample.'))
    }
  }
  
  df <- as.data.frame(obs)
  colnames(df) <- 1:50
  df$id <- 1:nrow(df)
  df$year <- year
  df$type <- type
  df <- df[ , c(51:53,1:50)]
  print('writing csv)')
  write.csv(df, file = file.path(out_path, paste0('samples_', year, '.csv')))
  df
}

#' plot_scores
#'
#' Superimposes sampling boxes over orthographic data.
#' Prints cover scores over each box.  
#' `scores` matches output format of plot_samples().  
#' `scores` must have same number of rows as `polys`, in the same order.
#'
#' @param in_path is a character vector containing the path to the orthoimagery
#' @param out_path is a character vector specifying the output directory for the plots
#' @param polys is a SpatialPolygonsDataFrame object (the sample boxes shapefile)
#' @param scores is a data.table recording cover score observations
#' @return prints rgb plots of `polys` to `out_path` with scores indicated by color
#' @importFrom magrittr %>%
#' @export
plot_scores <- function(in_path, 
                        out_path, 
                        scores,
                        polys = samples,
                        title = 'samples_')  {
  files <- get_ras(in_path)
  crs_ref <- sf::st_crs(raster::raster(file.path(in_path, files[1])))
  polys <- sf::st_transform(polys, crs_ref)
  scores <- as.matrix(scores[ , 4:53])
  
  for (i in seq_along(polys))	{
    sam <- fill_extent(polys[i, ], in_path)
    if (!is.null(sam)) {
      
      labcords <- matrix(0, nrow = 50, ncol = 2)
      pal <- get_palette(c('crimson', 'gold', 'forest'), .7)
      png(file.path(out_path, paste0(title ,i,'.png')),
          width = 36, height = 36, units = 'cm', res = 300)
      raster::plotRGB(sam, r = 1, g = 2, b = 3, main = paste0('Sample Site ',i))
      
      for (j in 1:50)  {
        box <- polys[[i]][[j]]
        box <- sf::st_sfc(sf::st_polygon(box), crs = crs_ref)
        area <- sf::st_coordinates(box)
        leg1 <- raster::pointDistance(area[1,1:2], area[2,1:2], lonlat = F)
        leg2 <- raster::pointDistance(area[2,1:2], area[3,1:2], lonlat = F)
        
        lines(area, col = get_palette('slate', .5), lwd = 3)
        if (leg1 > leg2) {
          lines(area[2:3, ], col = pal[scores[i, j]+1], lwd = 7)
          lines(area[4:5, ], col = pal[scores[i, j]+1], lwd = 7)
        } else {
          lines(area[1:2, ], col = pal[as.numeric(scores[i, j])+1], lwd = 7)
          lines(area[3:4, ], col = pal[as.numeric(scores[i, j])+1], lwd = 7)
        }
        labcords[j,] <- sf::st_coordinates(sf::st_centroid(box))
      }
      
      text(labcords, as.character(1:50), cex = 1.75, col = 'white')
      dev.off()
    }
  }
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
  year1_files <- get_ras(year1_path)
  year2_files <- get_ras(year2_path)
  for (i in seq_along(year1_files))  {
    year1 <- raster::raster(file.path(year1_path, year1_files[i]))
    setwd(year2_path)
    year2 <- raster::raster(file.path(year2_path, year2_files[i]))
    year1 <- raster::projectRaster(year1, year2)
    chng <- year2 - year1
    raster::writeRaster(chng, file.path(out_path, year1_files[i]), 'GTiff')
  }
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
                               lots = sf_lots, 
                               title = 'pred_chng_table.csv')  {
  og_dir <- getwd()
  files <- get_ras(chng_path)
  crs_ref <- sf::st_crs(raster::raster(file.path(chng_path, files[1])))[[1]]
  lots <- sf::st_transform(lots, crs_ref)
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
    ras <- raster::raster(file.path(chng_path, files[i]))
    vals <- raster::values(ras)
    chng[i] <- sum(vals[!is.na(vals)]) / 2
    area[i] <- length(vals[!is.na(vals)])
    cover[i] <- (sum(vals[!is.na(vals)])/2) / area[i]
    lot <- suppressWarnings(sf::st_crop(lots, ras))
    map_nos <- levels(factor(lot$MapTaxlot))
    mtls[i] <- squash(map_nos)
    perms[i] <- squash(permits$record_no[permits$maptaxlot %in% map_nos])
    
    ras <- raster::raster(file.path(year1_path, files[i]))
    vals <- raster::values(ras)
    vals[vals<0] <- 0
    vals[vals>2] <- 2
    area_yr1[i] <- length(vals[!is.na(vals)])
    cov_yr1[i] <- sum(vals[!is.na(vals)]) / 2
    
    ras <- raster::raster(file.path(year2_path, files[i]))
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
  write.csv(dt, title)
  print(paste0('Predicted cover extent year 1:  ', round(tot_yr1, 2)))
  print(paste0('Predicted cover extent year 2:  ', round(tot_yr2, 2)))
  print(paste0('Predicted cover change:  ', round(tot_chng, 2)))
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
pred_cover1 <- function(poly,
                        in_dir = NULL, 
                        out_dir = NULL, 
                        rgb_dir = NULL, 
                        cir_dir = NULL, 
                        buffer = sf_buff)  {
  
  if(!is.null(rgb_dir)) {
    r <- fill_extent(poly, rgb_dir)
    red <- raster::raster(r, layer = 1)
    grn <- raster::raster(r, layer = 2)
    blu <- raster::raster(r, layer = 3)
  }
  if(!is.null(cir_dir))  {
    r <- fill_extent(poly, cir_dir)
    nir <- raster::raster(r, layer = 1)
  } else {
    r <- fill_extent(poly, in_dir)
    nir <- raster::raster(r, layer = 4)
    red <- raster::raster(r, layer = 1)
    grn <- raster::raster(r, layer = 2)
    blu <- raster::raster(r, layer = 3)
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
  ras <- raster::mask(ras, as(poly, 'Spatial'))
  buffer <- match_crs(buffer, ras)
  ras <- raster::mask(ras, as(buffer, 'Spatial'))
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
                       buff = sf_buff)  {
  if (!is.null(in_path)) crs_ref <- get_crs(in_path)
  if (!is.null(rgb_path)) crs_ref <- get_crs(rgb_path)
  for (i in seq_along(polys))  {
    ras <- pred_cover1(poly = polys[i, ], 
                       in_dir = in_path, 
                       out_dir = out_path,
                       rgb_dir = rgb_path,
                       cir_dir = cir_path,
                       buffer = buff) 
    
    raster::writeRaster(ras, file.path(out_path, paste0('site_',i,'.tif')), 'GTiff', overwrite = T)
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
pred_cover_report <- function(in_path, taxlots = sf_lots)  {
  files <- get_ras(in_path)
  crs_ref <- sf::st_crs(raster::raster(file.path(in_path, files[1])))[[1]]
  taxlots <- sf::st_transform(taxlots, crs_ref)
  cover <- array(0, length(files))
  area <- array(0, length(files))
  mtls <- vector(length(files), mode = 'character')
  for (i in seq_along(files))  {
    ras <- raster::raster(file.path(in_path, files[i]))
    vals <- raster::values(ras)
    vals[vals<0] <- 0
    vals[vals>2] <- 2
    cover[i] <- sum(vals[!is.na(vals)])
    area[i] <- length(vals[!is.na(vals)])
    lot <- suppressWarnings(sf::st_crop(taxlots, ras))
    mtls[i] <- squash(levels(factor(lot$MapTaxlot)))
  }
  pct_cov <- mapply(function(a, b) (sum(a)/2) / b, a = cover, b = area)
  dt <- data.table::setDT(data.frame(lot = mtls, cover = pct_cov, area = area))
  write.csv(dt, 'pred_cov.csv')
  print(paste0('Predicted Cover Extent:  ', round(mean(!is.na(pct_cov))*100, 2), '%'))
  
}



#' sample_streams
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
                           buff = per_buff,
                           width = 50,
                           span = 75,
                           inc = 6,
                           steps = 25, 
                           thresh = 50)  {
  # convert sp files to sf
  if (!('sfc' %in% class(strms))) strms <- sf::st_as_sf(strms)
  if (!('sfc' %in% class(lots))) lots <- sf::st_as_sf(lots)
  if (!('sfc' %in% class(prc))) prc <- sf::st_as_sf(prc)
  if (!('sfc' %in% class(buff))) buff <- sf::st_as_sf(buff)
  crs_ref <- sf::st_crs(strms)
  
  # transform all data to common crs
  lots <- sf::st_transform(lots, crs_ref)
  prc <- sf::st_transform(prc, crs_ref)
  buff <- sf::st_transform(buff, crs_ref)
  
  # select only parts of the riparian buffer within lots of interest
  poly <- sf::st_intersection(sf::st_geometry(buff), sf::st_geometry(lots))
  prc <- prc[poly,1]
  
  # sample points from selected area
  prc_geom <- sf::st_geometry(prc)
  
  ct <- 1
  samples <- list()
  
  while (ct <= n) {
    print('making sample')
    print(ct)
    pt <- lapply(
      lapply(sample(1:length(prc_geom), 1),
             function(x) sf::st_sample(prc_geom[[x]], 1)
      ),
      function(y) sf::st_coordinates(y)
    )
    pt <- pt[[1]]
    
    # returns list of {list: sample points, list: stream coord matrix}
    pt <- find_line(pt, strms)
    center <- tryCatch(centerline(pt, span = 75, inc = 6, strms, thresh),
                       error = function(cond) {
                         message(cond)
                         return(NA)
                       })
    if (is.na(center[1])) next
    c_list <- tryCatch(corners(center, width),
                       error = function(cond) {
                         message(cond)
                         return(NA)
                       })
    if (is.na(c_list[1])) next
    c_buff <- c_list[[1]]
    c_corners <- c_list[[2]]
    
    lb <- tryCatch(short_side(c_corners[1, ], c_corners[2, ], c_buff),
                   error = function(cond) {
                     message(cond)
                     return(NA)
                   })
    if (is.na(lb[1])) next
    rev <- nrow(lb):1
    lb <- lb[order(rev), ]
    rb <- tryCatch(short_side(c_corners[3, ], c_corners[4, ], c_buff),
                   error = function(cond) {
                     message(cond)
                     return(NA)
                   })
    if (is.na(rb[1])) next
    lb_seg <- tryCatch(segment(lb, steps),
                       error = function(cond) {
                         message(cond)
                         return(NA)
                       })
    if (is.na(lb_seg[1])) next
    rb_seg <- tryCatch(segment(rb, steps),
                       error = function(cond) {
                         message(cond)
                         return(NA)
                       })
    if (is.na(rb_seg[1])) next
    boxes <- wrap(lb_seg, center)
    r_boxes <- wrap(center, rb_seg)
    bx_ln <- length(boxes)
    for (i in seq_along(r_boxes)) {
      boxes[[bx_ln + i]] <- r_boxes[[i]]
    }
    samples[[ct]] <- sf::st_multipolygon(boxes)
    ct <- ct + 1
  }
  samples <- sf::st_sfc(samples, crs = crs_ref)
}  



#' segment
#' 
#' Divide a line into even segments.
#' 
#' @param mat is the matrix of line coordinates
#' @param steps is the segment increment
#' @return matrix of point coordinates for segments
segment <- function(mat, steps) {
  print('segmenting')
  seg <- mat[1, ]
  ln <- sf::st_length(sf::st_linestring(mat)) / steps
  pt <- walk(mat, ln)
  seg <- rbind(seg, pt[[1]])
  if (steps > 1) {
    for (i in 2:steps) {
      mat <- orient(pt[[1]], mat, pt[[3]])
      pt <- walk(mat, ln)
      seg <- rbind(seg, pt[[1]])
    }
  }
  seg
}



#' short_side
#' 
#' Given two points on a polygon, return the shorter path between them.
#' 
#' @param x is a point with elements c(x, y)
#' @param y is a point with elements c(x, y)
#' @param circle is a matrix of polygon coordinates
#' @return the subset of `circle` containing the shorter path between `x` and `y`
short_side <- function(x, y, circle) {
  print('finding short side')
  difs <- raster::pointDistance(x[1:2], circle[, 1:2], lonlat = F)
  idx <- sum((1:length(difs)) * (difs == min(difs)))
  difs <- raster::pointDistance(y[1:2], circle[, 1:2], lonlat = F)
  idy <- sum((1:length(difs)) * (difs == min(difs)))
  mn <- min(idx, idy)
  mx <- max(idx, idy)
  
  m1 <- circle[mn:mx, ]
  m2 <- rbind(circle[mx:nrow(circle), ], circle[1:mn, ])
  ln1 <- sf::st_length(sf::st_linestring(m1))
  ln2 <- sf::st_length(sf::st_linestring(m2))
  
  if (ln1 < ln2) {
    return(m1)
  } else {
    return(m2)
  }
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



#' touches
#' 
#' one touches other
#' 
#' @param box1 is sf bbox
#' @param box2 is sf bbox
#' @return boolean touches
#' @export
touches <- function(box1, box2) {
  touches <- FALSE
  if (is_in(box1[c(1,2)], box2)) touches <- TRUE
  if (is_in(box1[c(1,4)], box2)) touches <- TRUE
  if (is_in(box1[c(3,2)], box2)) touches <- TRUE
  if (is_in(box1[c(3,4)], box2)) touches <- TRUE
  touches
}


#' thumbnails
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
#' @export
thumbnails <- function(in_path, out_path, polys = random_samples, output = 'tif')  {
  for (i in seq_along(polys)) {
    r <- fill_extent(polys[i, ], in_path, pad = .2)
    if (!is.null(r)) {
      print(paste0('filling ', i))
      if (output == 'tif') {
        raster::writeRaster(
          r, filename = file.path(out_path, paste0('sample_', i)),
          format = 'GTiff')
      }
      if (output == 'png') {
        polys <- sf::st_transform(polys, sf::st_crs(r))
        buff <- sf::st_transform(sf_buff, sf::st_crs(r))
        png(file.path(out_path, paste0('sample_', i, '.png')),
            height = 36, width = 36, units = 'cm', res = 300)
        raster::plotRGB(r)
        plot(polys[i,1], lwd = 2, col = get_palette('slate', .01), add = T)
        plot(buff, col = get_palette('ocean', .15), add = T)
        dev.off()
      }
    }
  }
}


#' traverse
#' 
#' Traverse down a stream network.
#' 
#' @param mat is the current oriented stream
#' @param dist is the distance to travel
#' @param mats is the stream network
#' @param thresh is the connection threshold
#' @return a list of c(next pt, distance remaining, orientation of travel)
#' @seealso orient
#' @seealso walk
traverse <- function(mat, dist, mats, thresh) {
  # print('calling traverse')
  next_pt <- walk(mat, dist)
  if (next_pt[[2]] > 0) {
    bear <- bearing(mat[(nrow(mat)-1), 1:2], mat[nrow(mat), 1:2])
    nbr <- junction(next_pt[[1]][3], next_pt[[1]], mats, thresh)
    nbr <- orient(next_pt[[1]], nbr, bear)
    nbr <- sf::st_segmentize(sf::st_sfc(sf::st_linestring(nbr)), 1)[[1]]
    next_pt <- traverse(nbr, next_pt[[2]], mats, thresh)
  }
  next_pt
}


#' triangle length
#' 
#' Returns the point a given distance and bearing from a coordinate pair.
#' 
#' @param pt is a coordinate pair
#' @param bear is a bearing in degrees
#' @param dist is a distance in meters
#' @return coordinates of the point `dist` from `pt` along bearing `bear`
#' @export
tri_length <- function(pt,bear,dist,north=TRUE)	{
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



#' uturn
#' 
#' given a bearing, return a 180 degree change in bearing
#' bearings in degrees
#' 
#' @param bear is a bearing in degrees
#' @return a 180 degree change in bearing from `bear`
uturn <- function(bear) {
  bear <- bear + 180
  if (bear > 360) bear <- bear - 360
  bear
}



#' walk
#' 
#' Traverse from top row to bottom row of matrix.  
#' If distance `dist` lies on matrix `mat`,  
#' return interpolated pt and remainder distance of zero.
#' If `dist` is farther than the end of `mat`, 
#' return end of mat as pt and remainder distance.
#' 
#' @param mat is a stream path
#' @param dist is the length to travel down the stream path
#' @return a list c(new position, distance remaining, orientation of travel)
walk <- function(mat, dist) {
  # print('calling walk')
  # interpolated distance between points in the stream line
  lens <- 0
  for (i in 2:nrow(mat)) {
    lens[i] <- sf::st_length(sf::st_sfc(sf::st_linestring(mat[(i-1):i, 1:2])))
  }
  # cumulative distance along the stream line
  difs <- cumsum(lens)
  if (max(difs) > dist) { # point lies on mat
    # print('point lies on mat')
    target <- difs - dist  # distance from target
    abv <- min(target[target > 0]) # closest pt past target
    abv <- sum((1:length(target) * (target == abv))) # convert to row id
    blw <- max(target[target < 0]) # closest pt before target
    blw <- sum((1:length(target) * (target == blw))) # convert to row id
    bear <- bearing(mat[blw, 1:2], mat[abv, 1:2])
    rem <- dist - difs[blw]
    pt <- tri_length(mat[(abv-1), ], bear, rem)
    pt <- c(pt, mat[(blw), -c(1:2)])
    rem <- 0
    
  } else { # return end of matrix at remainder distance
    # print('pt not on mat')
    pt <- mat[nrow(mat), ]
    bear <- bearing(mat[(nrow(mat)-1), 1:2], mat[nrow(mat), 1:2])
    rem <- difs[nrow(mat)]
  }
  list(pt, rem, bear)
}


#' within
#' 
#' small within big
#' 
#' @param small is sf bbox
#' @param big is sf bbox
#' @return boolean is in
#' @export
within <- function(small, big) {
  within <- FALSE
  if (is_in(box1[c(1,2)], box2) &
      is_in(box1[c(1,4)], box2) &
      is_in(box1[c(3,2)], box2) &
      is_in(box1[c(3,4)], box2)) within <- TRUE
  within
}


#' wrap
#' 
#' Wrap two matrices of equal row number into box polygons.
#' 
#' @param x is a matrix with cols c(x,y)
#' @param y is a matrix with cols c(x,y)
#' @return a list of polgyons
wrap <- function(x, y) {
  boxes <- list()
  for (i in 1:(nrow(x)-1)) {
    box <- matrix(0, nrow = 5, ncol = 2)
    box[1, ] <- x[i, 1:2]
    box[2, ] <- y[i, 1:2]
    box[3, ] <- y[i+1, 1:2]
    box[4, ] <- x[i+1, 1:2]
    box[5, ] <- x[i, 1:2]
    boxes[[i]] <- sf::st_polygon(list(box))
  }
  boxes
}
















































