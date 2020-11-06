

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
#' @importFrom magrittr %>%
#' @export


thumbnails <- function(in_path, out_path, polys = samples)  {
  files <- get_ras(in_path)
  crs_ref <- raster::crs(raster::raster(file.path(in_path, files[1])))
  polys <- sp::spTransform(polys, crs_ref)

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
      ras <- raster::stack(file.path(in_path, files[hit[1]]))
      ras <- raster::crop(ras, frame)
    }

    if (length(hit) > 1)	{
      for (k in 2:length(hit))	{
        r <- raster::stack(file.path(in_path, files[hit[k]]))
        r <- raster::crop(r, frame)
        ras <- raster::mosaic(ras, r, fun = max)
      }
    }

    ras <- raster::crop(ras, frame)
    raster::writeRaster(ras, filename = file.path(out_path, paste0('sample_', i)),
                        format = 'GTiff')
  }
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
  crs_ref <- raster::crs(raster::raster(file.path(in_path, files[1])))
  polys <- sp::spTransform(polys, crs_ref)
  scores <- as.matrix(scores[,4:53])
  
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
    
    
    png(file.path(out_path, paste0(title ,i,'.png')))
    raster::plotRGB(sam, r = 1, g = 2, b = 3, main = paste0('Sample Site ',i))
    area <- lapply(methods::slot(polys[i,], 'polygons'),
                   function(x) lapply(methods::slot(x, 'Polygons'),
                                      function(y) methods::slot(y, 'coords')))
    area <- area[[1]]
    
    for (j in 1:50)  {
      pal <- get_palette(c('crimson', 'gold', 'forest'), .7)
      lines(area[[j]], col = get_palette('slate', .5), lwd = 3)
      lines(area[[j]][1:2, ], col = pal[scores[i,j]+1], lwd = 7)
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
}



#' color_array
#'
#' extracts color statistics from raster image
#'
#' @param img is a raster stack of 4 bands rgbn
#' @return an array of summary color statistics
#' @export
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
#' @export
color_array_3band <- function(img) {
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
  df <- as.data.frame(car)
  colnames(df) <- c('red_mn', 'red_sd', 'red_sm',
                    'grn_mn', 'grn_sd', 'grn_sm',
                    'blu_mn', 'blu_sd', 'blu_sm')
  df
}



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








