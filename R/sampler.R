






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



#' set flag
#' 
#' 
#' 
#' @param pt is a coordinate pair (starting point)
#' @param mat is a matrix of coords representing a line with cols (x,y)
#' @param dist is a distance in meters
#' @param us is a logical boolean indicating whether distance is upstream or downstream
#' @return the point distance `dist` from starting point `pt` along line `mat`
#' @importFrom magrittr %>%

set_flag <- function(pt,mat,dist=75*0.3048,us=F,k=0)	{
  mat <- streamsnap(pt, mat)
  difs <- raster::pointDistance(pt, mat, lonlat=F)
  for (i in 1:nrow(mat))	{
    if (difs[i] == min(difs)) ptid <- i
  }
  
  if (us) mat <- mat[ptid:nrow(mat), ]
  if (!us) mat <- mat[ptid:1, ]
  
  so_far <- 0
  while (so_far <= dist)	{
    k <- k + 1
    rem <- dist - so_far
    so_far <- so_far + norm(matrix(mat[k+1,] - mat[k,], ncol = 2), type = '2')
  }
  if (nrow(mat) >= (k+1))	{
    bear <- get_bear(mat[k:(k+1),])
  }
  if (nrow(mat) == k)	{
    bear <- get_bear(mat[(k-1):k,])
  }
  
  pt <- tri_length(mat[k,], bear, rem)
  return(pt)
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


#' snap point to stream
#'
#' Given a coordinate pair `pt` and a matrix of coords `mat`, 
#' returns the point on `mat` nearest `pt`.
#'
#' @param pt is a coordinate pair (numeric vector)
#' @param mat is a matrix of coords representing the stream path
#' @return the point on the stream path closest to `pt`
#' @export

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



#' find line
#' 
#' finds the nearest line in a sp object, returns coordinate matrix
#' 
#' @param pt is a coordinate pair (x,y)
#' @param lines is a SpatialLinesDataFrame
#' @return coordinate matrix of nearest line to `pt` in `lines`
#' @importFrom raster pointDistance
#' @importFrom magrittr %>%

find_line <- function (pt, lines = prc_strms)  {
  lines <- coordinates(lines)
  lines <- lines[[1]]
  difs <- lapply(lines, function(a) min(apply(a, 1, function(b) pointDistance(b, pt, lonlat = F))))
  difs <- unlist(difs)
  id <- logindex(difs == min(difs))
  lines[[id]]
}

pt_to_line <- function(pt, line)  {
  difs <- apply(line, 1, function(a) pointDistance(a, pt, lonlat = F))
  nebr <- logindex(difs == min(difs))
  line[nebr,]
}

sample_streams <- function(n = 100, per = prc_per, strms = prc_strms)  {
  crs_ref <- raster::crs(strms)
  pts <- sp::spsample(per, 1000, 'random')
  pts <- sample(pts, n)
  pt_cords <- coordinates(pts)
  finds <- apply(pt_cords, 1, function(a) find_line(a, strms))
  pts <- matrix(0, nrow = length(finds), ncol = 2)
  for (i in seq_along(finds)) {
    pts[i,] <- pt_to_line(pt_cords[i,], finds[[i]])
  }
  boxes <- centerline(pts, finds)
}










centerline <- function(pts, strms) {
  boxes <- list()
  for(i in 1:nrow(pts))	{
    box <- matrix(0, ncol=2)
    try(	box <- box %>% rbind(pts[i,] %>% set_flag(strms[[i]], us=T)))
    for (j in 1:25)	{
      try(box <- box %>% rbind(box[nrow(box),] %>% set_flag(strms[[i]],dist= 6*0.3048)))
    }
    box <- box[-1,]
    boxes[[i]] <- box
  }
  
  short <- list()
  k <- 1
  for (i in 1:length(boxes))	{
    if(length(boxes[[i]])>1)	{
      short[[k]] <- boxes[[i]]
      k <- k+1
    }
  }
  short
}






sides <- function(mat) {
  lefts <- list()
  rights <- list()
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
  buff <- mat %>% Line %>% list %>% Lines(ID='first') %>% list %>% SpatialLines(crs_ref)
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
  
  #find length of 'left' side
  buff <- rbind(buff, buff)
  above <- TRUE
  done <- FALSE
  k <- luid
  while(above)	{
    k <- k + 1
    if (left_down[1]==buff[k,1] & left_down[2]==buff[k,2])	{
      seg <- buff[luid:k,]
      above <- FALSE
      done <- TRUE
    }
    if (right_up[1]==buff[k,1] & right_up[2]==buff[k,2]) above <- FALSE
  }
  k <- luid + nrow(buff)*0.5
  while(!above & !done)	{
    k <- k - 1
    if (left_down[1]==buff[k,1] & left_down[2]==buff[k,2])	{
      seg <- buff[(luid+nrow(buff)*0.5):k,]
      above <- TRUE
    }
  }
  seglen <- seg %>% Line %>% list %>% Lines(ID='first') %>% list %>% SpatialLines(crs_ref) %>% lineLength / 25
  
  #record 'left' side points
  left_side <- left_up %>% as.vector %>% as.matrix %>% t
  for (j in 1:24)	{
    left_side <- left_side %>% rbind(left_side[nrow(left_side), ] %>% set_pt(mat=seg, dist=seglen))
  }
  left_side <- left_side %>% rbind(left_down)
  lefts[[i]] <- left_side
  
  #find length of 'right' side
  above <- TRUE
  done <- FALSE
  k <- ruid
  while(above)	{
    k <- k + 1
    if (right_down[1] == buff[k, 1] & right_down[2] == buff[k, 2])	{
      seg <- buff[ruid:k, ]
      above <- FALSE
      done <- TRUE
    }
    if (left_up[1] == buff[k,1] & left_up[2] == buff[k,2]) above <- FALSE
  }
  k <- ruid + nrow(buff)*0.5
  while(!above & !done)	{
    k <- k - 1
    if (right_down[1] == buff[k,1] & right_down[2] == buff[k,2])	{
      seg <- buff[(ruid + nrow(buff) * 0.5):k, ]
      above <- TRUE
    }
  }
  seglen <- seg %>% Line %>% list %>% Lines(ID='first') %>% list %>% SpatialLines(crs_ref) %>% lineLength / 25
  
  #record 'right' side points
  right_side <- right_up %>% as.vector %>% as.matrix %>% t
  for (j in 1:24)	{
    right_side <- right_side %>% rbind(right_side[nrow(right_side), ] %>% set_pt(mat=seg, dist=seglen))
  }
  right_side <- right_side %>% rbind(right_down)
  rights[[i]] <- right_side
}

  
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

















