#' annual random monitoring samples
#' 
#' a spatial object containing the 54 sampling boxes
#' for monitoring the priority riparian corridor
#' 
#' @format an sf object with 54 multipolygons of 50 polygons each
'random_samples'

#' riparian buffer
#' 
#' the 50-ft riparian buffer cropped to the priority riparian corridor
#' 
#' @format an sf multipolygon
'sf_buff'

#' priority corridor taxlots
#' 
#' maptaxlots within the priority riparian corridor
#' 
#' @format an sf geometry collection of multipolygons with 630 obs. of 37 vars
'sf_lots'

#' perennial streams
#' 
#' perennial stream layer from NHD
#' 
#' @format an sf geometry collection of multilinestrings with 127 obs of 38 vars
'sf_perennial'

#' hydro-enforced drainage
#' 
#' hydro-enforced drainage layer from WSI
#' 
#' @format an sf multilinestring
'sf_streams'