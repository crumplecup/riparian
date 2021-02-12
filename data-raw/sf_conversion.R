# convert sp spatial objects to sf



sp_to_sf <- function(sp) {
  crs_ref <- sf::st_crs(sp)
  sf_list <- list()
  polys <- methods::slot(sp, 'polygons')

  for (i in seq_along(polys)) {
    p <- methods::slot(polys[[i]], 'Polygons')
    sf_polys <- list()
    for (j in seq_along(p)) {
      sf_polys[[j]] <- sf::st_polygon(list(methods::slot(p[[j]], 'coords')))
    }
    sf_list[[i]] <- sf::st_multipolygon(sf_polys)
  }
  sf::st_sfc(sf_list, crs = crs_ref)
}







random_samples <- sp_to_sf(samples)  # annual random monitoring boxes
sf_perennial <- sf::st_as_sf(prc_per)  # perennial streams (NHD)
sf_streams <- sf::st_as_sf(prc_strms)  # hydro-enforced drainage layer (WSI)
sf_lots <- sf::st_as_sf(prc_lots) # maptaxlots in the riparian corridor
sf_buff <- sf::st_as_sf(prc_buff) # riparian buffer for corridor

usethis::use_data(random_samples)
usethis::use_data(sf_perennial)
usethis::use_data(sf_streams)
usethis::use_data(sf_lots)
usethis::use_data(sf_buff)




