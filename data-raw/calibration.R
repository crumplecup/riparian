setwd('E:/Riparian')
write("TMP = 'E:/Riparian/temp'", file=file.path(Sys.getenv('R_USER'), '.Renviron'))
library(riparian)


samples <- readOGR('samples_nad83.shp')


bit <- get_cols_4band(samples[1,], 'E:/ortho2018')

















