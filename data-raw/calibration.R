setwd('E:/Riparian')
write("TMP = 'E:/Riparian/temp'", file=file.path(Sys.getenv('R_USER'), '.Renviron'))
library(riparian)
library(rgdal)
library(data.table)

samples <- readOGR('samples_nad83.shp')
res <- fread('samples2018.csv')                  #results of sampling observations


# subset results by year

res18 <- res[year == 2018]
res16 <- res[year == 2016]
res09 <- res[year == 2009]

# divide into sampling groups

res18_a <- res18[1:81]
res18_b <- res18[82:116]
res18_c <- res18[117:170]
res16_a <- res16[1:79]
res16_b <- res16[80:115]
res16_c <- res16[116:169]

res18_a$poly_id <- 1:81
res18_b$poly_id <- c(1:25,27,29:37)
res18_c$poly_id <- c(1:5,7:29,31:56)
res16_a$poly_id <- c(1:18,20:44,46:81)
res16_b$poly_id <- c(1:25,27:37)
res16_c$poly_id <- c(1:5,7:29,31:56)

#exclude final group from earlier groups

res18_a <- res18_a[!res18_a$poly_id %in% c(1,5:7,10:12,15,19,23,28:29,31:33,41:42,61,63,66,71:72,75:78)]
res18_b <- res18_b[!res18_b$poly_id %in% c(8,10,15,22,24:26,33:34)]
res18_c <- res18_c[!res18_c$poly_id %in% c(3,8,11,14,16,18:19,21,25,29:30,32,37,39,44:45,47:48,50:51,53)]
res16_a <- res16_a[!res16_a$poly_id %in% c(1,5:7,10:12,15,19,23,28:29,31:33,41:42,61,63,66,71:72,75:78)]
res16_b <- res16_b[!res16_b$poly_id %in% c(8,10,15,22,24:26,33:34)]
res16_c <- res16_c[!res16_c$poly_id %in% c(3,8,11,14,16,18:19,21,25,29:30,32,37,39,44:45,47:48,50:51,53)]



for (i in seq_along(samples))  {
  if (i == 1)  {
    dt <- get_cols_4band(samples[i,], 'E:/ortho2018')
  }
  if (i >1)  {
    dt <- rbind(dt, get_cols_4band(samples[i,], 'E:/ortho2018'))
  }
}

usethis::use_data(dt)


samples_09 <- readOGR('samples_ortho2009.shp')

for (i in seq_along(samples_09))  {
  if (i == 1)  {
    dt_09 <- get_cols_4band(samples_09[i,], 'E:/Riparian/ortho2009')
  }
  if (i >1)  {
    dt_09 <- rbind(dt_09, get_cols_4band(samples_09[i,], 'E:/Riparian/ortho2009'))
  }
}

usethis::use_data(dt_09)


groupA <- readOGR('groupA_ortho2018.shp')
groupA <- sp::geometry(groupA)
groupB <- readOGR('groupB_ortho2018.shp')
groupB <- sp::geometry(groupB)
groupC <- readOGR('groupC_ortho2018.shp')
groupC <- sp::geometry(groupC)

selA <- 1:81
selA <- selA[!selA %in% c(1,5:7,10:12,15,19,23,28:29,31:33,41:42,61,63,66,71:72,75:78)]

for (i in seq_along(selA))  {
  if (i == 1)  {
    dtA_18 <- get_cols_4band(groupA[selA[i],], 'E:/ortho2018')
  }
  if (i >1)  {
    dtA_18 <- rbind(dtA_18, get_cols_4band(groupA[selA[i],], 'E:/ortho2018'))
  }
}

usethis::use_data(dtA_18)


for (i in seq_along(selA))  {
  print(i)
  print(selA[i])
}


selA[i]

names(groupA)









