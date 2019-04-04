setwd('E:/Riparian')
write("TMP = 'E:/Riparian/temp'", file=file.path(Sys.getenv('R_USER'), '.Renviron'))
library(riparian)
library(rgdal)
library(data.table)
library(magrittr)

samples <- readOGR('samples_nad83.shp')
usethis::use_data(samples)

res <- fread('samples2018.csv')                  #results of sampling observations
samples2018 <- fread('samples2018.csv') 
colnames(samples2018) <- c('id','year','type',1:50)
usethis::use_data(samples2018, overwrite = T)

# subset results by year

res18 <- res[year == 2018]
res16 <- res[year == 2016]
res09 <- res[year == 2009]

# divide into sampling groups

res18_p <- res18[sample_id %in% 1:54]

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

#reshape values into vector for insertion in to data.table

vals18_p <- res18_p[,4:53] %>% unlist %>% 
  matrix(ncol=50) %>% t %>% 
  data.frame %>% setDT %>% unlist

vals18_a <- res18_a[,4:53] %>% unlist %>% 
  matrix(ncol=50) %>% t %>% 
  data.frame %>% setDT %>% unlist

vals18_b <- res18_b[,4:53] %>% unlist %>% 
  matrix(ncol=50) %>% t %>% 
  data.frame %>% setDT %>% unlist

vals18_c <- res18_c[,4:53] %>% unlist %>% 
  matrix(ncol=50) %>% t %>% 
  data.frame %>% setDT %>% unlist

dt18_p <- data.frame(rip = vals18_p) %>% cbind(dt)


mat <- matrix(c(rep(1,10), rep(2,10), rep(3,10)), ncol=3)
dt <- setDT(data.frame(mat))
unlist(dt)

for (i in seq_along(samples))  {
  if (i == 1)  {
    dt <- get_cols_4band(samples[i,], 'E:/ortho2018')
  }
  if (i >1)  {
    dt <- rbind(dt, get_cols_4band(samples[i,], 'E:/ortho2018'))
  }
}

dt18_p <- data.frame(rip = vals18_p) %>% cbind(dt)

usethis::use_data(dt18_p)
usethis::use_data(dt)



for (i in seq_along(polys))  {
  print(i)
  if (i == 1)  {
    dt_09 <- get_cols_4band(polys[i,], wd)
  }
  if (i >1)  {
    dt_09 <- rbind(dt_09, get_cols_4band(polys[i,], wd))
  }
}


vals09_p <- res09[,4:53] %>% unlist %>% 
  matrix(ncol=50) %>% t %>% 
  data.frame %>% setDT %>% unlist


dt09 <- data.frame(rip = vals09_p) %>% cbind(dt_09)

usethis::use_data(dt09)


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

dt18_a <- data.frame(rip = vals18_a) %>% cbind(dtA_18)

usethis::use_data(dt18_a)


# Group B 2018

selB <- 1:37
selB <- selB[!selB %in% c(8,10,15,22,24:26,33:34)]

for (i in seq_along(selB))  {
  if (i == 1)  {
    dtB_18 <- get_cols_4band(groupB[selB[i],], 'E:/ortho2018')
  }
  if (i >1)  {
    dtB_18 <- rbind(dtB_18, get_cols_4band(groupB[selB[i],], 'E:/ortho2018'))
  }
}

usethis::use_data(dtB_18)


#load previously scored 2016 data

ndvis <- read.csv('ndvi_samp.csv')
nirs <- read.csv('nir_samp.csv')
reds <- read.csv('red_samp.csv')
greens <- read.csv('grn_samp.csv')
blues <- read.csv('blu_samp.csv')

ndvis <- ndvis[,-1] %>% as.matrix
nirs <- nirs[,-1] %>% as.matrix
reds <- reds[,-1] %>% as.matrix
greens <- greens[,-1] %>% as.matrix
blues <- blues[,-1] %>% as.matrix

for (i in seq_along(selA))  {
  print(i)
  print(selA[i])
}

sampnos <- res_16[,1] %>% as.numeric
res_16 <- res_16[,-c(1:2)] %>% as.matrix

ndvi <- ndvis[1,]
nir <- nirs[1,]
red <- reds[1,]
green <- greens[1,]
blue <- blues[1,]
rip <- res_16[1,]

for (i in 2:nrow(ndvis))	{
  ndvi <- c(ndvi,ndvis[i,])
  nir <- c(nir,nirs[i,])
  red <- c(red,reds[i,])
  green <- c(green,greens[i,])
  blue <- c(blue,blues[i,])
  rip <- c(rip,res_16[i,])
}

selA[i]

names(groupA)

dt16 <- data.frame(rip=unlist(rip),ndvi=ndvi,red=red,grn=green,blu=blue) %>% setDT

usethis::use_data(dt16, overwrite = T)

dt_all <- dt16 %>% rbind(dt18_p %>% rbind(dt18_a))

usethis::use_data(dt_all, overwrite = T)

mod <- lm('rip ~ ndvi + red + grn + blu', data=dt_all)
usethis::use_data(mod, overwrite = T)



summary(mod)
plot(mod)

plot_pred_cover(samples[2,], 'E:/ortho2018')
