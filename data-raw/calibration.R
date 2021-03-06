

annual2018_color_array <- car
save(annual2018_color_array, file = 'annual2018_color_array.rds')

scores_random18 <- samples2018[year == 2018 & id %in% 1:54]
scores_random18 <- unlist(sub[, c(4:53)])
site_ids <- as.character(unlist(lapply(1:54, function(x) rep(x, 50))))


library(muddier)
df <- data.frame(scores = scores_random18)
df <- cbind(df, car)

df$site_id <- site_ids

df$par <- 1
df$par[df$scores == 0] <- 0

df$full <- 1
df$full[df$scores != 2] <- 0

lm18 <- lm(scores ~ red_mn + red_sd + red_sk + red_kr + 
             grn_mn + 
             blu_mn + # blu_sd + blu_sk + blu_kr + 
             nir_mn + 
             ndvi_mn, data = df)
summary(lm18)

lm18 <- lm(scores ~ site_id + red_mn + red_sd + red_sk + red_kr + 
            grn_mn +
            blu_mn + blu_sd + blu_sk + blu_kr +
           nir_mn +
             ndvi_mn + ndvi_sd + ndvi_sk + ndvi_kr
           , data = df)
summary(lm18)

bin18 <- glm(par ~ site_id + #red_mn + red_sd + red_sk + red_kr + 
#            grn_mn + 
            # blu_mn + # blu_sd + blu_sk + blu_kr + 
            # nir_mn + 
            # ndvi_mn
, data = df, family = 'binomial')
summary(bin18)


# set library path
.libPaths( c( 'P:/lib', .libPaths()) )

setwd('E:/Riparian')
write("TMP = 'E:/Riparian/temp'", file=file.path(Sys.getenv('R_USER'), '.Renviron'))
library(riparian)
library(rgdal)
library(data.table)
library(magrittr)

samples <- readOGR('samples_nad83.shp')
usethis::use_data(samples)

res <- fread('samples2018.csv')                  #results of sampling observations

samples2018 <- res
colnames(samples2018) <- c('id','year','type',1:50)
usethis::use_data(samples2018, overwrite = T)


# 50-foot stream buffer for prc lots
prc_lots <- readOGR('riparian_lots.shp')
prc_strms <- readOGR('ripov_streams.shp')
prc_buff <- rgeos::gBuffer(streams, width = 50)
prc_per <- readOGR('srpo_strms.shp')

setwd('E:/Riparian/riparian')
usethis::use_data(prc_lots)
usethis::use_data(prc_strms)
usethis::use_data(prc_buff)
usethis::use_data(prc_per)

# permit list 2013-2018
permits_13to18 <- fread('permits_13to18.csv')
setwd('E:/Riparian/riparian')
usethis::use_data(permits_13to18)


# subset results by year
res <- samples2018

res18 <- res[year == 2018]
res16 <- res[year == 2016]
res09 <- res[year == 2009]

# divide into sampling groups

res18_p <- res18[id %in% 1:54]

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

ortho18 <- '/media/crumplecup/catacomb/gis/benton_2018'

mat <- matrix(c(rep(1,10), rep(2,10), rep(3,10)), ncol=3)
dt <- setDT(data.frame(mat))
unlist(dt)
dt18_p <- data.frame(rip = vals18_p) #%>% cbind(dt)

for (i in seq_along(samples))  {
  if (i == 1)  {
    dt <- get_cols_4band(samples[i,], ortho18)
  }
  if (i >1)  {
    dt <- rbind(dt, get_cols_4band(samples[i,], '/media/crumplecup/catacomb/gis/benton_2018'))
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

data(dt09, package = 'riparian')
dt_all <- dt_all %>% rbind(dt09)
usethis::use_data(dt_all, overwrite = T)

mod18 <- lm('rip ~ ndvi + red + grn + blu', data=dt_all)
summary(mod18)
usethis::use_data(mod18, overwrite = T)



sub <- samples2018[id %in% 1:54 & year == 2018, ] 
sub <- samples2018[id %in% 1:54 & year == 2016, ] 
ar <- array(t(sub[ , 4:53]), c(Reduce('*', dim(sub[, 4:53]))))
labs <- vector(length(ar), mode = 'character')
labs[ar == 0] <- 'bg'
labs[ar == 1] <- 'pc'
labs[ar == 2] <- 'fc'

setwd('E:/Riparian')
dir.create('rip_4band')
dir.create('rip_3band')
in_dir <- 'S:/maps/OrthoPhotos/Hexagon2016/RGB'
out_dir <- ('E:/Riparian/rip_3band/')
in_dir <- 'E:/ortho2018'
out_dir <- ('E:/Riparian/rip_4band/')


make_3band <- function(polys, keys, in_path, out_path, lab = 'id') {
  files <- get_rasters(in_path)
  polys <- match_crs(polys, raster::raster(file.path(in_path,files[1])))
  k <- 1
  for (i in seq_along(polys)) {
    boxs <- lapply(methods::slot(polys[i, ], 'polygons'),
                   function(x) methods::slot(x, 'Polygons'))[[1]]
    ras <- stack_extent(polys[i, ], in_path)
    ras <- raster::subset(ras, 1:3)
    
    for (j in seq_along(boxs)) {
      box <- spatialize(boxs[[j]], raster::crs(ras))
      m <- raster::crop(raster::mask(ras, box), raster::extent(box))
      nm <- paste(lab, i, j, keys[k], sep = '_')
      raster::writeRaster(m, paste0(out_path, nm), 'HFA')
      k <- k + 1
    }
  }
}

# make_3band(samples, labs, in_dir, out_dir, 'samples2018')



make_4band <- function(polys, keys, in_path, out_path, lab = 'id',
                       rgb_path = NULL, cir_path = NULL) {
  if (!is.null(rgb_path)) {
    rgb_files <- get_rasters(rgb_path)
    cir_files <- get_rasters(cir_path)
    polys <- match_crs(polys, raster::raster(file.path(rgb_path,rgb_files[1])))
    
  }
  if (is.null(cir_path)) {
    files <- get_rasters(in_path)
    polys <- match_crs(polys, raster::raster(file.path(in_path, files[1])))
  }
  k <- 1
  for (i in seq_along(polys)) {
    boxs <- lapply(methods::slot(polys[i,], 'polygons'),
                   function(x) methods::slot(x, 'Polygons'))[[1]]
    if (!is.null(rgb_path)) {
      rgb_ras <- stack_extent(polys[i, ], rgb_path)
      cir_ras <- raster::subset(stack_extent(polys[i, ], cir_path), 1)
      ras <- raster::stack(rgb_ras, cir_ras)
    }
    if (is.null(rgb_path)) {
      ras <- stack_extent(polys[i, ], in_path)
    }
    for (j in seq_along(boxs)) {
      box <- spatialize(boxs[[j]], raster::crs(ras))
      m <- raster::crop(raster::mask(ras, box), raster::extent(box))
      nm <- paste(lab, i, j, keys[k], sep = '_')
      raster::writeRaster(m, paste0(out_path, nm), 'HFA')
      k <- k + 1
    }
  }
}



# make_4band(samples, labs, in_dir, out_dir, 'samples2016')

sub <- samples2018[id %in% 1:54 & year == 2016, ] 
ar <- array(t(sub[ , 4:53]), c(Reduce('*', dim(sub[, 4:53]))))
labs <- vector(length(ar), mode = 'character')
labs[ar == 0] <- 'bg'
labs[ar == 1] <- 'pc'
labs[ar == 2] <- 'fc'

in_dir <- 'S:/maps/OrthoPhotos/Hexagon2016/RGB'
out_dir <- ('E:/Riparian/rip_3band/')

make_3band(samples, labs, in_dir, out_dir, 'samples2016')

sub <- samples2018[id %in% 1:54 & year == 2018, ] 
ar <- array(t(sub[ , 4:53]), c(Reduce('*', dim(sub[, 4:53]))))
labs <- vector(length(ar), mode = 'character')
labs[ar == 0] <- 'bg'
labs[ar == 1] <- 'pc'
labs[ar == 2] <- 'fc'

in_dir <- 'E:/ortho2018'
make_3band(samples, labs, in_dir, out_dir, 'samples2018')

out_dir <- ('E:/Riparian/rip_4band/')
make_4band(samples, labs, in_dir, out_dir, 'samples2018')


sub <- samples2018[id %in% 1:54 & year == 2016, ] 
ar <- array(t(sub[ , 4:53]), c(Reduce('*', dim(sub[, 4:53]))))
labs <- vector(length(ar), mode = 'character')
labs[ar == 0] <- 'bg'
labs[ar == 1] <- 'pc'
labs[ar == 2] <- 'fc'

rgb_dir <- 'S:/maps/OrthoPhotos/Hexagon2016/RGB'
cir_dir <- 'S:/maps/OrthoPhotos/Hexagon2016/CIR'

make_4band(samples, labs,
           rgb_path = rgb_dir,
           cir_path = cir_dir,
           out_path = out_dir, 
           lab = 'samples2016')


sub <- samples2018[id %in% 1:54 & year == 2009, ] 
ar <- array(t(sub[ , 4:53]), c(Reduce('*', dim(sub[, 4:53]))))
labs <- vector(length(ar), mode = 'character')
labs[ar == 0] <- 'bg'
labs[ar == 1] <- 'pc'
labs[ar == 2] <- 'fc'

in_dir <- 'E:/ortho2009'
out_dir <- ('E:/Riparian/rip_3band/')

make_3band(samples, labs, in_dir, out_dir, 'samples2009')


sub <- samples2018[year == 2018, ] 
sub <- sub[grep('p', sub$id), ]
p_ids <- 0
for (i in 1:nrow(sub)) {
  p_ids[i] <- as.numeric(strsplit(sub$id[i], 'p')[[1]][2])
}
p_sites <- readOGR('groupB_ortho2018.shp')
p_sites <- p_sites[p_ids, ]
ar <- array(t(sub[ , 4:53]), c(Reduce('*', dim(sub[, 4:53]))))
labs <- vector(length(ar), mode = 'character')
labs[ar == 0] <- 'bg'
labs[ar == 1] <- 'pc'
labs[ar == 2] <- 'fc'

in_dir <- 'E:/ortho2018'
out_dir <- ('E:/Riparian/rip_3band/')

make_3band(p_sites, labs, in_dir, out_dir, 'p_sites2018')

out_dir <- ('E:/Riparian/rip_4band/')
make_4band(p_sites, labs, in_dir, out_dir, 'p_sites2018')


sub <- samples2018[year == 2016, ] 
sub <- sub[grep('p', sub$id), ]
ar <- array(t(sub[ , 4:53]), c(Reduce('*', dim(sub[, 4:53]))))
labs <- vector(length(ar), mode = 'character')
labs[ar == 0] <- 'bg'
labs[ar == 1] <- 'pc'
labs[ar == 2] <- 'fc'

in_dir <- 'S:/maps/OrthoPhotos/Hexagon2016/RGB'
out_dir <- ('E:/Riparian/rip_3band/')

make_3band(p_sites, labs, in_dir, out_dir, 'p_sites2016')

out_dir <- ('E:/Riparian/rip_4band/')
make_4band(p_sites, labs, 
           rgb_path = rgb_dir,
           cir_path = cir_dir, 
           out_path = out_dir, 
           lab = 'p_sites2016')



sub <- samples2018[year == 2018, ] 
sub <- sub[grep('s', sub$id), ]
s_ids <- 0
for (i in 1:nrow(sub)) {
  s_ids[i] <- as.numeric(strsplit(sub$id[i], 's')[[1]][2])
}
s_sites <- readOGR('groupC_ortho2018.shp')
s_sites <- s_sites[s_ids, ]
ar <- array(t(sub[ , 4:53]), c(Reduce('*', dim(sub[, 4:53]))))
labs <- vector(length(ar), mode = 'character')
labs[ar == 0] <- 'bg'
labs[ar == 1] <- 'pc'
labs[ar == 2] <- 'fc'

in_dir <- 'E:/ortho2018'
out_dir <- ('E:/Riparian/rip_3band/')

make_3band(s_sites, labs, in_dir, out_dir, 's_sites2018')

out_dir <- ('E:/Riparian/rip_4band/')
make_4band(s_sites, labs, in_dir, out_dir, 's_sites2018')


sub <- samples2018[year == 2016, ] 
sub <- sub[grep('s', sub$id), ]
ar <- array(t(sub[ , 4:53]), c(Reduce('*', dim(sub[, 4:53]))))
labs <- vector(length(ar), mode = 'character')
labs[ar == 0] <- 'bg'
labs[ar == 1] <- 'pc'
labs[ar == 2] <- 'fc'

in_dir <- 'S:/maps/OrthoPhotos/Hexagon2016/RGB'
out_dir <- ('E:/Riparian/rip_3band/')

make_3band(s_sites, labs, in_dir, out_dir, 's_sites2016')

out_dir <- ('E:/Riparian/rip_4band/')
make_4band(s_sites, labs, 
           rgb_path = rgb_dir,
           cir_path = cir_dir, 
           out_path = out_dir, 
           lab = 's_sites2016')




base_dir <- 'E:/Riparian/ml'
dir.create(base_dir)

train_dir <- file.path(base_dir, 'train')
validation_dir <- file.path(base_dir, 'validation')
test_dir <- file.path(base_dir, 'test')
dir.create(train_dir)
dir.create(validation_dir)
dir.create(test_dir)

train_bg_dir <- file.path(train_dir, 'bg')
train_pc_dir <- file.path(train_dir, 'pc')
train_fc_dir <- file.path(train_dir, 'fc')
dir.create(train_bg_dir)
dir.create(train_pc_dir)
dir.create(train_fc_dir)

validation_bg_dir <- file.path(validation_dir, 'bg')
validation_pc_dir <- file.path(validation_dir, 'pc')
validation_fc_dir <- file.path(validation_dir, 'fc')
dir.create(validation_bg_dir)
dir.create(validation_pc_dir)
dir.create(validation_fc_dir)

test_bg_dir <- file.path(test_dir, 'bg')
test_pc_dir <- file.path(test_dir, 'pc')
test_fc_dir <- file.path(test_dir, 'fc')
dir.create(test_bg_dir)
dir.create(test_pc_dir)
dir.create(test_fc_dir)

in_dir <- 'E:/Riparian/rip_4band'
files <- list.files(in_dir) %>% sample
files_n <- length(files)
qtr_n <- ceiling(files_n / 4)
test_files <- files[1:qtr_n]
val_files <- files[(qtr_n + 1):(2*qtr_n)]
train_files <- files[(2*qtr_n + 1):files_n]

fnames <- train_files[grep('bg', train_files)]
file.copy(file.path(in_dir, fnames), file.path(train_bg_dir))
fnames <- train_files[grep('pc', train_files)]
file.copy(file.path(in_dir, fnames), file.path(train_pc_dir))
fnames <- train_files[grep('fc', train_files)]
file.copy(file.path(in_dir, fnames), file.path(train_fc_dir))

fnames <- val_files[grep('bg', val_files)]
file.copy(file.path(in_dir, fnames), file.path(validation_bg_dir))
fnames <- val_files[grep('pc', val_files)]
file.copy(file.path(in_dir, fnames), file.path(validation_pc_dir))
fnames <- val_files[grep('fc', val_files)]
file.copy(file.path(in_dir, fnames), file.path(validation_fc_dir))

fnames <- test_files[grep('bg', test_files)]
file.copy(file.path(in_dir, fnames), file.path(test_bg_dir))
fnames <- test_files[grep('pc', test_files)]
file.copy(file.path(in_dir, fnames), file.path(test_pc_dir))
fnames <- test_files[grep('fc', test_files)]
file.copy(file.path(in_dir, fnames), file.path(test_fc_dir))

library(keras)
datagen <- image_data_generator(rescale = 1/255)

train_generator <- flow_images_from_directory(
  train_dir,
  datagen,
  target_size = c(150, 150),
  batch_size = 20,
  class_mode = 'categorical'
)

validation_generator <- flow_images_from_directory(
  validation_dir,
  datagen,
  target_size = c(150, 150),
  batch_size = 20,
  class_mode = 'categorical'
)

batch <- generator_next(train_generator)

model <- keras_model_sequential() %>%
  layer_conv_2d(filters = 32, kernel_size = c(3, 3), activation = "relu",
                input_shape = c(150, 150, 4)) %>%
  layer_max_pooling_2d(pool_size = c(2, 2)) %>%
  layer_conv_2d(filters = 64, kernel_size = c(3, 3), activation = "relu") %>%
  layer_max_pooling_2d(pool_size = c(2, 2)) %>%
  layer_conv_2d(filters = 128, kernel_size = c(3, 3), activation = "relu") %>%
  layer_max_pooling_2d(pool_size = c(2, 2)) %>%
  layer_conv_2d(filters = 128, kernel_size = c(3, 3), activation = "relu") %>%
  layer_max_pooling_2d(pool_size = c(2, 2)) %>%
  layer_flatten() %>%
  layer_dense(units = 512, activation = "relu") %>%
  layer_dense(units = 3, activation = "softmax")
summary(model)

model %>% compile(
  loss = "categorical_crossentropy",
  optimizer = optimizer_rmsprop(lr = 1e-4),
  metrics = c("acc")
)

history <- model %>% fit_generator(
  train_generator,
  steps_per_epoch = 100,
  epochs = 4,
  validation_data = validation_generator,
  validation_steps = 50
)



in_dir <- 'E:/ortho2018'

tensor_4band <- function(polys, keys, in_path, frame = c(55, 55)) {
  files <- get_rasters(in_path)
  ar <- array(1, c(length(labs), frame, 4))
  k <- 1
  for (i in seq_along(polys)) {
    poly <- match_crs(samples[1,], raster::raster(file.path(in_path,files[1])))
    boxs <- lapply(methods::slot(poly, 'polygons'),
                   function(x) methods::slot(x, 'Polygons'))[[1]]
    ras <- stack_extent(poly, in_path)
    
    for (j in seq_along(boxs)) {
      box <- spatialize(boxs[[j]], raster::crs(ras))
      m <- raster::crop(raster::mask(ras, box), raster::extent(box))
      for (b in 1:4) {
        r <- raster::raster(m, layer = b)
        mr <- matrix(r, ncol = ncol(r), byrow = TRUE)
        mr[is.na(mr)] <- 255
        mr <- mr / 255
        ar[k, 1:nrow(mr), 1:ncol(mr), b] <- mr
      }
      k <- k + 1
    }
  }
  return(list(ar, keys))
}

rip_4band <- tensor_4band(samples, labs, in_dir)

model <- keras_model_sequential() %>%
  layer_conv_2d(filters = 32, kernel_size = c(3, 3), activation = "relu",
                input_shape = c(55, 55, 4)) %>%
  layer_max_pooling_2d(pool_size = c(2, 2)) %>%
  layer_conv_2d(filters = 64, kernel_size = c(3, 3), activation = "relu") %>%
  layer_max_pooling_2d(pool_size = c(2, 2)) %>%
  layer_conv_2d(filters = 128, kernel_size = c(3, 3), activation = "relu") %>%
  layer_max_pooling_2d(pool_size = c(2, 2)) %>%
  layer_conv_2d(filters = 128, kernel_size = c(3, 3), activation = "relu") %>%
  layer_max_pooling_2d(pool_size = c(2, 2)) %>%
  layer_flatten() %>%
  layer_dense(units = 512, activation = "relu") %>%
  layer_dense(units = 3, activation = "softmax")
summary(model)

model %>% compile(
  loss = "categorical_crossentropy",
  optimizer = optimizer_rmsprop(lr = 1e-4),
  metrics = c("acc")
)

rip_4band[[2]] <- ar

files_n <- length(rip_4band[[2]])
qtr_n <- ceiling(files_n / 4)
test_ids <- 1:qtr_n
val_ids <- (qtr_n + 1):(2*qtr_n)
train_ids <- (2*qtr_n + 1):files_n

gen <- function(index, batch_size = 20) {
  ids <- sample(index, batch_size)
  samps <- rip_4band[[1]]
  keys <- rip_4band[[2]]
  list(samps[ids, , , ], keys[ids])
}

batch <- gen(train_ids)

sum(!rip_4band[[2]] == labs)

history <- model %>% fit_generator(
  gen(train_ids),
  steps_per_epoch = 100,
  epochs = 30,
  validation_data = gen(val_ids),
  validation_steps = 50
)

id_bin <- to_categorical(rip_4band[[2]])
net <- model %>% fit(rip_4band[[1]][1:1200,,,], id_bin[1:1200,], 
                     epochs = 10, batch_size = 256)


setwd('E:/Riparian/riparian')
usethis::use_data(rip_4band)




sub <- samples2018[year == 2018, ] 
sub <- sub[grep('o', sub$id), ]
oids <- sub$id[-55]

sub <- samples2018[year == 2018, ] 
sub <- sub[grep('s', sub$id), ]
sids <- sub$id
s_ids <- 0
for (i in 1:nrow(sub)) {
  s_ids[i] <- as.numeric(strsplit(sub$id[i], 's')[[1]][2])
}
s_sites <- readOGR('groupC_ortho2018.shp')
s_sites <- s_sites[s_ids, ]

sub <- samples2018[year == 2018, ] 
sub <- sub[grep('p', sub$id), ]
pids <- sub$id
p_ids <- 0
for (i in 1:nrow(sub)) {
  p_ids[i] <- as.numeric(strsplit(sub$id[i], 'p')[[1]][2])
}
p_sites <- readOGR('groupB_ortho2018.shp')
p_sites <- p_sites[p_ids, ]

olen <- length(oids)
plen <- length(pids)

crs_ref <- raster::crs(samples)
p_sites <- spTransform(p_sites, crs_ref)
s_sites <- spTransform(s_sites, crs_ref)

op <- slot(samples, 'polygons')
pp <- slot(p_sites, 'polygons')
sp <- slot(s_sites, 'polygons')
ss <- list()

for (i in seq_along(op)) {
  slot(op[[i]], 'ID') <- oids[i]
  ss[[i]] <- op[[i]]
}

for (i in seq_along(pp)) {
  slot(pp[[i]], 'ID') <- pids[i]
  ss[[i + olen]] <- pp[[i]]
}

for (i in seq_along(sp)) {
  slot(sp[[i]], 'ID') <- sids[i]
  ss[[i + olen + plen]] <- sp[[i]]
}

ids <- data.frame(ID = c(oids, pids, sids))
rownames(ids) <- c(oids, pids, sids)

polys <- SpatialPolygons(ss, proj4string = crs_ref)
supersample <- SpatialPolygonsDataFrame(polys, data = ids)
setwd('E:/Riparian/riparian')
usethis::use_data(supersample, overwrite = T)

thumdir <- 'E:/Riparian/thumbnails'
thumbnails('E:/ortho2018', thumdir, supersample)

setwd(thumdir)
file.rename(list.files(getwd()), 
            paste0('samples2018_', 1:length(supersample)))

thumbnails('S:/maps/OrthoPhotos/Hexagon2016/RGB', thumdir, supersample)
file.rename(list.files(getwd()), 
            paste0('rgb2016_', 1:length(supersample)))

thumbnails('S:/maps/OrthoPhotos/Hexagon2016/CIR', thumdir, supersample)
file.rename(list.files(getwd()), 
            paste0('cir2016_', 1:length(supersample)))

thumbnails('E:/ortho2009', thumdir, supersample)
file.rename(list.files(getwd()), 
            paste0('samples2009_', 1:length(supersample)))

setwd('/home/crumplecup/work')
load('rip_4band_df.rds')
load('rip_lmod.rds')
rip_lmod <- lm(as.character(lmod[1,1]), data = df)

load('rip_lmod3.rds')
rip_lmod3 <- lm(as.character(lmod[1,1]), data = df)

load('rip_bin1.rds')
dt <- df
dt$cov1 <- dt$cov
dt$cov1[dt$cov1 == 2] <- 1
dt <- dt[ , c(17,2:16)]
rip_bin1 <- glm(as.character(bin1[1,1]), family = 'binomial', data = dt)
load('rip_bin2.rds')
dt <- df
dt$cov1 <- dt$cov - 1
dt <- dt[dt$cov > 0, ]
dt <- dt[ , c(17,2:16)]
rip_bin2 <- glm(as.character(bin2[1,1]), family = 'binomial', data = dt)

setwd('/home/crumplecup/work/riparian')
usethis::use_data(rip_lmod)
usethis::use_data(rip_lmod3)
usethis::use_data(bad_box)


setwd('/home/crumplecup/work')
load('sam.rds')
bad_box <- sam[15, ]


rasdir <- '/media/crumplecup/Seagate Backup Plus Drive/gis/ortho2018'
slcdir <- '/home/crumplecup/work/slices'
plot_samples(rasdir, slcdir, sam[1:14,], 2018, 'random', 'lm')



