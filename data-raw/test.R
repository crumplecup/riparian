
ortho20 <- '/media/crumplecup/catacomb/gis/benton20_3in'
ortho20b <- '/media/crumplecup/catacomb/gis/benton20_9in'
ortho18 <- '/media/crumplecup/catacomb/gis/benton_2018'

sam_path <- '/home/crumplecup/work/samples'

# annual thumbnail paths
thumb_20 <- file.path(sam_path, 'annual2020/thumb20_3in')
thumb_20b <- file.path(sam_path, 'annual2020/thumb20_9in')
thumb_18 <- file.path(sam_path, 'annual2020/thumb18')

# plot output paths
out_20 <- file.path(sam_path, 'permit2020/2020_3in')
out_20b <- file.path(sam_path, 'permit2020/2020_9in')
out_18 <- file.path(sam_path, 'permit2020/2018')

# write thumbnails to drive
thumbnails(ortho20, thumb_20, random_samples)
thumbnails(ortho20b, thumb_20b, random_samples)
thumbnails(ortho18, thumb_18, random_samples)

# plot sample boxes
plot_samples(thumb_20, out_20, random_samples, 2020, method = 'lm3')
plot_samples(thumb_20b, out_20b, random_samples, 2020, method = 'lm3')
plot_samples(thumb_18, out_18, random_samples, 2018, method = 'lm')

library(data.table)
perms <- fread('permits20.csv')
perm_mtls <- perms$maptaxlot
perm_lots <- sf_lots[sf_lots$MapTaxlot %in% perm_mtls, ]

permit_samples <- sample_streams(75, lots = perm_lots)
# screen out bad boxes [omitted]
# save(permit_samples, file = 'permit_samples.rds')

# thumbnail paths
thumb20 <- file.path(sam_path, 'permit2020/2020/thumbs_3in')
thumb20b <- file.path(sam_path, 'permit2020/2020/thumbs_9in')
thumb18 <- file.path(sam_path, 'permit2020/2018/thumbs')

# sample plot paths
out18 <- file.path(sam_path, 'permit2020/2018')
out20 <- file.path(sam_path, 'permit2020/2020/3in')
out20b <- file.path(sam_path, 'permit2020/2020/9in')

# write thumbnail tifs to drive
thumbnails(ortho20, thumb20, permit_samples)
thumbnails(ortho20b, thumb20b, permit_samples)
thumbnails(ortho18, thumb18, permit_samples)

# plot samples (also produces csv)
# missing samples for the 3in ortho, use the 9in csv
plot_samples(thumb20, out20, permit_samples, 2020, method = 'lm3')
plot_samples(thumb20b, out20b, permit_samples, 2020, method = 'lm3')
plot_samples(thumb18, out18, permit_samples, 2018, method = 'lm')


test <- file.path(sam_path, 'permit2020/test')
setwd('/home/crumplecup/work/samples/permit2020')
setwd(test)

# correct auto-filled csvs before running
obs20 <- fread('samples_2020.csv')
obs20 <- obs20[, -1]
obs18 <- fread('samples_2018.csv')
obs18 <- obs18[, -1]

# plot scores (already done by plot_samples)
# plot_scores(thumb18, test, obs18, permit_samples, title = 'permits18_')
build_change_table(rbind(obs18, obs20), 2018, 2020, permit_samples, perm_lots, perms)

# plot % cover, prints stats to console
plot_cover(obs18, title = 'permits20_cover18')
plot_cover(obs20, title = 'permits20_cover20')

# plot % change, prints stats to console
plot_change(rbind(obs18, obs20), 2018, 2020, title = 'permits20_change18to20')

# use to print thumbnails of permit taxlots
# helpful to determine if tree removal occurs inside a lot with a permit
thumbnails(ortho20, test, perm_lots, 'png')

# 4-band data needed (not available for 2020)
# pred18 <- file.path(test, 'pred_cover18')
# pred_cover(perm_lots, in_path = ortho18, pred18)
# pred_cover_report(pred18)

# for testing load the 2009 data
# pred09 <- file.path(test, 'pred_cover09')
# pred_cover(perm_lots, in_path = '/media/crumplecup/catacomb/gis/ortho2009/', pred09)
# pred_cover_report(pred09)

# chng_path <- file.path(test, 'chng_09to18')
# pred_change(pred09, pred18, chng_path)
# pred_change_report(chng_path, pred09, pred18, perms)
