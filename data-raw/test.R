
inspect_samples(permit_samples, ortho20, '/home/erik/output/inspect/', permits = perm_sub)

blen <- 10
bit <- rep(5, blen)
tbit <- stats::fft(bit) / blen
stats::fft(tbit, inverse = T)
byt <- 1:blen
tbyt <- stats::fft(byt) / blen
ibyt <- stats::fft(tbyt, inverse = T) / blen
tb <- tbyt + tbit
itb <- stats::fft(tb, inverse = T)
bat <- c(5, 5, 5, 5, 5)
tbat <- stats::fft(bat) / 5
tb <- tbit * tbat
itb <- stats::fft(tb, inverse = T) / 5
plot(abs(itb))

rec <- data.table::fread('/home/erik/output/df20k_1000.csv')
recf <- data.table::fread('/home/erik/output/ff20k_1000.csv')
recg <- data.table::fread('/home/erik/output/fg20k_1000.csv')
sum(nrow(rec) + nrow(recf) + nrow(recg))
# AD test
plot(rec$input, rec$ad, pch = 20, col = get_palette('charcoal'),
     xlab = 'rate', ylab = 'test statistic')
points(recf$input, recf$ad, col = get_palette('crimson'), pch = 20)
points(recg$input, recg$ad, col = get_palette('ocean'), pch = 20)

# Chi-squared
plot(rec$input, rec$chs, pch = 20, col = get_palette('charcoal'),
     xlab = 'rate', ylab = 'test statistic')
points(recf$input, recf$chs, col = get_palette('crimson'), pch = 20)
points(recg$input, recg$chs, col = get_palette('ocean'), pch = 20)

# KS test
plot(rec$input, rec$ks, pch = 20, col = get_palette('charcoal'),
     xlab = 'rate', ylab = 'test statistic')
points(recg$input, recg$ks, col = get_palette('ocean'), pch = 20)
points(recf$input, recf$ks, col = get_palette('crimson'), pch = 20)
legend('topright', legend = c('debris flows', 'fines', 'gravels'), 
       fill = get_palette(c('charcoal', 'crimson', 'ocean'), .8))

# KS test sub
subd <- rec[rec$ks < .20, ]
subf <- recf[recf$ks < .20, ]
subg <- recg[recg$ks < .20, ]

plot(subd$input, subd$ks, pch = 20, col = get_palette('charcoal'),
     xlab = 'rate', ylab = 'test statistic')
points(subg$input, subg$ks, col = get_palette('ocean'), pch = 20)
points(subf$input, subf$ks, col = get_palette('crimson'), pch = 20)
legend('topright', legend = c('debris flows', 'fines', 'gravels'), 
       fill = get_palette(c('charcoal', 'crimson', 'ocean'), .8))

# Kuiper test
plot(rec$input, rec$kp, pch = 20, col = get_palette('charcoal'),
     xlab = 'rate', ylab = 'test statistic')
points(recg$input, recg$kp, col = get_palette('ocean'), pch = 20)
points(recf$input, recf$kp, col = get_palette('crimson'), pch = 20)


d <- 0.1 / 0.237
pwr <- seq(0, .99, .01)
n <- unlist(lapply(pwr, function(x) power.t.test(d= d, power = x)$n))

png('monitoring_stat_pwr.png', height = 17, width = 25, units = 'cm', res = 300)
plot(n, pwr, type = 'l', col = get_palette('ocean'), lwd = 3,
     xlab = '# of samples', ylab = 'statistical power',
     main = 'power to detect mean cover change of 10% in priority corridors')
abline(v = 54, lty = 2, col = get_palette('coral'), lwd = 2)
pts <- matrix(c(54, 0.59), ncol = 2)
abline(v = 89, lty = 2, col = get_palette('charcoal'), lwd = 2)
pts <- rbind(pts, c(89, 0.8))
abline(v = 119, lty = 2, col = get_palette('gold'), lwd = 2)
pts <- rbind(pts, c(119, 0.9))
text(60, .59, '54')
text(95, .8, '89')
text(126, .9, '119')
points(pts, col = get_palette(c('coral', 'charcoal', 'gold'), .7))
legend('bottomright', legend = c('59% power (2018 level)', '80% power', '90% power'),
       fill = get_palette(c('coral', 'charcoal', 'gold'), .7))
dev.off()

plot(sf::st_zm(random_samples[33,]))



load('/home/erik/data/permit_samples_144_final.rds')
load('/home/erik/data/perm_sub.rds')
perm_lots <- sf_lots[sf_lots$MapTaxlot %in% perm_sub$maptaxlot, ]
plot(perm_lots[, 1])

ortho20 <- '/media/erik/catacomb/gis/benton20_3in'
ortho20b <- '/media/crumplecup/catacomb/gis/benton20_9in'
ortho18 <- '/media/crumplecup/catacomb/gis/benton_2018'

sam_path <- '/home/erik/output/thumbs'
# dir.create('/home/erik/output/perm_plots')
perm_path <- '/home/erik/output/perm_plots'

# annual thumbnail paths
thumb_20 <- file.path(sam_path, 'annual2020/thumb20_3in')
thumb_20b <- file.path(sam_path, 'annual2020/thumb20_9in')
thumb_18 <- file.path(sam_path, 'annual2020/thumb18')

# plot output paths
out_20 <- file.path(sam_path, 'permit2020/2020_3in')
out_20b <- file.path(sam_path, 'permit2020/2020_9in')
out_18 <- file.path(sam_path, 'permit2020/2018')

# write thumbnails to drive
thumbnails(ortho20, sam_path, permit_samples)
thumbnails(ortho20b, thumb_20b, random_samples)
thumbnails(ortho18, thumb_18, random_samples)

# plot sample boxes
plot_samples(sam_path, '/home/erik/output', permit_samples, 2020, method = 'lm3')
plot_samples(thumb_20, out_20, random_samples, 2020, method = 'lm3')
plot_samples(thumb_20b, out_20b, random_samples, 2020, method = 'lm3')
plot_samples(thumb_18, out_18, random_samples, 2018, method = 'lm')

library(data.table)
perms <- data.table::fread('/home/erik/data/permits20.csv')
perm_mtls <- perms$maptaxlot
perm_lots <- sf_lots[sf_lots$MapTaxlot %in% perm_mtls, ]

perm_mtls <- perm_sub$maptaxlot
perm_lots <- sf_lots[sf_lots$MapTaxlot %in% perm_mtls, ]
class(perm_lots)
thumbnails(ortho20, '/home/erik/output', perm_lots, output = 'png')


permit_samples <- sample_streams(3, lots = perm_lots)
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
