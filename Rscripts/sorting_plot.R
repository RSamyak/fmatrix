
stop() ## load sorting_script.R first

normF <- function(Fmat, p = 1){
  sum(abs(Fmat**p))**(1/p)
}


exttreelength <- function(Fmat){
  tree <- mytree_from_F(Fmat)

  n <- ncol(Fmat) + 1

  external <- temp$edge[,1] <= n | temp$edge[,2] <= n

  sum(tree$edge.length[external])
}



kingman_prob <- sapply(F.list, beta_BF_likelihood, betas = 0)

l1kingman_mean <- F.list[[brute.mean.weighted(F.list, F.list, dist = "l1", weights = kingman_prob)]]
l2kingman_mean <- F.list[[brute.mean.weighted(F.list, F.list, dist = "l2", weights = kingman_prob)]]



encodkingman_mean <- F.list[[brute.mean.weighted(encod.list, encod.list, dist = "l1", weights = kingman_prob)]]

temp_encod <- my_encod(encodkingman_mean)

sorted.encoddist2kingmean.row <- sapply(sorted.encod.list.row, function(x){sum(abs(x - temp_encod)/2**(-(1:length(x))))})
sorted.encoddist2kingmean.col <- sapply(sorted.encod.list.col, function(x){sum(abs(x - temp_encod)/2**(-(1:length(x))))})


######


which(check_same)

plotF.list(sorted.list.row[1:7])
plotF.list(sorted.list.col[1:7])
plotF.list(sorted.list.row[16:21])
plotF.list(sorted.list.col[16:21])
plotF.list(sorted.list.row[(479:488)[-3]])
plotF.list(sorted.list.col[(479:488)[-3]])




######
sorted.norm.row <- sapply(sorted.list.row, normF)
sorted.norm.col <- sapply(sorted.list.col, normF)


sorted.norm2.row <- sapply(sorted.list.row, normF, p = 2)
sorted.norm2.col <- sapply(sorted.list.col, normF, p = 2)


sorted.numcherries.row <- sapply(sorted.list.row, num_cherries)
sorted.numcherries.col <- sapply(sorted.list.col, num_cherries)


sorted.exttreelength.row <- sapply(sorted.list.row, exttreelength)
sorted.exttreelength.col <- sapply(sorted.list.col, exttreelength)


sorted.l1dist2kingmean.row <- sapply(sorted.list.row, distance_Fmat, Fmat2 = l1kingman_mean, dist = "l1")
sorted.l1dist2kingmean.col <- sapply(sorted.list.col, distance_Fmat, Fmat2 = l1kingman_mean, dist = "l1")

sorted.l2dist2kingmean.row <- sapply(sorted.list.row, distance_Fmat, Fmat2 = l2kingman_mean, dist = "l2")
sorted.l2dist2kingmean.col <- sapply(sorted.list.col, distance_Fmat, Fmat2 = l2kingman_mean, dist = "l2")


dev.off()

if(! dir.exists("plots")) dir.create("plots")

png("plots/norm.png", width = 2000, height = 1000)
par(mfrow = c(1,2))
plot(sorted.norm.row)
points(which(check_same), sorted.norm.row[check_same], col = "red")
plot(sorted.norm.col)
points(which(check_same), sorted.norm.col[check_same], col = "red")
dev.off()

png("plots/norm2.png", width = 2000, height = 1000)
par(mfrow = c(1,2))
plot(sorted.norm2.row)
points(which(check_same), sorted.norm2.row[check_same], col = "red")
plot(sorted.norm2.col)
points(which(check_same), sorted.norm2.col[check_same], col = "red")
dev.off()


png("plots/numcherries.png", width = 2000, height = 1000)
par(mfrow = c(1,2))
plot(sorted.numcherries.row)
points(which(check_same), sorted.numcherries.row[check_same], col = "red")
plot(sorted.numcherries.col)
points(which(check_same), sorted.numcherries.col[check_same], col = "red")
dev.off()

png("plots/exttreelength.png", width = 2000, height = 1000)
par(mfrow = c(1,2))
plot(sorted.exttreelength.row)
points(which(check_same), sorted.exttreelength.row[check_same], col = "red")
plot(sorted.exttreelength.col)
points(which(check_same), sorted.exttreelength.col[check_same], col = "red")
dev.off()

png("plots/l1dist2kingmean.png", width = 2000, height = 1000)
par(mfrow = c(1,2))
plot(sorted.l1dist2kingmean.row)
points(which(check_same), sorted.l1dist2kingmean.row[check_same], col = "red")
plot(sorted.l1dist2kingmean.col)
points(which(check_same), sorted.l1dist2kingmean.col[check_same], col = "red")
dev.off()

png("plots/l2dist2kingmean.png", width = 2000, height = 1000)
par(mfrow = c(1,2))
plot(sorted.l2dist2kingmean.row)
points(which(check_same), sorted.l2dist2kingmean.row[check_same], col = "red")
plot(sorted.l2dist2kingmean.col)
points(which(check_same), sorted.l2dist2kingmean.col[check_same], col = "red")
dev.off()

png("plots/enco/ddist2kingmean.png", width = 2000, height = 1000)
par(mfrow = c(1,2))
plot(sorted.encoddist2kingmean.row)
points(which(check_same), sorted.encoddist2kingmean.row[check_same], col = "red")
plot(sorted.encoddist2kingmean.col)
points(which(check_same), sorted.encoddist2kingmean.col[check_same], col = "red")
dev.off()

