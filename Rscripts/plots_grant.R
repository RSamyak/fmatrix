
library(fmatrix)
data("simulations-data")
#means_list_large_n <- lapply(iso.Fmat.list, sa_mean, temp.schedule = "exp", alpha = .9995)
#png("means_list_large_n.png", height = 500, width = 2500)
#plotF.list(means_list_large_n[length(means_list_large_n):1], ncol = 5)
#temp <- phylodyn:::gen_Fmat(rcoal(8, br = 1))
#plotF(temp, node.labels = T,tip.labels = T)
#mymatching(mytree_from_F(temp))

betas = c(seq(-1,0, by = .2), seq(0, 1, by = .2)[-1], seq(1, 10, by = 1,)[-1], seq(10, 100, by = 10)[-1])

F.lists <- F.list9
kingprob <- sapply(F.lists, kingman_likelihood)
betaprob <- sapply(F.lists, beta_BF_likelihood, betas = betas)

kingtemp <- fmatrix:::frechet_var_pop(F.lists, kingprob)
  kingtemp$meandist
  kingtemp$var
  kingtemp$fmean
  #plotF(kingtemp$fmean,edge.width=2)

  dist.kingman<-rep(0,length(F.lists))
  for (j in 1:length(F.lists)){
    dist.kingman[j]<-sqrt(sum((F.lists[[j]]-kingtemp$fmean)^2))
  }

  dist.caterpilar<-rep(0,length(F.lists))
  dist.balanced<-rep(0,length(F.lists))

  for (j in 1:length(F.lists)){
    dist.caterpilar[j]<-sqrt(sum((F.lists[[j]]-F.lists[[length(F.lists)]])^2))
  }
  balanced.index<-which.max(dist.caterpilar)

  for (j in 1:length(F.lists)){
    dist.balanced[j]<-sqrt(sum((F.lists[[j]]-F.lists[[balanced.index]])^2))
  }

  dist.kingman.sign<-dist.kingman
  for (j in 1:length(dist.kingman.sign)){
    if (dist.balanced[j]<dist.caterpilar[j]){
      dist.kingman.sign[j]<--dist.kingman.sign[j]
    }
  }

  #dist.kingman<-dist.kingman/max(dist.kingman)

  dist.kingman.sort<-sort(dist.kingman.sign,index.return=TRUE)



#####Figure for paper
par(mfrow=c(3,5))
plot(seq(1,length(dist.kingman)),kingprob[dist.kingman.sort$ix],main="Tajima/Yule",type="S",ylab="Probability",xlab="Trees")
bb<-.5
plot(seq(1,length(dist.kingman)),exp(-bb*sort(dist.kingman.sign)^2)/sum(exp(-bb*sort(dist.kingman.sign)^2)),main="Exp-Dist center Tajima",xlab="Trees",ylab="Probability",type="l")
plot(sort(dist.kingman.sign),exp(-bb*sort(dist.kingman.sign)^2)/sum(exp(-bb*sort(dist.kingman.sign)^2)),main="Exp-Dist center Tajima",xlab="Signed Distance",ylab="Probability",type="l")

hist.king<-tapply(kingprob, dist.kingman, FUN=sum)
kb<-round(as.numeric(names(hist.king)),2)
names(hist.king)<-kb
barplot(hist.king,main="Tajima, n=9",ylab="Probability",xlab="Distance to the Mean")

hist.king2<-tapply(kingprob, dist.kingman.sign, FUN=sum)
kb<-round(as.numeric(names(hist.king2)),2)
names(hist.king2)<-kb
barplot(hist.king2,main="Tajima, n=9",ylab="Probability",xlab="Signed Distance")

j<-2
temp <- fmatrix:::frechet_var_pop(F.lists, betaprob[j, ])

temp$meandist
temp$var
temp$fmean


dist.beta<-rep(0,length(F.lists))
for (j in 1:length(F.lists)){
  dist.beta[j]<-sqrt(sum((F.lists[[j]]-temp$fmean)^2))
}

dist.beta.sign<-dist.beta
for (j in 1:length(dist.beta.sign)){
  if (dist.balanced[j]<dist.caterpilar[j]){
    dist.beta.sign[j]<--dist.beta.sign[j]
  }
}

dist.beta.sort<-sort(dist.beta.sign,index.return=TRUE)
j<-2
plot(seq(1,length(dist.beta)),betaprob[j,dist.beta.sort$ix],main="Beta b=-0.8",type="S",ylab="Probability",xlab="Trees")

bb<-.5
plot(seq(1,length(dist.beta)),exp(-bb*sort(dist.beta.sign)^2)/sum(exp(-bb*sort(dist.beta.sign)^2)),main="Exp-Dist center Beta b=-.4",xlab="Trees",ylab="Probability",type="l")

plot(sort(dist.beta.sign),exp(-bb*sort(dist.beta.sign)^2)/sum(exp(-bb*sort(dist.beta.sign)^2)),main="Exp-Dist center Beta b=-.4",xlab="Signed Distance",ylab="Probability",type="l")

hist.beta1<-tapply(betaprob[j,],dist.beta,FUN=sum)
xb<-round(as.numeric(names(hist.beta1)),2)
names(hist.beta1)<-xb
barplot(hist.beta1,main="Beta b=-0.8, n=9",ylab="Probability",xlab="Distance to the Mean")

hist.beta2<-tapply(betaprob[j,],dist.beta.sign,FUN=sum)
xb<-round(as.numeric(names(hist.beta2)),2)
names(hist.beta2)<-xb
barplot(hist.beta2,main="Beta b=-0.8, n=9",ylab="Probability",xlab="Singed Distance")


j<-29
temp <- fmatrix:::frechet_var_pop(F.lists, betaprob[j, ])

temp$meandist
temp$var
temp$fmean


dist.beta<-rep(0,length(F.lists))
for (j in 1:length(F.lists)){
  dist.beta[j]<-sqrt(sum((F.lists[[j]]-temp$fmean)^2))
}

dist.beta.sign<-dist.beta
for (j in 1:length(dist.beta.sign)){
  if (dist.balanced[j]<dist.caterpilar[j]){
    dist.beta.sign[j]<--dist.beta.sign[j]
  }
}

dist.beta.sort<-sort(dist.beta.sign,index.return=TRUE)
j<-29
plot(seq(1,length(dist.beta)),betaprob[j,dist.beta.sort$ix],main="Beta b=100",type="S",ylab="Probability",xlab="Trees")


bb<-.5
plot(seq(1,length(dist.beta)),exp(-bb*sort(dist.beta.sign)^2)/sum(exp(-bb*sort(dist.beta.sign)^2)),main="Exp-Dist center Beta b=100",xlab="Trees",ylab="Probability",type="l")
plot(sort(dist.beta.sign),exp(-bb*sort(dist.beta.sign)^2)/sum(exp(-bb*sort(dist.beta.sign)^2)),main="Exp-Dist center Beta b=100",xlab="Signed Distance",ylab="Probability",type="l")

hist.beta1<-tapply(betaprob[j,],dist.beta,FUN=sum)
xb<-round(as.numeric(names(hist.beta1)),2)
names(hist.beta1)<-xb
barplot(hist.beta1,main="Beta b=100, n=9",ylab="Probability",xlab="Distance to the Mean")

hist.beta2<-tapply(betaprob[j,],dist.beta.sign,FUN=sum)
xb<-round(as.numeric(names(hist.beta2)),2)
names(hist.beta2)<-xb
barplot(hist.beta2,main="Beta b=-0.8, n=9",ylab="Probability",xlab="Singed Distance")





#####Figure for grant

par(mfrow=c(2,4))
n <- ncol(kingtemp$fmean) + 1
tree<-tree_from_F(rbind(rep(0,n),cbind(rep(0,n-1),kingtemp$fmean)),seq(1,n-1))
ape::plot.phylo(ape::ladderize(tree), direction = "downwards", show.tip.label = FALSE, edge.width=2)


plot(seq(1,length(dist.kingman)),sort(kingprob),main="Tajima/Yule",type="S",ylab="Probability",xlab="",cex.lab=1.6,cex.main=2)
#plot(seq(1,length(dist.kingman)),kingprob[dist.kingman.sort$ix],main="Tajima/Yule",type="S",ylab="Probability",xlab="Trees")
bb<-.2
plot(seq(1,length(dist.kingman)),exp(-bb*sort(dist.kingman.sign)^2)/sum(exp(-bb*sort(dist.kingman.sign)^2)),main="Exp-Dist",xlab="",ylab="Probability",type="l",cex.lab=1.6,cex.main=2)


hist.king2<-tapply(kingprob, dist.kingman, FUN=sum)
kb<-round(as.numeric(names(hist.king2)),2)
names(hist.king2)<-kb
barplot(hist.king2,main="Tajima, n=9",ylab="Probability",xlab="",cex.lab=1.6,cex.main=2)



j<-29
temp <- fmatrix:::frechet_var_pop(F.lists, betaprob[j, ])

temp$meandist
temp$var
temp$fmean

tree<-tree_from_F(rbind(rep(0,n),cbind(rep(0,n-1),temp$fmean)),seq(1,n-1))
ape::plot.phylo(ape::ladderize(tree), direction = "downwards", show.tip.label = FALSE, edge.width=2)

dist.beta<-rep(0,length(F.lists))
for (j in 1:length(F.lists)){
  dist.beta[j]<-sqrt(sum((F.lists[[j]]-temp$fmean)^2))
}

dist.beta.sign<-dist.beta
for (j in 1:length(dist.beta.sign)){
  if (dist.balanced[j]<dist.caterpilar[j]){
    dist.beta.sign[j]<--dist.beta.sign[j]
  }
}

dist.beta.sort<-sort(dist.beta.sign,index.return=TRUE)
j<-29
plot(seq(1,length(dist.beta)),sort(betaprob[j,]),main="Beta b=100",type="S",ylab="Probability",xlab="Trees",cex.lab=1.6,cex.main=2)

#plot(seq(1,length(dist.beta)),betaprob[j,dist.beta.sort$ix],main="Beta b=100",ylab="Probability",xlab="Trees")

bb<-.2
plot(seq(1,length(dist.beta)),exp(-bb*sort(dist.beta.sign)^2)/sum(exp(-bb*sort(dist.beta.sign)^2)),main="Exp-Dist",xlab="Trees",ylab="Probability",type="l",cex.lab=1.6,cex.main=2)


hist.beta2<-tapply(betaprob[j,],dist.beta,FUN=sum)
xb<-round(as.numeric(names(hist.beta2)),2)
names(hist.beta2)<-xb
barplot(hist.beta2,main="Beta b=100, n=9",ylab="Probability",xlab="Distance to Mean",cex.lab=1.6,cex.main=2)



##more attemps, this time sorting always the same way

par(mfrow=c(3,3))
plot(seq(1,length(dist.kingman)),kingprob[dist.kingman.sort$ix],main="Tajima/Yule",type="S",ylab="Probability",xlab="Trees")
bb<-.5
plot(seq(1,length(dist.kingman)),exp(-bb*sort(dist.kingman.sign)^2)/sum(exp(-bb*sort(dist.kingman.sign)^2)),main="Exp-Dist center Tajima",xlab="Trees",ylab="Probability",type="l")


hist.king2<-tapply(kingprob, dist.kingman.sign, FUN=sum)
kb<-round(as.numeric(names(hist.king2)),2)
names(hist.king2)<-kb
barplot(hist.king2,main="Tajima, n=9",ylab="Probability",xlab="Signed Distance")

j<-2
temp <- fmatrix:::frechet_var_pop(F.lists, betaprob[j, ])

temp$meandist
temp$var
temp$fmean


dist.beta<-rep(0,length(F.lists))
for (j in 1:length(F.lists)){
  dist.beta[j]<-sqrt(sum((F.lists[[j]]-temp$fmean)^2))
}

dist.beta.sign<-dist.beta
for (j in 1:length(dist.beta.sign)){
  if (dist.balanced[j]<dist.caterpilar[j]){
    dist.beta.sign[j]<--dist.beta.sign[j]
  }
}

dist.beta.sort<-sort(dist.beta.sign,index.return=TRUE)
j<-2
plot(seq(1,length(dist.beta)),betaprob[j,dist.kingman.sort$ix],main="Beta b=-0.8",type="S",ylab="Probability",xlab="Trees")

bb<-.5
plot(seq(1,length(dist.beta)),exp(-bb*(dist.beta.sign[dist.kingman.sort$ix])^2)/sum(exp(-bb*sort(dist.beta.sign)^2)),main="Exp-Dist center Beta b=-.4",xlab="Trees",ylab="Probability",type="l")

hist.beta2<-tapply(betaprob[j,],dist.kingman.sign,FUN=sum)
xb<-round(as.numeric(names(hist.beta2)),2)
names(hist.beta2)<-xb
barplot(hist.beta2,main="Beta b=-0.8, n=9",ylab="Probability",xlab="Singed Distance")


j<-29
temp <- fmatrix:::frechet_var_pop(F.lists, betaprob[j, ])

temp$meandist
temp$var
temp$fmean


dist.beta<-rep(0,length(F.lists))
for (j in 1:length(F.lists)){
  dist.beta[j]<-sqrt(sum((F.lists[[j]]-temp$fmean)^2))
}

dist.beta.sign<-dist.beta
for (j in 1:length(dist.beta.sign)){
  if (dist.balanced[j]<dist.caterpilar[j]){
    dist.beta.sign[j]<--dist.beta.sign[j]
  }
}

dist.beta.sort<-sort(dist.beta.sign,index.return=TRUE)
j<-29
plot(seq(1,length(dist.beta)),betaprob[j,dist.kingman.sort$ix],main="Beta b=100",type="S",ylab="Probability",xlab="Trees")

bb<-.5
plot(seq(1,length(dist.beta)),exp(-bb*(dist.beta.sign[dist.kingman.sort$ix])^2)/sum(exp(-bb*sort(dist.beta.sign)^2)),main="Exp-Dist center Beta b=100",xlab="Trees",ylab="Probability",type="l")

hist.beta2<-tapply(betaprob[j,],dist.kingman.sign,FUN=sum)
xb<-round(as.numeric(names(hist.beta2)),2)
names(hist.beta2)<-xb
barplot(hist.beta2,main="Beta b=100, n=9",ylab="Probability",xlab="Singed Distance")

