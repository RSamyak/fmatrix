library(phylodyn)
library("phangorn")
library("mvtnorm")
library("ape")
library(fmatrix)


# We then set seed and sample size (n), mutation rate (mu). Samp_times and n below mirrors Exp n=70 in Cappello et al. (2020) [see table 1]

set.seed(154)
n<-50
samp_times = 0
mu=22 #mutation rate



treeP<-rcoal(n)
plot(ladderize(treeP), direction = "downwards")
plotF(gen_Fmat(treeP))

data1<-simulate_data(mu,treeP)
saveRDS(treeP, file = "temp.RDS")
phylodyn:::beastfile(data1, "temp.fasta")


#####

set.seed(1540)

n<-50
samp_times = 0
mu=15 #mutation rate


treeP2<-rcoal(n)
plot(ladderize(treeP2), direction = "downwards")
plotF(gen_Fmat(treeP2))

saveRDS(treeP2, file = "temp2.RDS")
data2<-simulate_data(mu,treeP2)
phylodyn:::beastfile(data2, "temp2.fasta")


#####

set.seed(8964)

n<-50
samp_times = 0
mu=15 #mutation rate


treeP3<-rcoal(n)
plot(ladderize(treeP3), direction = "downwards")
plotF(gen_Fmat(treeP3))

saveRDS(treeP3, file = "temp3.RDS")
data3<-simulate_data(mu,treeP3)
phylodyn:::beastfile(data3, "temp3.fasta")


#####

set.seed(1)

n<-50
samp_times = 0
mu=300 #mutation rate


treeP4<-rcoal(n)
plot(ladderize(treeP4), direction = "downwards")
plotF(gen_Fmat(treeP4))

saveRDS(treeP4, file = "temp4.RDS")
data4<-simulate_data(mu,treeP4)
phylodyn:::beastfile(data4, "temp4.fasta")

