library(ape)

rm(list=ls())

source("Gen_Sim_funs.R")

####Example
n <- 8
res<-construct_trees(n)
nrow(res$res)


#Indeed for n=5, there are 5 trees
#we construct the Fmatrices

F.list<-list()
for (j in 1:nrow(res$res)){
  F.list[[j]]<- Fmat_from_construct(n, res$res[j,], res$F.matrixRel)
}


n.tr <- length(F.list)



# Calculate pairwise distances brute force
Dmat.l1 <- matrix(NA, nrow=n.tr, ncol=n.tr)

for (i in 1:n.tr) {
  for(j in 1:i){
    Dmat.l1[i,j] <- Dmat.l1[j,i] <- sum(abs(F.list[[i]] - F.list[[j]]))
  }
}



plot(table(Dmat.l1[Dmat.l1>0])/2,ylab="frequency",main=paste("Histogram of pairwise d1, n=", n))

mean<-apply(Dmat.l1,1,mean)
min.ind <- which.min(mean)

dist_to_mean <- Dmat.l1[,min.ind]

plot(table(dist_to_mean),ylab="frequency",main=paste("Histogram of d1 to the mean, n=",n),xlab="d(T,T_mean)")


# Kingman model prob
kingman.prob <- sapply(F.list, FUN = kingman_likelihood)

##### SET beta HERE TO WHICHEVER ONE YOU WANT
beta.prob <- sapply(F.list, FUN = beta_BF_likelihood, betas = 0)


weights <- beta.prob


##### NOTE: Assumes that only one value was passed to betas parameter in defining beta.prob
wDmat.l1 <- Dmat.l1 %*% diag(weights)

meanw<-apply(wDmat.l1,1,mean)
wmin.ind <- which.min(meanw)

plotF(F.list[[wmin.ind]])
plotF(F.list[[min.ind]])

wdist_to_mean <- wDmat.l1[,wmin.ind]


# plot(density(dist_to_mean),xlab="distance to mean", ylab="Probability",main="Uniform vs Kingman",ylim=c(0,.15))
# points(tableW,type="l",col="red")

######################################################################

mybigF <- matrix(0, 7,7)
diag(mybigF) <- 2:8
mybigF[2:7,1] <- c(1,0,0,0,0,0)
mybigF[3:7,2] <- c(2,1,1,1,1)
mybigF[4:7,3] <- c(3,2,1,1)
mybigF[5:7,4] <- c(4,3,3)
mybigF[6:7,5] <- c(5,4)
mybigF[7,6] <- 6

plotF(mybigF)
num_cherries(mybigF)

tmp <- unlist(lapply(F.list, beta_BF_likelihood, betas = -1))

tmp2 <- unlist(lapply(F.list, kingman_likelihood))




######################################################################
## Testing 




test <- meanF(F.list, weights = beta.prob)

# test <- meanF(F.list[sample(1:length(F.list), size=rpois(1, lambda=20))])

is.Fmat(round(test))

test2 <- closest_Fmat(test, F.list)
test3 <- closest_Fmat(round(test), F.list)

identical(test2, F.list[[wmin.ind]])
identical(test2, test3)

test2 - test3
test2 - F.list[[wmin.ind]]
test3 - F.list[[wmin.ind]]

plotF(test2)
plotF(test3)
plotF(F.list[[wmin.ind]])

where.Fmat.debug(test, tol = 1e-6)

test4 <- nearby_Fmat(test)
test5 <- nearby_Fmat(round(test))

plotF(test4)
plotF(test5)



################################################ SIMULATIONS

l <- 6

len1 <- len2 <- c()

out <- matrix(NA, nrow=0, ncol=l+3)


n.sim <- 100
size <- 50

counts <- rep(0, length = l)

set.seed(123)
for(i in 1:n.sim){
  ind <- sample(1:n.tr, size = size, replace = T, prob = beta.prob)
  sampleF <- F.list[ind]
  
  target <- F.list[[brute.mean(sampleF, F.list)]]
  
  test <- meanF(sampleF)
  
  test1 <- closest_Fmat(test, F.list)
  test2 <- closest_Fmat(round(test), F.list)
  test3 <- nearby_Fmat(test)
  test4 <- nearby_Fmat(round(test))
  lis1 <- nearby_Fmats(test)
  lis2 <- nearby_Fmats(round(test))
  test5 <- closest_Fmat(test, lis1)
  test6 <- closest_Fmat(test, lis2)
  
  len1 <- c(len1, length(lis1))
  len2 <- c(len2, length(lis2))
  
  counts <- counts + c(identical(target, test1),
                       identical(target, test2),
                       identical(target, test3),
                       identical(target, test4),
                       identical(target, test5),
                       identical(target, test6))
  cat(i, fill=T)
}

out <- rbind(out, c(n.sim, size, NA, counts/n.sim))
out


################################ Constructing F.lists and means.list, etc
# 
# means2.list <- list()
# 
# for(ind in 11:5){
#   n <- ind
#   res<-construct_trees(n)
#   
#   
#   #Indeed for n=5, there are 5 trees
#   #we construct the Fmatrices
#   
#   F.list<-list()
#   for (j in 1:nrow(res$res)){
#     F.list[[j]]<- Fmat_from_construct(n, res$res[j,], res$F.matrixRel)
#   }
#   
#   
#   n.tr <- length(F.list)
#   
#   min.ind <- 1
#   
#   Fmean <- F.list[[min.ind]]
#   
#   min.dist <- mean(sapply(F.list, function(x) sum(abs(x - Fmean))))
#   
#   for (i in 1:n.tr){
#     dist <- mean(sapply(F.list, function(x){sum(abs(x - F.list[[i]]))}))
#     
#     if(dist<= min.dist){cat(i, " ")}
#     if(dist == min.dist){cat("!")}
#     
#     if(dist < min.dist){
#       min.ind <- i
#       Fmean <- F.list[[min.ind]]
#       min.dist <- dist
#     }
#   }
#   
#   # Calculate pairwise distances brute force
#   # Dmat.l1 <- matrix(NA, nrow=n.tr, ncol=n.tr)
#   # 
#   # for (i in 1:n.tr) {
#   #   for(j in 1:i){
#   #     Dmat.l1[i,j] <- Dmat.l1[j,i] <- sum(abs(F.list[[i]] - F.list[[j]]))
#   #   }
#   # }
#   # 
#   # 
#   # mean<-apply(Dmat.l1,1,mean)
#   # min.ind <- which.min(mean)
#   
#   means2.list <- append(means2.list, list(F.list[[min.ind]]))
# }
# means.list
# 
# 
# list.F.lists <- list()
# for(ind in 5:11){
#   n <- ind
#   res<-construct_trees(n)
#   
#   
#   #Indeed for n=5, there are 5 trees
#   #we construct the Fmatrices
#   
#   F.list<-list()
#   for (j in 1:nrow(res$res)){
#     F.list[[j]]<- Fmat_from_construct(n, res$res[j,], res$F.matrixRel)
#   }
#   
#   list.F.lists <- append(list.F.lists, list(F.list))
# }
# 
# F.list5 <- list.F.lists[[1]]
# F.list6 <- list.F.lists[[2]]
# F.list7 <- list.F.lists[[3]]
# F.list8 <- list.F.lists[[4]]
# F.list9 <- list.F.lists[[5]]
# F.list10 <- list.F.lists[[6]]
# F.list11 <- list.F.lists[[7]]
# 
# D.list5 <- lapply(F.list5, Dmat_from_Fmat)
# D.list6 <- lapply(F.list6, Dmat_from_Fmat)
# D.list7 <- lapply(F.list7, Dmat_from_Fmat)
# D.list8 <- lapply(F.list8, Dmat_from_Fmat)
# D.list9 <- lapply(F.list9, Dmat_from_Fmat)
# D.list10 <- lapply(F.list10, Dmat_from_Fmat)
# D.list11 <- lapply(F.list11, Dmat_from_Fmat)
# 
# 
# 
# identical(F.list10[[1000]], Fmat_from_Dmat(D.list10[[1000]]))
# 
# Fmat_from_Dmat(D.list10[[1000]])

########################### For constructing distmatrices.RData
# for(ind in 5:9){
#   D.list <- get(paste("D.list", ind, sep=""))
#   n.tr <- length(D.list)
#   
#   D.l2 <- D.l1 <- matrix(NA, n.tr, n.tr)
#   
#   for(i in 1:n.tr){
#     for(j in 1:i){
#       D.l1[i,j] <- D.l1[j,i] <- sum(abs(D.list[[i]] - D.list[[j]]))
#       D.l2[i,j] <- D.l2[j,i] <- sqrt(sum((D.list[[i]] - D.list[[j]])**2))
#     }
#   }
#   
#   assign(paste("Dmat-D.l1.",ind,sep=""), D.l1)
#   assign(paste("Dmat-D.l2.",ind,sep=""), D.l2)
# }
# 
# save(list = c(paste("Dmat-D.l1.",5:9,sep=""),paste("Dmat-D.l2.",5:9,sep="")), file = "distmatrices-D.RData")


########################### 
o <- which(apply(Dmat.l2.8,1,function(x){sum(x <= 1) > 1}) == 0)

o <- which(Dmat.l1.8[,71] <= 2)

plotF.list(F.list8[o])


############################

n.sim <- 1000



F.list <- F.list9
n.tr <- length(F.list)

beta.prob <- sapply(F.list, FUN = beta_BF_likelihood, betas = 10^5)
out <- c()
for(i in 1:n.sim){
  f1 <- F.list[sample(1:n.tr, 2, prob = beta.prob, replace = T)]
  new <- sum(abs(f1[[1]] - f1[[2]]))
  out <- append(out, new)
}
hist(out)
