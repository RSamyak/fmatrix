

# Simulated annealing to obtain a global optimiser of Sample Frechet mean

############################################
#### Load the dependencies #################
############################################


library(ape)
# devtools::install_github("JuliaPalacios/phylodyn")
library(phylodyn)

# source("Gen_Sim_funs.R")
path <- "./R/"
for(file in list.files(path)) source(file.path(path, file))

############################################


############################################
#### Create sample of F-matrices ###########
############################################
load("data/F-lists.RData")

F.list <- F.list10

size <- 1000 # sample size

n.tr <- length(F.list)
ind <- sample(1:n.tr, size = size, replace = T)
sampleF <- F.list[ind]



x <- F.list[[brute.mean(sampleF = sampleF, totalF.list = F.list, dist = "l2")]]
############################################



############################################
#### Chain Simulation functions ############
############################################

energy <- function(state, sample=sampleF, d="l2"){
  state.F <- gen_Fmat(tree_from_matching(state))
  
  if(d == "l1") distances <- sapply(sampleF, function(x){sum(abs(x - state.F))})
  if(d == "l2") distances <- sapply(sampleF, function(x){sum((x - state.F)**2)})
  
  return(sum(distances))
}

probab <- function(e1, e2, temp, eps = 0.0001){
  if(e2 < e1) return(1)
  if(e1 == e2) cat("!")
  return(exp(-(e2-e1+eps)/temp))
}


############################################
#### Chain parameters ######################
############################################


n <- 300

n.sample <- 30
sampleF <- list()
for(i in 1:n.sample){
  sampleF[[i]] <- rFmat(n)
}


init <- rcoal(n, br = 1)

config <- matching(init)

e1 <- energy(config) # note sample has been set to default

# plot(init, direction = "downwards")
# tiplabels()
# nodelabels



max.iter <- 5000

temp <- 1000
alpha <- 0.995

track <- T
if(track){
  track.config <- matrix(NA, nrow = max.iter, ncol = length(config))
  track.proposal <- matrix(NA, nrow = max.iter, ncol = length(config))
  track.energy.config <- rep(NA, max.iter)
  track.energy.proposal <- rep(NA, max.iter)
  track.temp <- rep(NA, max.iter)
  track.acceptance <- rep(NA, max.iter)
  track.prob <- rep(NA, max.iter)
}



for(i in 1:max.iter){
  new <- proposal_matching(config)
  e2 <- energy(new)
  
  prob <- probab(e1, e2, temp)
  
  if(track){
    track.config[i,] <- config
    track.proposal[i,] <- new
    track.energy.config[i] <- e1
    track.energy.proposal[i] <- e2
    track.temp[i] <- temp
    track.prob[i] <- prob
  }
  
  if(runif(1) <= prob){
    config <- new
    e1 <- e2
    if(track){
      track.acceptance[i] <- 1
    }
  } else{
    if(track){
      track.acceptance[i] <- 0
    }
  }
  
  if(track & (i > 500)) if(sum(track.acceptance[(i-100):i])==0) break
  
  temp <- alpha*temp
  cat(i,"\n")
}

# max.iter <- 100
# 
# set.seed(1)
# init <- rcoal(n, br = 1)
# config <- matching(init)
# ma <- config
# for(i in 1:max.iter){
#   # x <- .Random.seed
#   
#   # save(list=c("x", "ma"), file = "temp.RData")
#   
#   ma <- proposal(ma)
#   
#   cat(ma, "\n")
#   cat(i, energy(ma), "\n")
# }
# 



# proposal(config)

png("png/track.temp.png")
plot(track.temp)
dev.off()


plot(track.acceptance)

plot(track.prob)

plot(track.energy.proposal - track.energy.config)

plotF(x)
plotF(gen_Fmat(tree_from_matching(config)))

energy(config)
energy(matching(mytree_from_F(x)))

library(ggplot2)

ggplot(data = NULL) + 
  geom_line(aes(x = 1:max.iter, y = track.energy.proposal-track.energy.config))

ggplot(data = NULL) + 
  geom_line(aes(x = 1:max.iter, y = track.energy.config))




## Points:
## 
## - Tempering schedule - how to initialise?
## - Multiple initialisations
## - Keep track of minimal energy, last entry may only be local minimum
## - neighbours/proposal seems to work fine?
## 
## . Profiled code, energy / gen_Fmat is most expensive 
## . compare code in terms of memory/time to nearby_Fmats procedure
## 










