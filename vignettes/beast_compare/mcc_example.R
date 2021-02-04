library(fmatrix)
library(phylodyn)
library("phangorn")
library("mvtnorm")
library("ape")
library(dplyr)

root <- rprojroot::has_file("DESCRIPTION")
setwd(root$find_file("vignettes/beast_compare/"))

dataname <- "temp"
# dataname <- "temp2"
# dataname <- "temp3"
# dataname <- "temp4"

trees_list <- ape::read.nexus(paste0(dataname, ".trees"))

Fmat_list <- lapply(trees_list, gen_Fmat)


mcc_tree2 <- phangorn::mcc(trees_list)
mcc_Fmat2 <- gen_Fmat(mcc_tree2)
mcc_Fmat2 %>% plotF

mcc_tree <- ape::read.nexus(paste0(dataname, ".mcc"))
mcc_Fmat <- gen_Fmat(mcc_tree)
is.Fmat(mcc_Fmat)
plotF(mcc_Fmat)

matrix_F <- matrix_list(Fmat_list)

energy <- function(state, matr=matrix_F, d="l2", type = 0){
  Fmat <- gen_Fmat(tree_from_matching(state))
  n <- ncol(Fmat) + 1
  state.F <- reduce_to_vector(Fmat)

  if(d == "l1") distances <-
    sum(abs( (matr - state.F) %*% diag(1/scaling_vector(n, type = type)) ))/n
  if(d == "l2") distances <-
    sqrt(sum(( (matr - state.F) %*% diag(1/scaling_vector(n, type = type)) )**2)/n)

  return(sum(distances))
}

probab <- function(e1, e2, temp, eps = 0.0001){
  if(e2 < e1) return(1)
  # if(e1 == e2) cat("!")
  return(exp(-(e2-e1+eps)/temp))
}


# schedule <- "lin"
# schedule <- "exp"
schedule <- "log"

dist_type <- 1

n <- 50

# set.seed(4564)
init <- rcoal(n, br = 1)

config <- matching(init)

e1 <- energy(config, type = dist_type) # note sample has been set to default

# plot(init, direction = "downwards")
# tiplabels()
# nodelabels



max.iter <- 25000

init.temp <- 1000
temp <- init.temp
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


start.time <- Sys.time()
for(i in 1:max.iter){
  new <- proposal_matching(config)
  e2 <- energy(new, type = dist_type)

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

  if(schedule == "exp") temp <- alpha*temp
  if(schedule == "log") temp <- init.temp/(1+log(1+i))
  if(schedule == "lin") temp <- init.temp/(1+i)
  # cat(i,"\n")
}

time2 <- Sys.time() - start.time
time2

plot(track.energy.config)

track.config[which.min(track.energy.config), ] %>% tree_from_matching %>% gen_Fmat %>% plotF
mcc_Fmat %>% plotF

energy(matching(mcc_tree2))
min(track.energy.config, na.rm = TRUE)

save(track.config, track.proposal, track.energy.config, track.energy.proposal, track.prob, track.temp, track.acceptance, file = paste(dataname, schedule, "schedule.RData", sep = "_"))



#
# max.iter2 <- 25000
# offset.temp <- nrow(track.config)
#
# alpha <- 0.995
#
# track2 <- T
# if(track2){
#   track2.config <- matrix(NA, nrow = max.iter2, ncol = length(config))
#   track2.proposal <- matrix(NA, nrow = max.iter2, ncol = length(config))
#   track2.energy.config <- rep(NA, max.iter2)
#   track2.energy.proposal <- rep(NA, max.iter2)
#   track2.temp <- rep(NA, max.iter2)
#   track2.acceptance <- rep(NA, max.iter2)
#   track2.prob <- rep(NA, max.iter2)
# }
#
#
# start.time <- Sys.time()
# for(i in 1:max.iter2){
#   new <- proposal_matching(config)
#   e2 <- energy(new)
#
#   prob <- probab(e1, e2, temp)
#
#   if(track2){
#     track2.config[i,] <- config
#     track2.proposal[i,] <- new
#     track2.energy.config[i] <- e1
#     track2.energy.proposal[i] <- e2
#     track2.temp[i] <- temp
#     track2.prob[i] <- prob
#   }
#
#   if(runif(1) <= prob){
#     config <- new
#     e1 <- e2
#     if(track2){
#       track2.acceptance[i] <- 1
#     }
#   } else{
#     if(track2){
#       track2.acceptance[i] <- 0
#     }
#   }
#
#   if(track2 & (i > 500)) if(sum(track2.acceptance[(i-100):(i)])==0) break
#
#   temp <- init.temp/(1+log(1+i+offset.temp))
#   # cat(i,"\n")
# }
#
# track.config <- rbind(track.config, track2.config)
# track.proposal <- rbind(track.proposal, track2.proposal)
# track.energy.proposal <- c(track.energy.proposal, track2.energy.proposal)
# track.energy.config <- c(track.energy.config, track2.energy.config)
# track.acceptance <- c(track.acceptance, track2.acceptance)
# track.temp <- c(track.temp, track2.temp)
# track.prob <- c(track.prob, track2.prob)
