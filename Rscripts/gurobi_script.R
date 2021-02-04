library(fmatrix)
library(phangorn)
library(dplyr)
library(gurobi)
library(slam)

index <- function(i, j){
  if(j > i) stop("i < j")
  (i-1)*(i-1+1)/2 + j
}

inv.index <- function(k, eps = 1e-6){

  i <- ceiling(sqrt(2*k + 0.25) - 0.5 - eps)

  j <- k - ((i-1)*(i-1+1)/2)

  return(c(i, j))
}


do_everything_gurobi <- function(nn, B = 1000){

  F.list <- list()
  for(i in 1:B){
    F.list <- append(F.list, list(phylodyn:::gen_Fmat(rcoal(nn, br = 1))))
  }
    # F.list <- F.list9

  qmat <- meanF(F.list)


  n <- ncol(F.list[[1]]) + 1

  NVar <- n*(n-1)/2


  # temp <- simple_triplet_zero_matrix(nrow = NConst, ncol = NVar)


  myQ <- simple_triplet_diag_matrix(rep(1, NVar))
  myobj <- as.vector(t(qmat)[upper.tri(t(qmat), diag = TRUE)])



  ## diagonal
  myA1 <- slam::simple_triplet_zero_matrix(ncol = NVar, nrow = n-1)
  myrhs1 <- rep(NA, n-1)
  mysense1 <- rep('=', n-1)
  j <- 1
  for(i in 1:(n-1)){
    myA1[j, index(i,i)] <- 1
    myrhs1[j] <- i+1
    j <- j+1
  }

  ## subdiagonal
  myA2 <- slam::simple_triplet_zero_matrix(ncol = NVar, nrow = n-2)
  myrhs2 <- rep(NA, n-2)
  mysense2 <- rep('=', n-2)
  j <- 1
  for(i in 1:(n-2)){
    myA2[j, index(i+1,i)] <- 1
    myrhs2[j] <- i
    j <- j+1
  }

  ## first col condn
  myA3 <- slam::simple_triplet_zero_matrix(ncol = NVar, nrow = 2*(n-1 - 3 + 1))
  myrhs3 <- rep(NA, 2*(n-1 - 3 + 1))
  mysense3 <- rep('<', 2*(n-1 - 3 + 1))
  j <- 1
  for(i in 3:(n-1)){
    myA3[j, index(i,1)] <- 1
    myA3[j, index(i-1, 1)] <- -1
    myrhs3[j] <- 0
    j <- j+1

    myA3[j, index(i,1)] <- -1
    myA3[j, index(i-1, 1)] <- 1
    myrhs3[j] <- 1
    j <- j+1
  }


  ## other col condn
  myA4 <- slam::simple_triplet_zero_matrix(ncol = NVar, nrow = 5*(n-4)*(n-3)/2)
  myrhs4 <- rep(NA, 5*(n-4)*(n-3)/2)
  mysense4 <- rep('<', 5*(n-4)*(n-3)/2)
  j <- 1
  for(i in 4:(n-1)){
    for(k in 2:(i-2)){
      myA4[j, index(i-1,k)] <- 1
      myA4[j, index(i, k)] <- -1
      myrhs4[j] <- 1
      j <- j+1

      myA4[j, index(i, k-1)] <- 1
      myA4[j, index(i, k)] <- -1
      myrhs4[j] <- 0
      j <- j+1

      myA4[j, index(i, k-1)] <- 1
      myA4[j, index(i-1, k)] <- 1
      myA4[j, index(i-1, k-1)] <- -1
      myA4[j, index(i, k)] <- -1
      myrhs4[j] <- 1
      j <- j+1

      myA4[j, index(i, k)] <- 1
      myA4[j, index(i-1, k)] <- -1
      myrhs4[j] <- 0
      j <- j+1

      myA4[j, index(i, k)] <- 1
      myA4[j, index(i, k-1)] <- -1
      myA4[j, index(i-1, k)] <- -1
      myA4[j, index(i-1, k-1)] <- 1
      myrhs4[j] <- 0
      j <- j+1
    }
  }

  myA <- rbind(myA1, myA2, myA3, myA4)
  myrhs <- c(myrhs1, myrhs2, myrhs3, myrhs4)
  mysense <- c(mysense1, mysense2, mysense3, mysense4)
  model <- list()

  model$modelsense <- 'min'
  model$vtype <- rep('I', NVar)

  model$obj     <- -2*myobj
  model$Q     <- myQ

  model$A <- myA
  model$rhs   <- myrhs

  model$sense <- mysense

  result <- gurobi(model)

 return(result)
}


nn_range <- 5:21

time <- rep(NA, length(nn_range))
result_list <- list()

for(i in 1:length(nn_range)){
  start_time <- Sys.time()
  result_list[[i]] <- do_everything_gurobi(nn_range[i])
  time[i] <- as.difftime(Sys.time() - start_time, units = "secs")
}

time[16:17] <- time[16:17]*60

library(tidyverse)
ggplot(data = NULL) +
  geom_point(aes(x = nn_range, y = time)) +
  geom_line(aes(x = nn_range, y = time)) +
  scale_y_log10() +
  theme_bw() +
  labs(x = "n", y = "Computing time (seconds)")
