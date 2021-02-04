
n <- 15
N <- 100


library(fmatrix)
gen_Fmat <- phylodyn:::gen_Fmat
library(tidyverse)


F.list <- list()
for(ii in 1:N){
  F.list[[ii]] <- phylodyn:::gen_Fmat(rcoal(n, br = 1))
}



matrix_F <- matrix_list(F.list)

energy <- function(state, matr=matrix_F, d="l2"){

  state.F <- reduce_to_vector(state)

  if(d == "l1") distances <- sum(abs(matr - state.F))/n
  if(d == "l2") distances <- sqrt(sum((matr - state.F)**2)/n)

  return(sum(distances))
}


testF2 <- fmatrix:::sa_mean(F.list)
testF2 %>% plotF



Fm <- meanF(F.list)

Dm <- fmatrix:::Dmat_from_Fmat(Fm)

# DmRd <- round(Dm) ### Do special rounding here

Dsm <- apply(Dm, 2, diff)
for(i in 1:nrow(Dsm)){
  Dsm[i, i+1] <- 0
}

apply(Dsm, 2, sum)

ret <- matrix(0, nrow(Dm), ncol(Dm))
ret[1, ] <- Dm[1, ]

Dst <- rep(NA, nrow(Dm))

for(j in 2:nrow(Dm)){
  tmpret <- ret[j-1,]
  tmpret[j] <- 2

  Dst[j] <- which.max(abs(Dsm[j-1,]))

  tmpret[Dst[j]] <- tmpret[Dst[j]] - 1

  ret[j, ] <- tmpret
}

testF <- fmatrix:::Fmat_from_Dmat(ret)

plotF(testF)


sqrt(sum((testF - testF2)**2))

energy(testF)
energy(testF2)


