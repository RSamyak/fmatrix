#### Temporary file to check Beta model Fmat entry probs from mean_prob_entries.R


kingprob <- sapply(F.list9, kingman_likelihood)
mean_kingman <- meanF(F.list9, weights = kingprob)

mean_kingman

meanD <- fmatrix:::Dmat_from_Fmat(mean_kingman)

betaprob <- sapply(F.list9, beta_BF_likelihood, betas = pi)
meanFb <- meanF(F.list9, weights = betaprob)
meanDb <- fmatrix:::Dmat_from_Fmat(meanFb)

betaprob8 <- sapply(F.list8, beta_BF_likelihood, betas = pi)
meanFb8 <- meanF(F.list8, weights = betaprob8)
meanDb8 <- fmatrix:::Dmat_from_Fmat(meanFb8)

betaprob7 <- sapply(F.list7, beta_BF_likelihood, betas = pi)
meanFb7 <- meanF(F.list7, weights = betaprob7)
meanDb7 <- fmatrix:::Dmat_from_Fmat(meanFb7)


bet <- pi


probcheck <- function(kk, nn, dd = 2, weights = kingprob){
  out <- sapply(D.list9, function(u){u[nn, kk] == dd})
  sum(out * weights)
}



ret1 <- ret2 <- ret0 <- matrix(NA, 8,8)
for(i in 1:8){
  for(j in 1:8){
    ret0[i,j] <- probcheck(j,i,0, betaprob)
    ret1[i,j] <- probcheck(j,i,1, betaprob)
    ret2[i,j] <- probcheck(j,i,2, betaprob)

  }
}
#
# ret1 * ((1:8 * c(1,(1:7))) %*% t(rep(1,8)))/2
#
# tmp <- matrix(0,8,8)
# for(n in 1:8){
#   for(j in 1:8){
#     tmp[n,j] <- j*(n-j)
#   }
# }
#
# tmp[upper.tri(tmp)] <- 0
# tmp

i <- 3; j <- 2; ret1[i, j]; 2*(i-2)*beta(pi+i-2,pi+3)/beta(pi+1,pi+1)

i <- 4; j <- 2; tmp <- ret1[i,j] * beta(pi+1, pi + 1)

mycheck <- function(tmp, tol = 1e-6, ilimit = c(-5,5), jlimit = c(-5,5)){

  for(ii in ilimit[1]:ilimit[2]){
    for(jj in jlimit[1]:jlimit[2]){
      suppressWarnings({
        guess <- beta(pi + ii, pi + jj)
      })
      if(is.na(guess)) next
      test <- tmp / guess
      if(abs(test - round(test)) < tol) return(c(ii, jj))
    }
  }
  return(-1)
}

mycheck(ret1[4,2] * beta(pi+1, pi+1), tol = 1e-6, ilimit = c(-2,6), jlimit = c(-2,6))


myfun <- function(n, bet = 0){

  sum <- 0

  for(k in 1:max(1,n-2)){
    if(n-k-2 >= 0){
      for(j in 0:(n-k-2)){
        add <- choose(n-k-2, j) * (beta(bet+j+3, bet+n-j-2) * beta(bet+j+2,bet+1)) / (beta(bet + 1, bet + 1)**2)
        if(is.na(add)) browser()
        sum <- sum + add
      }
    }
  }


  sum <- 4*sum
  return(sum)
}

myfun(6, bet=pi)
ret1
