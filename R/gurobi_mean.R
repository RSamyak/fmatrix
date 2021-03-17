
index <- function(i, j){
  if(j > i) stop("i < j")
  (i-1)*(i-1+1)/2 + j
}

inv.index <- function(k, eps = 1e-6){

  i <- ceiling(sqrt(2*k + 0.25) - 0.5 - eps)

  j <- k - ((i-1)*(i-1+1)/2)

  return(c(i, j))
}


#' @export gurobi_kingman_mean
gurobi_kingman_mean <- function(n, qmat = kingman_m(n)){

  # qmat <- kingman_m(n)

  #install.packages('/Library/gurobi911/mac64/R/gurobi_9.1-1_R_4.0.2.tgz', repos=NULL)
  #install.packages('slam')
  #library(gurobi)
  #library(slam)

  # n <- ncol(F.list[[1]]) + 1

  NVar <- n*(n-1)/2


  # temp <- simple_triplet_zero_matrix(nrow = NConst, ncol = NVar)


  myQ <- slam::simple_triplet_diag_matrix(rep(1, NVar))
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

  print(result$objval)
  print(result$x)

  out <- result$x

  outF <- matrix(0, n-1, n-1)
  outF[upper.tri(outF, diag = T)] <- out
  outF <- t(outF)

  return(outF)
}
