#' @param Dmat
#'
#' @param export Fmat_from_Dmat
Fmat_from_Dmat <- function(Dmat){
  n <- ncol(Dmat) + 1

  Fmat <- matrix(0, n-1, n-1)
  Fmat[, 1] <- Dmat[, 1]

  if(n>=2) {
    for(i in 2:(n-1)){
      if(i == n-1) {
        Fmat[n-1, n-1] <- sum(Dmat[n-1, 1:(n-1)])
        break
      }
      Fmat[i:(n-1), i] <- apply(Dmat[i:(n-1), 1:i], 1, sum)
    }

  }
  return(Fmat)
}

#' @param Fmat
#'
#' @param export Dmat_from_Fmat
Dmat_from_Fmat <- function(Fmat){
  n <- ncol(Fmat) + 1

  Dmat <- matrix(0, n-1, n-1)
  Dmat[, 1] <- Fmat[, 1]

  if(n>=2) for(i in 2:(n-1)){
    Dmat[i:(n-1), i] <- Fmat[i:(n-1), i] - Fmat[i:(n-1), i-1]
  }
  return(Dmat)
}
