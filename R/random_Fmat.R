#' @param n
#'
#' @export
rFmat <- function(n){

  # Placeholder function for an actual rFmat

  qmat <- matrix(NA, nrow = n-1, ncol = n-1)


  qmat[upper.tri(qmat)] <- 0

  diag(qmat) <- 2:n

  qmat[lower.tri(qmat)] <- rbinom(sum(lower.tri(qmat)), n, prob = .5)


  for(i1 in 1:(n-2)){
    qmat[i1+1, i1] <- i1
  }


  ret <- nearby_Fmat(qmat)

  # if(is.Fmat(ret)) return(ret)
  # else(stop("nearby_Fmat not working as expected.  You should never see this message."))

  return(ret)
}
