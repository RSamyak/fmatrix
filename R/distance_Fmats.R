#' Evaluate Distance between Fmatrices
#'
#' Evaluate Distance between Fmatrices
#'
#' @param Fmat1
#' @param Fmat2
#' @param dist
#'
#' @export

distance_Fmat <- function(Fmat1, Fmat2, dist = c("l2", "l1"), scale = NULL){

  n <- ncol(Fmat1) + 1

  if(is.null(scale)) scale = rep(1, n - 1)
  else if(scale == 1) scale = 2:n

  dist <- match.arg(dist)

  # if(is.Fmat(Fmat1[[1]]))

  ## Vectorise this thing

  if(dist == "l1"){
    return(sum(abs( diag(1/scale) %*% (Fmat1 - Fmat2) ) ))
  }

  if(dist == "l2"){
    return(sqrt(sum(( diag(1/scale) %*% (Fmat1 - Fmat2) )**2)))
  }

  warning("dist not valid, returning error.")
  return(-1)
}


#' Evaluate weighted Distance between Fmatrices
#'
#' Evaluate weighted Distance between Fmatrices
#'
#' @param wFmat1
#' @param wFmat2
#' @param dist
#'
#' @export

distance_Fmat_weighted <- function(wFmat1, wFmat2, dist = c("l2", "l1")){

  dist <- match.arg(dist)

  # if(is.Fmat(Fmat1[[1]]))

  ## Vectorise this thing

  if(dist == "l1"){
    return(sum(abs(wFmat1$f * wFmat$w - wFmat2$f * wFmat2$w)))
  }

  if(dist == "l2"){
    return(sqrt(sum((wFmat1$f * wFmat$w - wFmat2$f * wFmat2$w)**2)))
  }

  warning("dist not valid, returning error.")
  return(-1)
}
