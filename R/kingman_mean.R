#' Compute the Kingman M matrix
#'
#' Compute the barycentre of the F-matrices based on a formula
#' under the Kingman model
#'
#' @param n            The number of leaves in the model
#' @param entry        If NULL, full matrix is returned.
#'                     Else, needs to be a 2-tuple (i, j)
#'                     for the desired entry.
#'
#'
#' @export kingman_m
kingman_m <- function(n = NULL, entry = NULL){

  if(is.null(n) & is.null(entry))
    stop("at least one of n and entry is required.")

  if(!is.null(entry)) {
    stopifnot(length(entry) == 2)

    i <- entry[1]
    j <- entry[2]

    if(j > i) return(0)

    return(j*(j+1)/i)
  }

  ## Under the Kingman model,
  ## ret[i, j] = j*(j+1) / i

  columns <- sapply(1:(n-1), function(j) {j*(j+1)})

  ret <- matrix(rep(columns, n-1), n-1, byrow = TRUE)
  ret <- sweep(ret, MARGIN = 1, STATS = 1:(n-1), FUN = "/")

  ret[upper.tri(ret)] <- 0

  return(ret)
}


#' Compute the 'most unbalanced' F-matrix
#'
#' Compute the F-fmatrix of the most unbalanced tree
#'
#' @param n            The number of leaves
#'
#'
#' @export unb_Fmat
unb_Fmat <- function(n = NULL){

  ret <- matrix(0, n-1, n-1)

  Dmat <- ret
  for(j in 1:nrow(Dmat)){
    Dmat[j:(n-1), j] <- c(2, rep(1, n-j-1))
  }

  ret <- Fmat_from_Dmat(Dmat)

  return(ret)
}


#' Compute the 'most balanced' F-matrix
#'
#' Compute the F-fmatrix of the most balanced tree
#'
#' @param n            The number of leaves
#'
#'
#' @export bal_Fmat
bal_Fmat <- function(n = NULL){

  ret <- matrix(0, n-1, n-1)

  for(j in 1:nrow(ret)){
    ret[j:(n-1), j] <- c((j+1):1, rep(0, n))[1:(n-j)]
  }

  return(ret)
}
