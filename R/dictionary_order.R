#'  Compare in dictionary order
#'
#'  Compare x and y in dictionary order
#'
#'  @param x
#'  @param y
#'
#'  @export `%<=%`
`%<=%` <- function(x, y){
  n <- min(length(x), length(y))

  for(i in 1:n){
    if(x[i] > y[i]) return(FALSE)
    if(x[i] < y[i]) return(TRUE)
  }

  if(length(y) < length(x)) return(FALSE)
  return(TRUE)
}

#'  Compare in dictionary order in reverse
#'
#'  Compare x and y in dictionary order reverse
#'
#'  @param x
#'  @param y
#'

#' @export `%>=%`
`%>=%` <- function(x, y){
  `%<=%`(y, x)
}



#' Evaluate signed distance of F-matrix
#'
#' Evaluate signed distance of F-matrix
#'
#' @param Fmat
#' @param dist
#' @param precomp     Default is NULL.
#'                    Else, list of Fbal, Funb, Fkin
#'                    which are precomputed
#'
#' @export signed_dist
signed_dist <- function(Fmat, d = c("l2", "l1"),
                        precomp = precomp_bal_unb_kin(ncol(Fmat) + 1)){

  d <- match.arg(d)

  n <- ncol(Fmat) + 1


    Fbal <- precomp$Fbal
    Funb <- precomp$Funb
    Fkin <- precomp$Fkin



  dist <- distance_Fmat(Fmat, Fkin, dist = d)

  #TODO: Equality is not dealt with nicely here.
  sign <- ifelse(distance_Fmat(Fmat, Funb, dist = d) <=
                   distance_Fmat(Fmat, Fbal, dist = d),
                 -1, 1)

  return(sign*dist)
}

#' Precompute matrices for signed distance
#'
#' Precompute matrices for signed distance
#'
#' @param n
precomp_bal_unb_kin <- function(n){
  Fbal <- bal_Fmat(n)
  Funb <- unb_Fmat(n)
  Fkin <- kingman_m(n)

  return(list(Fbal = Fbal,
              Funb = Funb,
              Fkin = Fkin))
}

#' Sort a list of F-matrices
#'
#' Sort a list of F-matrices
#'
#' @param F.list
#'
#'
#' @export sort_Flist
sort_Flist <- function(F.list,
                       n = ncol(F.list[[1]]) + 1,
                       return_dist = FALSE){

  precomp <- precomp_bal_unb_kin(n)

  signed_distances <- sapply(F.list, signed_dist, precomp = precomp)

  o <- order(signed_distances)

  if(return_dist){
    ret <- list(F.list = F.list[o],
                distances = signed_distances)
  } else{
    ret <- F.list[o]
  }

  return(ret)
}



#'  Compare in signed distance order
#'
#'  Compare x and y in signed distance order
#'
#'  @param x
#'  @param y
#'
#'  @export `%<d%`
`%<d%` <- function(x, y, ...){

  f <- function(u){signed_dist(u, ...)}

  if(f(y) < f(x)) return(FALSE)
  return(TRUE)
}

#'  Compare in signed distance order in reverse
#'
#'  Compare x and y in signed distance order in reverse
#'
#'  @param x
#'  @param y
#'

#' @export `%>d%`
`%>d%` <- function(x, y, ...){
  `%<d%`(y, x, ...)
}

