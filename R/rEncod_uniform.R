



switchSignsForAlternatingPerm <- function(X, raw = FALSE){
  flags <- rep(c(FALSE, TRUE), length.out = length(X))

  Y <- X
  Y[flags] <- 1 - X[flags]
  if(raw){
    return(Y)
  } else{
    return(rank(Y))
  }
}


rUnifAlternatingPerm <- function(n, max.iter = 1000, raw = FALSE){

  ### Generates random uniform alternating
  ### permutation according to Marchal (2012)

  mysin <- function(x){sin(x*pi/2)}

  counter <- 0

  while(TRUE){

    counter <- counter + 1
    if(counter > max.iter) {return(NA)}

    U <- runif(n)

    X <- rep(NA, n)

    X[1] <- U[1]

    for(i in 1:(n-1)){
      X[i+1] = 1 - (2/pi)*asin(U[i+1] * mysin(X[i]))
    }

    alpha <- mysin(X[n])/mysin(X[1])

    thresh <- 1/(alpha + (1/alpha))


    t <- runif(1)

    if(t < thresh){
      return(switchSignsForAlternatingPerm(X, raw = raw))
    }
    else if(t < 2*thresh){
      return(switchSignsForAlternatingPerm(rev(X), raw = raw))
    }
  }

}


relativeComplement <- function(perm){
  ordered <- sort(perm)
  rev(ordered)[rank(perm)]
}


nestledTransform <- function(perm){
  if(length(perm) > 1)
    if (which.min(perm) < which.max(perm)){
    perm <- relativeComplement(perm)
  }
  perm
}


nestling <- function(perm){
  ## Recursive function to do the nestling of a permutation
  perm <- nestledTransform(perm)

  if(length(perm) == 0) {return(NA)}

  index <- which.min(perm)

  first_part <- nestledTransform(perm[0:(index-1)])

  suffix <- 0
  if(index + 1 <= length(perm)){suffix <- (index+1):length(perm)}
  second_part <- nestledTransform(perm[suffix])

  return(list(nestling(first_part),
              min(perm),
              nestling(second_part)))
}


tabulate_nestling <- function(nestled, tab = matrix(nrow = 0, ncol=2), parent = 0){
  ## Recursive function to tabulate the nodes and their parents

  if(identical(NA, nestled)){return(tab)}

  return(rbind(tabulate_nestling(nestled[[1]], tab = tab, parent = nestled[[2]]),
               matrix(c(nestled[[2]], parent), ncol = 2),
               tabulate_nestling(nestled[[3]], parent = nestled[[2]])))

}

#' Generate a random encoding
#'
#' Generate a random encoding uniformly
#'
#' @param n      n+1 is the number of tips of the tree
#' @param type   flag to switch between isochronous and hetero
#'                 chronous versions of encodings
#'
#' @export rEncod.single.uniform
rEncod.single.uniform <- function(n,
                           type = c("encod", "hEncod")){

  ### Uses bijection to Alternating permutations
  ### from Donaghey (1974)

  type <- match.arg(type)

  perm <- rUnifAlternatingPerm(n-1)

  tab <- tabulate_nestling(nestling(perm))

  ## We add +1 to match the notation for labelling in encod
  ##
  encod <- tab[order(tab[,1]), 2] + 1

  ret <- NA
  if(type == "encod"){
    ret <- encod[-1]
  }
  if(type == "hEncod"){
    ret <- fmatrix:::preprocess_hEncod(encod)
  }

  return(ret)
}
