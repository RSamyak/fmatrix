#' @param upperrow
#' @param lowerrow
find_drop <- function(upperrow, lowerrow){
  check <- lowerrow == (upperrow - 1)
  out <- which(check)[1]
  return(out + 1)
}


#' @param Fmat
#'
#' @export my_encod
my_encod <- function(Fmat){
  n <- ncol(Fmat) + 1
  nn <- n - 1

  out <- c()

  for(ind1 in 1:(nn-1)){
    out <- append(out, find_drop(Fmat[ind1, ], Fmat[ind1 + 1, ]))
  }
  return(out)
}


#' @param myencod
#'
#' @export Fmat_from_myencod
Fmat_from_myencod <- function(myencod){

  if(myencod[1] == 1){myencod <- myencod[-1]}

  n <- length(myencod) + 2

  ret <- matrix(0, n-1, n-1)

  diag(ret) <- 2:n

  for(i in 1:(n-1)){
    if (i == 1) {
      next
    }

    if (i == 2) {
      ret[i, 1:(i-1)] <- 1; next
    }

    dropped <- myencod[i-1] - 1

    ret[i, 1:(i-1)] <- c(ret[i-1, 0:(dropped - 1)], ret[i-1, (dropped:(i-1))] - 1)


  }
  return(ret)
}


#' @param myencod
#'
#' @export is.myencod
is.myencod <- function(string, digits.tol = 6, debug = FALSE){

  string <- round(string, digits = digits.tol)

  error_ndigits <- floor(log10(length(string))) + 1

  ret <- 1

  for(i in 1:length(string)){
    if(! string[i] == round(string[i]) ) { ret <- (-2 - (i / 10**error_ndigits)); break }
    if(! string[i] <= (i+1) |
       ! string[i] >= 2) { ret <- (-3 - (i / 10**error_ndigits)); break }
  }
  if(! all(table(string) <= 2)) { ret <- (-4) }

  if(debug) return(ret)

  return(ret > 0)
}

sample.vec <- function(x, ...) x[sample(length(x), ...)]

#' @param myencod
#' @param check          Do we check if the input myencod is valid?
#'                           Default is FALSE
#'
#' @export proposal_myencod
proposal_myencod <- function(myencod, check = FALSE, digits.tol = 6){

  myencod <- round(myencod, digits = digits.tol)

  if(check) stopifnot(is.myencod(myencod))

  n <- length(myencod) + 2

  if(n <= 3) return(myencod)

  ret <- myencod

  i <- sample.vec(seq(2, n-2),1) # choose an index to change

  eligible <- 2:(i+1)

  while(length(eligible) > 0){
    proposed <- sample.vec(eligible, 1)

    if(sum( myencod[-i] == proposed ) >= 2){
      eligible <- setdiff(eligible, proposed)
      next
    }
    else{
      ret[i] <- proposed

      if(length(ret) != length(myencod)) browser()

      # if(check) stopifnot(is.myencod(ret))
      if(check) {if(!is.myencod(ret)) browser()}
      return(ret)
    }
  }

  stop("Invalid status, eligible entries exhausted without finding valid proposal.")
  return(-1)
}


