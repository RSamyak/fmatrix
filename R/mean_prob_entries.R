

kingman_meanF_entry <- function(i, j){
  ## This function returns the mean F-matrix entry value
  ## under the Tajima / Kingman / Yule model
  ## at row i and column j
  stopifnot(j <= i)
  j * (j+1) / i
}

kingman_meanD_entry <- function(i, j){
  ## This function returns the mean D-matrix entry value
  ## under the Tajima / Kingman / Yule model
  ## at row i and column j
  stopifnot(j <= i)
  2 * j / i
}

kingman_probD_entry <- function(i, j, dd = 2){
  ## This function returns the probability
  ## under the Tajima / Kingman / Yule model
  ## of D[i, j] being equal to dd
  stopifnot(j <= i)

  if(dd == 2){
    return(j * (j-1) / (i * (i-1)))
  }

  if(dd == 1){
    return(2 * j * (i-j) / (i * (i-1)))
  }


  if(dd == 0){
    return((i - j) * (i-j-1)/ (i * (i-1)))
  }

  return(0)
}


beta_probD_entry <- function(i, j, bet = 0, dd = 2){
  ## This function returns the probability
  ## under the BF Beta model
  ## of D[i, j] being equal to dd
  stopifnot(j <= i)

  if(dd == 2){
    return(0)
  }

  if(dd == 1){

    ret <- 2*beta(bet + i, bet + 1)/beta(bet + 1, bet + 1)
    return(ret)
  }


  if(dd == 0){
    return(0)
  }

  return(0)

}

