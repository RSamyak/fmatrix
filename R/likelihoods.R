#' @param Fmat
#'
#' @export kingman_likelihood

kingman_likelihood <- function(Fmat){
  ## Simple function to calculate probability of F-matrix under Kingman Model

  n <- ncol(Fmat) + 1
  c <- num_cherries(Fmat)
  fact <- (2^(n-1))/factorial(n-1)
  return(fact*2^(-c))
}


## Functions to calculate likelihood of F-matrix under Blum-Francois' Beta-splitting model (only beta_BF_likelihood to be exported, rest are internal)

#' @param num
#' @param denom
beta_BF_expression <- function(num, denom){
  k <- length(num)
  if(!length(denom)==k) stop("same length needed")

  t <- function(x){
    out <- 1
    for(ind in 1:k){
      out <- out * (x+num[ind])/(2*x + denom[ind])
    }
    return(out)
  }
  return(t)
}




#' @param myencod
#' @param node
child_seq <- function(myencod, node){
  nn <- length(myencod) + 1

  left <- 2

  childL <- c()
  childR <- c()

  for (ind1 in (node-1):(nn-1)){

    if(myencod[ind1] == node){
      if(left == 2) childR <- append(childR, ind1 + 2)
      if(left == 1) childL <- append(childL, ind1 + 2)

      left <- left - 1
    }
    else if(myencod[ind1] %in% childL){
      childL <- append(childL, ind1 + 2)
    } else if(myencod[ind1] %in% childR){
      childR <- append(childR, ind1 + 2)
    }
  }
  return(list(childL=childL
              , childR=childR
  ))
}

#' @param myencod
#' @param beta
beta_BF_numden_from_encod <- function(myencod, beta = 0){

  n <- length(myencod) + 2

  num <- c()
  den <- c()

  for(ind1 in 4:n){ # ind1 is size of tree

    for(ind2 in (ind1-1):2 ){ # ind2 looks for ancestry

      childs <- child_seq(myencod[1:(ind1-2)], ind2)
      childL <- childs$childL
      childR <- childs$childR

      if(ind1 %in% childL){
        l <- length(childL)
        r <- length(childR) + 1
        if(!(l==1 & r == 1)){
          num <- append(num, l)
          den <- append(den, l+r)
        }
      } else if(ind1 %in% childR){
        l <- length(childL) + 1
        r <- length(childR)
        if(!(l==1 & r == 1)){
          num <- append(num, r)
          den <- append(den, l+r)
        }
      }

    }

  }
  return(list(num=num,den=den))
}

#' @param Fmat
#' @param betas
#'
#' @export beta_BF_likelihood
beta_BF_likelihood <- function(Fmat, betas=0){
  #### Note this is Blum and Francois (2006) Beta-splitting model, not Aldous; as described in Sainudiin, Veber (2016)
  #### Only export this function, the rest are called

  n <- ncol(Fmat) + 1

  if(n <= 3) return(rep(1,length(betas)))

  numden <- beta_BF_numden_from_encod(my_encod(Fmat))
  num <- numden$num
  den <- numden$den

  func <- beta_BF_expression(num, den)
  return(sapply(betas, FUN = func))
}



# # Aldous Beta-splitting model - no explicit likelihood, what to do?
# beta_aldous_likelihood <- function(Fmat, betas=0){
#   n <- ncol(Fmat) + 1
#
#   if(n <= 3) return(rep(1,length(betas)))
#
#
# }
