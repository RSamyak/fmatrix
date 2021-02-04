#' @param Fmat
#' @param tol
#'
#' @export is.Fmat.debug
is.Fmat.debug <- function(Fmat, tol=0){
  ## Function to test if a given matrix is a valid F-matrix, returns a negative value if not

  n <- ncol(Fmat) + 1
  if(! nrow(Fmat) == (n-1)) return(-0.1)
  if(! all( abs(Fmat[upper.tri(Fmat)] - 0) <= tol)) return(-2)
  if(! all( abs(diag(Fmat) - 2:n) <= tol)) return(-3)

  for(i1 in 1:(n-2)){
    if(!abs(Fmat[i1+1,i1] - i1) <= tol) return(-4)
  }

  if(n >=4 ){
    for(i1 in 3:(n-1)){
      if(!all(any(abs(Fmat[i1, 1] - Fmat[i1-1,1]) <= tol,
                  abs(Fmat[i1, 1] - Fmat[i1-1,1] + 1) <= tol
      ),
      Fmat[i1, 1] >= -tol
      )) return(-5)

      if(i1 >=4){
        for(i2 in 2:(i1 - 2)){
          if(!all(any(abs(Fmat[i1, i2] - Fmat[i1-1,i2]) <= tol,
                      abs(Fmat[i1, i2] - Fmat[i1-1,i2] + 1) <= tol),
                  Fmat[i1, i2] >= Fmat[i1, i2-1] - tol,
                  Fmat[i1-1, i2] - Fmat[i1, i2] >= Fmat[i1-1, i2-1] - Fmat[i1, i2-1] - tol
          )) return(-6)
        }
      }
    }
  }

  return(T)
}

#' @param Fmat
#' @param tol
#'
#' @export is.Fmat
is.Fmat <- function(Fmat){
  return(is.Fmat.debug(Fmat)>0)
}



#' @param Fmat
#' @param tol
#'
#' @export where.Fmat.debug
where.Fmat.debug <- function(Fmat, tol=0){
  n <- ncol(Fmat) + 1
  if(! nrow(Fmat) == (n-1)) return(-0.1)
  if(! all( abs(Fmat[upper.tri(Fmat)] - 0) <= tol)) return(-2)
  if(! all( abs(diag(Fmat) - 2:n) <= tol)) return(-3)

  for(i1 in 1:(n-2)){
    if(!abs(Fmat[i1+1,i1] - i1) <= tol) return(-4)
  }

  if(n >= 4){

    flag <- matrix(NA, n-1, n-1)

    for(i1 in 3:(n-1)){
      if(!all(any(abs(Fmat[i1, 1] - Fmat[i1-1,1]) <= tol,
                  abs(Fmat[i1, 1] - Fmat[i1-1,1] + 1) <= tol
      ),
      Fmat[i1, 1] >= -tol
      )) {
        flag[i1, 1] <- F
      } else flag[i1, 1] <- T

      if(i1 >=4){
        for(i2 in 2:(i1 - 2)){
          if(!all(any(abs(Fmat[i1, i2] - Fmat[i1-1,i2]) <= tol,
                      abs(Fmat[i1, i2] - Fmat[i1-1,i2] + 1) <= tol
          ),
          Fmat[i1, i2] >= Fmat[i1, i2-1] - tol,
          Fmat[i1-1, i2] - Fmat[i1, i2] >= Fmat[i1-1, i2-1] - Fmat[i1, i2-1] - tol
          )) {
            flag[i1, i2] <- F
          } else flag[i1, i2] <- T
        }
      }
    }
  }

  return(flag)
}
