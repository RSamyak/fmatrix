
#' @param qmat
#' @param rd
#' @param tol
#' @param ret.flag
#'
#' @export
nearby_Fmat <- function(qmat, rd = .25, tol = 1e-5, ret.flag = F){
  ## Candidate Fmat
  if(is.Fmat(qmat)) return(qmat)
  n <- ncol(qmat) + 1

  if(! nrow(qmat) == (n-1)) stop("Number of rows and columns do not match")

  Fmat <- matrix(NA, n-1, n-1)
  flag <- matrix(NA, n-1, n-1)

  if(! all(abs(qmat[upper.tri(qmat)]-0) < tol)) stop("Upper triangular part of qmat is not zero")

  Fmat[upper.tri(Fmat)] <- 0

  if(! all(abs(diag(qmat) - 2:n) < tol) ) stop("Diagonal of qmat improper")

  diag(Fmat) <- 2:n

  for(i1 in 1:(n-2)){
    if(abs(qmat[i1+1,i1] - i1) > tol) stop("Subdiagonal of qmat improper")
    Fmat[i1+1, i1] <- i1
  }

  if(n >= 4){

    for(i1 in 3:(n-1)){

      current_cell_res <- nearby_Fmat_check_neighbour_col1(current_cell = qmat[i1, 1],
                                                           up = Fmat[i1-1, 1],
                                                           rd = rd,
                                                           tol = tol)
      Fmat[i1, 1] <- current_cell_res$Fmat_current_cell
      flag[i1, 1] <- current_cell_res$flag_current_cell

      if(i1 >=4){
        for(i2 in 2:(i1 - 2)){

          current_cell_res <- nearby_Fmat_check_neighbour(current_cell = qmat[i1, i2],
                                                          up = Fmat[i1-1, i2],
                                                          left = Fmat[i1, i2-1],
                                                          upnleft = Fmat[i1-1, i2-1],
                                                          leftflag = flag[i1, i2-1],
                                                          rd = rd,
                                                          tol = tol)

          Fmat[i1, i2] <- current_cell_res$Fmat_current_cell
          flag[i1, i2] <- current_cell_res$flag_current_cell


        }
      }
    }
  }

  if(ret.flag){
    return(list(Fmat, flag))
  } else return(Fmat)
}
