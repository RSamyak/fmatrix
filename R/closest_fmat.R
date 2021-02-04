#' @param qmat
#' @param F.list
#'
#' @export

closest_Fmat <- function(qmat, F.list){
  ## Closest Fmat, given full F.list

  if(is.Fmat(qmat)) return(qmat) ## Note this is assuming qmat is in F.list
  return(F.list[[brute.mean(list(qmat),F.list)]])
}
