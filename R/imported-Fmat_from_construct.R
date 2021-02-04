#' @param n
#' @param output
#' @param relpos
#'
#' @export
Fmat_from_construct <- function(n, output, relpos){
  #This function converts the output of the construct_trees function to

  # Expected syntax:
  # F.list<-list()
  # for (j in 1:nrow(res$res)){
  #   F.list[[j]]<- Fmat_from_construct(n, res$res[j,], res$F.matrixRel)
  # }
  out <- matrix(0,nrow=n-1,ncol=n-1)
  diag(out) <- seq(2,n)
  for (i in 1:length(output)){
    coord <- relpos[relpos[,1]==i,]
    out[(coord[3]+1),coord[2]] <- output[i]
  }
  return(out)
}
