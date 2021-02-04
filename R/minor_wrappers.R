#' @param Fmat
#'
#' @export mytree_from_F
mytree_from_F <- function(Fmat, coal.times = NULL){

  n <- ncol(Fmat) + 1
  if(is.null(coal.times)) coal.times <- seq(1, n-1)

  tree_from_F(rbind(rep(0,n),cbind(rep(0,n-1),Fmat)), coal.times)
}


#' @param Fmat
#'
#' @export num_cherries
num_cherries <- function(Fmat){
  #Simple function to return number of cherries

  n <- dim(Fmat)[1]
  c<-sum(diff(Fmat[n,])==2)
  return(c)
}

#' @param sampleF
#' @param totalF.list
#' @param dist
#'
#' @export brute.mean
brute.mean<-function(sampleF,totalF.list, dist=c("l1","l2")){
  ##Brute force frechet mean

  dist <- match.arg(dist)
  dist.matrix<-matrix(0,nrow=length(totalF.list),ncol=length(sampleF))
  if(dist == "l1") for (j in 1:length(totalF.list)){
    for (i in 1:length(sampleF)){
      dist.matrix[j,i]<-sum(abs(totalF.list[[j]]-sampleF[[i]]))
    }
  }
  if(dist == "l2") for (j in 1:length(totalF.list)){
    for (i in 1:length(sampleF)){
      dist.matrix[j,i]<-sum((totalF.list[[j]]-sampleF[[i]])**2)
    }
  }
  total<-apply(dist.matrix,1,sum)
  return(which.min(total))
}

#' @param sampleF
#' @param totalF.list
#' @param dist
#'
#' @export brute.mean.weighted
brute.mean.weighted<-function(sampleF, totalF.list, dist=c("l1","l2"), weights){
  ##Brute force frechet mean  (weighted)

  dist <- match.arg(dist)
  dist.matrix<-matrix(0,nrow=length(totalF.list),ncol=length(sampleF))
  if(dist == "l1") for (j in 1:length(totalF.list)){
    for (i in 1:length(sampleF)){
      dist.matrix[j,i]<-sum(abs(totalF.list[[j]]-sampleF[[i]]))
    }
  }
  if(dist == "l2") for (j in 1:length(totalF.list)){
    for (i in 1:length(sampleF)){
      dist.matrix[j,i]<-sum((totalF.list[[j]]-sampleF[[i]])**2)
    }
  }

  dist.matrix <- dist.matrix %*% diag(weights)

  total<-apply(dist.matrix,1,sum)
  return(which.min(total))
}


#' @param Fmat
#'
#' @export plotF
plotF <- function(Fmat, node.labels = NULL, tip.labels = NULL, ...){
  ## Wrapper to plot tree given F-matrix

  n <- ncol(Fmat) + 1
  tree<-tree_from_F(rbind(rep(0,n),cbind(rep(0,n-1),Fmat)),seq(1,n-1))
  par(mar = c(2,2,2,2))
  ape::plot.phylo(ape::ladderize(tree), direction = "downwards", show.tip.label = FALSE, ...)
  if(!is.null(node.labels)) {
    if(is.atomic(node.labels)) {if(node.labels) do.call(nodelabels, list(text = 2:n))}
    else do.call(nodelabels, node.labels)
  }

  if(!is.null(tip.labels)) {
    if(is.atomic(tip.labels)) {if(tip.labels) do.call(tiplabels, list(text = rep(0,n)))}
    else do.call(tiplabels, tip.labels)
  }
}


#' @param tree
#'
#' @export plotT
plotT <- function(tree, ...){
  ## Wrapper to plot phylo
  par(mar = c(2,2,2,2))
  ape::plot.phylo(tree,
                  # direction = "downwards",
                  show.tip.label = FALSE, ...)
  ape::axisPhylo()
}



#' @param F.list
#'
#' @export plotF.list

plotF.list <- function(F.list, nrow = NULL, ncol = NULL, colours = NULL, ...){
  n <- length(F.list)

  set <- FALSE

  if(is.null(nrow) & is.null(ncol)){
    t <- ceiling(sqrt(n))
    ncol <- t
    nrow <- ceiling(n/t)

    set <- TRUE
  } else{
    if(is.null(nrow)) nrow <- ceiling(n / ncol)
    if(is.null(ncol)) ncol <- ceiling(n / nrow)

    set <- TRUE
  }

  if(! set){
    if(nrow * ncol <= n) warning("nrow x ncol <= length(F.list)")
  }

  par(mfrow=c(nrow,ncol), mar = c(1,1,1,1))
  for(i in 1:n){
    if(!is.null(colours)) {
      plotF(F.list[[i]], edge.color = colours[i], ...)
    }
    else{
      plotF(F.list[[i]], ...)
    }
  }
  # dev.off()
}

#' @param F.list
#' @param folder
#'
#' @export plotF.list.png
plotF.list.png <- function(F.list, folder = "new", ...){
  if(!dir.exists(file.path("png", folder))) dir.create(file.path("png", folder), recursive = TRUE)
  for(i in 1:length(F.list)){
    png(paste0("png/", folder, "/file",i,".png"))
    plotF(F.list[[i]], ...)
    dev.off()
  }
}


#' @param F.list
#' @param weights
#'
#' @export meanF

meanF <- function(F.list, weights=rep(1,length(F.list))){
  ## Naive weighted mean function for F-matrices, output need not be F-matrix, but rounded F-matrix is close?

  n <- ncol(F.list[[1]]) + 1

  out <- matrix(0, n-1, n-1)
  for (i in 1:length(F.list)){
    out <- out + weights[i]*F.list[[i]]
  }
  return(out/sum(weights))
}

###' @param tree
# forget_lengths <- function(tree){
#   #### gen_Fmat for forgetting edge lengths
#
#   m <- max(tree$edge)
#   n <- m - tree$Nnode
#
#   o <- which(tree$edge[,2] <= n)
#   tree$edge.length[o] <- m+1-tree$edge[o,1]
#   tree$edge.length[-o] <- tree$edge[-o,2] - tree$edge[-o,1]
#
#   return(tree)
# }

###' @param tree
#'
# @export gen_Fmat_unweighted
# gen_Fmat_unweighted <- function(tree){
#   gen_Fmat(forget_lengths(tree))
# }


#' @param F.list
#'
#' @export matrix_list
matrix_list <- function(F.list, diag = FALSE){
  #### Better storage for distance than list?
  sapply(F.list, function(u){
    as.vector(u[lower.tri(u, diag = diag)])
  })
}


#' @param Fmat
#'
#' @export reduce_to_vector
reduce_to_vector <- function(Fmat, diag = FALSE){
  as.vector(Fmat[lower.tri(Fmat, diag = diag)])
}


#' Create a scaling vector for distances
#'
#' Create a scaling vector for distances
#'
#' @param n
#' @param type
#'
#' @export scaling_vector
scaling_vector <- function(n, type = 0){
  if(type == 0){
    t1 <- matrix(1, n-1, n-1)

    ret <- reduce_to_vector(t1)
  }

  if(type == 1){
    t1 <- 2:n %*% t(rep(1, n-1))

    ret <- reduce_to_vector(t1)
  }

  if(type == 2){
    t1 <- rep(1, n-1) %*% t(2:n)

    ret <- reduce_to_vector(t1)
  }

  return(ret)
}
