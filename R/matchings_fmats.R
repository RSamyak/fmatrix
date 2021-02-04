#' @param tree
#'
#' @export matching
matching<-function(tree){
  #### Converting phylo tree to and from perfect matchings
  o <- order(tree$edge[,1], decreasing = T)
  return(tree$edge[o, 2])
}

do_mymatching <- function(matching, undo = FALSE){

  n <- length(matching)/2 + 1

  if(undo){
    ret <- matching
    ret[ret != 0] <- ret[ret != 0] + n - 1
    ret[ret == 0] <- 1:n
  }

  else{
    ret <- replace(matching, matching <= n, 0)
    ret[ret > n] <- ret[ret > n] - n + 1
  }

  ret
}



#' @param config
#'
#' @export mytree_from_mymatching
mytree_from_mymatching <- function(config){
  tree_from_matching(do_mymatching(config, undo = TRUE))
}

#' @param tree
#'
#' @export mymatching
mymatching <- function(tree){
  do_mymatching(matching(tree), undo = FALSE)
}



#' @param config
#'
#' @export tree_from_matching
tree_from_matching <- function(config){

  n<-length(config)/2 + 1

  tree<-rcoal(n,br=1)

  parents<-rep(seq(n+1,2*n-1),each=2)

  tree$edge[,1]<-rev(parents)
  tree$edge[,2]<-config

  m <- max(tree$edge)
  n <- m - tree$Nnode

  o <- which(tree$edge[,2] <= n)
  tree$edge.length[o] <- m+1-tree$edge[o,1]
  tree$edge.length[-o] <- tree$edge[-o,2] - tree$edge[-o,1]
  return(tree)
}

#' @param config
#'
#' @export proposal_matching
proposal_matching <- function(config){
  n <- length(config)/2 + 1

  ret <- config

  i <- sample(seq(1,n-2),1) # choose a pair to change

  selected <- config[2*i + 1:2]

  eligible <- selected[selected > (2*n - i) | selected <= n] # either nodes which are lower or leaves

  if(length(eligible)==1) sampleR <- eligible
  else sampleR <- sample(eligible, size=1)

  sampleL <- sample(config[2*i - 1:0], size=1)

  place1 <- which(ret == sampleL)
  place2 <- which(ret == sampleR)

  ret[place1] <- sampleR
  ret[place2] <- sampleL

  return(ret)
}


