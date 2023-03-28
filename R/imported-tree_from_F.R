#' Generate tree in phylo form from F-matrix
#' Update: March 2023. Label assignment is uniform
#' This functions generates a tree in phylo format to plot from an F-matrix
#'
#' @param matF
#' @param coal_times
#'
#' @export
  tree_from_F <- function(matF, coal_times){
  #generate an ape tree (Jaime's code)
  #F is the actual form used in the code that differs from the paper's notation
  
  n= dim(matF)[1]
  edge=matrix(rep(0,6*n-6),ncol=3)
  edge[,2]= (1:(2*n-2))
  vintages = c()
  times=c(rep(0,n),coal_times)
  for (j in n:2){
    new_node = 2*n-j+1
    setleaves<-intersect(which(edge[,1]==0),1:n)
    if (length(setleaves)>1){
    next_leaf = sample(setleaves,2)}else{
      next_leaf<-setleaves
    }
    F_difference = rev(matF[,j]-matF[,j-1])
    
    if (F_difference[1]==2){
      edge[next_leaf,1]=c(new_node,new_node)
      vintages=c(vintages, new_node)
    }
    else if (F_difference[1]==1){
      selected_vintage = which(F_difference == 2)[1]+n-1
      edge[selected_vintage,1]=new_node
      if (length(next_leaf)==1){whon<-next_leaf}else{
      whon<-sample(as.vector(next_leaf),1)}
      #print(whon)
      edge[whon,1]=new_node
      vintages = c(vintages[vintages != selected_vintage],new_node)
    }
    else {
      selected_vintage1 =which(F_difference == 1)[1]+n-1
      selected_vintage2 =which(F_difference == 2)[1]+n-1
      edge[selected_vintage1,1]=new_node
      edge[selected_vintage2,1]=new_node
      vintages = vintages[vintages!=selected_vintage1]
      vintages = vintages[vintages!=selected_vintage2]
      vintages<-c(vintages,new_node)
    }
  }
  #edge[1:n,]=edge[order(edge[1:n,2]),]
  
  for (j in 1:(2*n-2)) {
    #I don't understand this
    edge[j,3]=times[edge[j,1]]-times[edge[j,2]]
  }
  edge[,1]=3*n-edge[,1]
  edge[-(1:n),2]=3*n-edge[-(1:n),2]
  
  final_tree=rcoal(n,br=coal_times)
  final_tree$edge=edge[,-3]
  final_tree$edge.length=edge[,3]
  final_tree$Nnode=n-1
  class(final_tree) <- "phylo"
  final_tree <- reorder(final_tree,"postorder")
  #final_tree$edge[final_tree$edge[, 2] <= n, 2] <- 1:n
  return(final_tree)
}

# tree_from_F <- function(matF, coal_times){
#   #generate an ape tree (Jaime's code)
#   #F is the actual form used in the code that differs from the paper's notation

#   n= dim(matF)[1]
#   edge=matrix(rep(0,6*n-6),ncol=3)
#   edge[,2]= (1:(2*n-2))
#   vintages = c()
#   times=c(rep(0,n),coal_times)
#   for (j in n:2){
#     new_node = 2*n-j+1
#     next_leaf = intersect(which(edge[,1]==0),1:n)[1]
#     F_difference = rev(matF[,j]-matF[,j-1])

#     if (F_difference[1]==2){
#       edge[next_leaf:(next_leaf+1),1]=new_node
#       vintages=c(vintages, new_node)
#     }
#     else if (F_difference[1]==1){
#       selected_vintage = which(F_difference == 2)[1]+n-1
#       edge[selected_vintage,1]=new_node
#       edge[next_leaf,1]=new_node
#       vintages = c(vintages[vintages != selected_vintage],new_node)
#     }
#     else {
#       selected_vintage1 =which(F_difference == 1)[1]+n-1
#       selected_vintage2 =which(F_difference == 2)[1]+n-1
#       edge[selected_vintage1,1]=new_node
#       edge[selected_vintage2,1]=new_node
#       vintages = vintages[vintages!=selected_vintage1]
#       vintages = vintages[vintages!=selected_vintage2]
#       vintages<-c(vintages,new_node)
#     }
#   }
#   #edge[5:8,2]=c(6,7,8,5)
#   edge[1:n,]=edge[order(edge[1:n,2]),]
#   #edge=edge[order(edge[,1]),]

#   for (j in 1:(2*n-2)) {
#     #I don't understand this
#     edge[j,3]=times[edge[j,1]]-times[edge[j,2]]
#   }
#   edge[,1]=3*n-edge[,1]
#   edge[-(1:n),2]=3*n-edge[-(1:n),2]

#   final_tree=rcoal(n,br=coal_times)
#   final_tree$edge=edge[,-3]
#   final_tree$edge.length=edge[,3]
#   final_tree$Nnode=n-1
#   class(final_tree) <- "phylo"
#   final_tree <- reorder(final_tree,"postorder")
#   final_tree$edge[final_tree$edge[, 2] <= n, 2] <- 1:n
#   return(final_tree)
# }
