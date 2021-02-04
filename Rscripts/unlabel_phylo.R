unlabel_phylo <- function(phylo){
  phylo$tip.label <- paste0('t', 1:length(phylo$tip.label))
  phylo
}

unlabel_multiPhylo <- function(multiPhylo){
  out <- lapply(multiPhylo, unlabel_phylo)
  class(out) <- class(multiPhylo)
  out
}

temp2 <- unlabel_multiPhylo(temp$data_list$trees_sublist)


library(phangorn)
densiTree(temp2)
