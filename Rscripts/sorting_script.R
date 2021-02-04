library(fmatrix)
library(tidyverse)
library(ape)
library(phylodyn)

F.list <- F.list9

for(i in 1:length(F.list)){
  attr(F.list[[i]], "orig_index") <- i
}

encod.list <- lapply(F.list, my_encod)

for(i in 1:length(encod.list)){
  attr(encod.list[[i]], "orig_index") <- i
}


test1 <- F.list[[100]]
test2 <- F.list[[200]]


`%<myencod%` <- function(x, y) `%<=%`(my_encod(x), my_encod(y))
`%>myencod%` <- function(x, y) `%<myencod%`(y, x)

`%<col%` <- function(x, y) `%<=%`(as.vector(x[lower.tri(x)]), as.vector(y[lower.tri(y)]))
`%>col%` <- function(x, y) `%<col%`(y, x)

`%<row%` <- function(x, y) `%<=%`(as.vector(t(x)[upper.tri(x)]), as.vector(t(y)[upper.tri(y)]))
`%>row%` <- function(x, y) `%<row%`(y, x)


test1 %<myencod% test2
test1 %>myencod% test2
test1 %<col% test2
test1 %>col% test2
test1 %<row% test2
test1 %>row% test2


my_identical <- function(x, y, ...) {
  attributes(x) <- attributes(y) <- NULL
  identical(x, y, ...)
}


myrank <- function(obj, sorted.list){

  my_identical <- function(x, y, ...) {
    attributes(x) <- attributes(y) <- NULL
    identical(x, y, ...)
  }

  for(i in 1:length(sorted.list)){
    if(my_identical(obj, sorted.list[[i]])) return(i)
  }
  return(NA)
}

sorted.encod.list.encod <- bubbleSort(encod.list, `%lessthan%` = `%<=%`)
sorted.list.encod <- lapply(sorted.encod.list.encod, function(x) {
  ret <- Fmat_from_myencod(x)
  attr(ret, "orig_index") <- attr(x, "orig_index")
  return(ret)
}
)

sorted.list.row <- bubbleSort(F.list, `%lessthan%` = `%<row%`)
sorted.list.col <- bubbleSort(F.list, `%lessthan%` = `%<col%`)
# sorted.list.encod.2 <- bubbleSort(F.list, `%lessthan%` = `%<myencod%`)


check_same <- rep(FALSE, length(F.list))
for(i in 1:length(F.list)){
  if(my_identical( sorted.list.row[[i]], sorted.list.col[[i]] )) check_same[i] <- TRUE
}
which(check_same)



# myfun__test <- function(F.list, `%<1%` = `%<row%`, `%<2%` = `%<col%`){
#   for(i in 1:length(F.list)){
#     for(j in 1:i){
#       if(xor(F.list[[i]] %<1% F.list[[j]], F.list[[i]] %<2% F.list[[j]])) {
#         cat("phew!", fill = TRUE); return("diff")
#       }
#     }
#   }
#   return("same")
# }
#
# myfun__test(F.list9, `%<myencod%`, `%<row%`) ## same




# for(i in 1:length(sorted.encod.list)){
#   if(! my_identical(sorted.list.encod[[i]], sorted.list.encod.2[[i]])) {
#     cat(i); break
#   }
#   if(i == length(sorted.encod.list)) cat("done!")
# }


sorted.matching.list <- lapply(sorted.list.encod, function(x){
  ret <- matching(mytree_from_F(x))
  attr(ret, "orig_index") <- attr(x, "orig_index")
  return(ret)
}
)


orig_indices.row <- orig_indices.col <- orig_indices.encod <- c()
for(i in 1:length(F.list)){
  orig_indices.col[i] <- attr(sorted.list.col[[i]], "orig_index")
  orig_indices.row[i] <- attr(sorted.list.row[[i]], "orig_index")
  orig_indices.encod[i] <- attr(sorted.list.encod[[i]], "orig_index")
}


# i <- 1
# plotF(sorted.list.encod[[i]]); i <- i+1


# orig_indices_encod <- c()
# for(i in 1:length(sorted.encod.list.encod)){
#   orig_indices_encod[i] <- attr(sorted.encod.list.encod[[i]], "orig_index")
# }
# identical(orig_indices.encod, orig_indices_encod) ## TRUE



sorted.encod.list.row <- lapply(sorted.list.row, function(x){
  ret <- my_encod(x)
  attr(ret, "orig_index") <- attr(x, "orig_index")
  return(ret)
}
)

sorted.encod.list.col <- lapply(sorted.list.col, function(x){
  ret <- my_encod(x)
  attr(ret, "orig_index") <- attr(x, "orig_index")
  return(ret)
}
)


inv_orig_indices.row <- Matrix::invPerm(orig_indices.row)
inv_orig_indices.col <- Matrix::invPerm(orig_indices.col)
inv_orig_indices.encod <- Matrix::invPerm(orig_indices.encod)


#############################################################
{
set.seed(1)
n <- 9
init <- rcoal(n, wt = 1)

F.list <- get(paste0("F.list", n))

beta.prob <- sapply(F.list, beta_BF_likelihood, betas = 5)

n.sample <- 30
sampleF <- list()
for(i in 1:n.sample){
  sampleF[[i]] <- F.list[[sample(length(F.list), 1, prob = beta.prob)]]
}


temp <- 1000
alpha <- 0.995

config <- my_encod(gen_Fmat(init))


energy <- function(state, sample=sampleF, d="l2"){

  state.F <- Fmat_from_myencod(state)
  # state.F <- gen_Fmat(tree_from_matching(state))

  if(d == "l1") distances <- sapply(sampleF, function(x){sum(abs(x - state.F))})
  if(d == "l2") distances <- sapply(sampleF, function(x){sum((x - state.F)**2)})

  return(sum(distances))
}

probab <- function(e1, e2, temp, eps = 0.0001){
  if(e2 < e1) return(1)
  # if(e1 == e2) cat("!")
  return(exp(-(e2-e1+eps)/temp))
}



track_ranks <- c()

e1 <- energy(config)

N <- 5000
for(i in 1:N){
  new <- proposal_myencod(config)
  thisrank <- myrank(config, sorted.encod.list.encod)

  track_ranks <- append(track_ranks, thisrank)

  e2 <- energy(new)

  prob <- probab(e1, e2, temp)


  if(runif(1) <= prob){
    config <- new
    e1 <- e2

  }

  temp <- alpha*temp
}


} ## Annealing on my_encod, dist = Fmat
plot(track_ranks)  # this is tracking ranks according to %<myencod%

track_ranks_row <- sapply(track_ranks, function(x){
  inv_orig_indices.row[orig_indices.encod[x]]
})
plot(track_ranks_row) # this is tracking ranks according to %<row%

track_ranks_col <- sapply(track_ranks, function(x){
  inv_orig_indices.col[orig_indices.encod[x]]
})
plot(track_ranks_col)  # this is tracking ranks according to %<col%



#############################################################
{
set.seed(1)
n <- 9
init <- rcoal(n, wt = 1)

beta.prob <- sapply(F.list9, beta_BF_likelihood, betas = 5)

n.sample <- 30
sampleF <- list()
for(i in 1:n.sample){
  sampleF[[i]] <- F.list9[[sample(length(F.list9), 1, prob = beta.prob)]]
}


temp <- 1000
alpha <- 0.995

config <- matching(init)

energy <- function(state, sample=sampleF, d="l2"){

  # state.F <- Fmat_from_myencod(state)
  state.F <- gen_Fmat(tree_from_matching(state))

  if(d == "l1") distances <- sapply(sampleF, function(x){sum(abs(x - state.F))})
  if(d == "l2") distances <- sapply(sampleF, function(x){sum((x - state.F)**2)})

  return(sum(distances))
}

probab <- function(e1, e2, temp, eps = 0.0001){
  if(e2 < e1) return(1)
  # if(e1 == e2) cat("!")
  return(exp(-(e2-e1+eps)/temp))
}


track_ranks_matching <- c()

temp_myrank <- function(config){
  myrank(my_encod(gen_Fmat(tree_from_matching(config))), sorted.encod.list.encod)
}


N <- 5000
for(i in 1:N){
  new <- proposal_matching(config)

  thisrank <- temp_myrank(config)
  track_ranks_matching <- append(track_ranks_matching, thisrank)




  e2 <- energy(new)

  prob <- probab(e1, e2, temp)


  if(runif(1) <= prob){
    config <- new
    e1 <- e2

  }

  temp <- alpha*temp
}

} ## Annealing on matching, dist = Fmat
plot(track_ranks_matching)
length(unique(track_ranks))
length(unique(track_ranks_matching))


#############################################################
{
set.seed(1)
n <- 9
init <- rcoal(n, wt = 1)

config <- matching(init)


track_ranks_matching <- c()

temp_myrank <- function(config){
  myrank(my_encod(gen_Fmat(tree_from_matching(config))), sorted.encod.list.encod)
}


N <- 5000
for(i in 1:N){
  new <- proposal_matching(config)

  thisrank <- temp_myrank(config)
  track_ranks_matching <- append(track_ranks_matching, thisrank)



  config <- new
}
} ## Random sampling on matching
plot(track_ranks_matching)

######
{
  set.seed(1)
  n <- 9
  init <- rcoal(n, wt = 1)

  config <- my_encod(gen_Fmat(init))


  track_ranks_encod <- c()


  N <- 5000
  for(i in 1:N){
    new <- proposal_myencod(config)

    thisrank <- myrank(config, sorted.encod.list.encod)
    track_ranks_encod <- append(track_ranks_encod, thisrank)



    config <- new
  }
} ## Random sampling on my_encod
plot(track_ranks_encod)

track_ranks_col <- sapply(track_ranks_encod, function(x){
  inv_orig_indices.col[orig_indices.encod[x]]
})
plot(track_ranks_col)  # this is tracking ranks according to %<col%



