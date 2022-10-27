

expand_strings <- function(single){
  if(single[1] != 1){single <- c(1, single)}
  
  n <- length(single) + 1
  
  counts <- table(single)
  names <- as.numeric(names(counts))
  
  othernames <- setdiff(1:n, names)
  othercounts <- rep(0, length(othernames))
  counts <- c(counts, othercounts)
  names <- c(names, othernames)
  
  counts <- counts[order(names)]
  
  allowed <- (1:n)[counts < 2][-1]

  lapply(allowed, function(i){append(single, i)[-1]})
}

expand_strings_all <- function(smaller_list){
  ret <- list()
  
  for(i in 1:length(smaller_list)){
    ret <- append(ret, expand_strings(smaller_list[[i]]))
  }
  
  ret
  
}

construct_strings <- function(n){
  
  this <- list(c(1))
  
  if(n == 2){
    return(this)
  }
  
  for(i in 3:n){
    this <- expand_strings_all(this)
  }
  
  this
}

# Ns <- 3:12
# 
# prev <- tlist2 <- list(c(1))
# 
# for(n in Ns){
#   
#   this <- expand_strings_all(prev)
#   
#   assign(paste0("tlist", n), this)
#   
#   prev <- this
# }
# 
# rm(prev, this, n)
#
# save(list = ls()[grep("^tlist",ls())], file = "tlists.RData")


