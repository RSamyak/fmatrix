library(fmatrix)

weight_matr <- function(u.t){

  n <- length(u.t) - 1

  (rep(1, n) %*% t(u.t[-length(u.t)]) ) - ((u.t[-1]) %*% t(rep(1,n)))

}

matr_list <- function(list, info.list){

  k <- length(list)

  sapply(1:k, function(u){
    info <- info.list[[u]]
    Fmat <- list[[u]]
    ret <- Fmat[-1, -1]
    weights <- weight_matr(info$t)

    reduce_to_vector(ret * weights, diag = TRUE)
  })

}

# matrix_F <- matr_list(Fmat_sublist, Finfo_sublist)
# matrix_F <- matr_list(Fmat_sublist[1:10], Finfo_sublist[1:10])

energy <- function(state, matr = matrix_F, d="l2"){

  n <- ncol(matr)

  tr <- tree_from_hEncod(state)
  tr.dat <- phylodyn:::gen.tr.data(tr)

  state.F <- tr.dat$Fmat[-1, -1] * weight_matr(tr.dat$u.info$t)
  state.F <- reduce_to_vector(state.F, diag = TRUE)

  if(d == "l1") distances <-
    sum(abs( (matr - state.F)
             # %*% diag(1/scaling_vector(n, type = type))
    ))/n
  if(d == "l2") distances <-
    sqrt(sum(( (matr - state.F)
               # %*% diag(1/scaling_vector(n, type = type))
    )**2)/n)

  return(mean(distances))
}

# energy(hEncod_sublist[[5]])

probab <- function(e1, e2, temp, eps = 0.0001){
  if(e2 < e1) return(1)
  # if(e1 == e2) cat("!")
  return(exp(-(e2-e1+eps)/temp))
}
