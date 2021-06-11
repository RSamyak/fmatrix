#' @param tr
#' @param tol
#'
#'
#' @export gen_hEncod

gen_hEncod <- function(tr, tol = 6) {
  ## Modified from gen_Fmat code from juliapr/phylodyn by RSamyak

  if (class(tr) != "phylo") {
    stop("The input tree must be a phylo object")
  }

  summarize_phylo_times <- fmatrix:::summarize_phylo(tr, tol = 10**(-tol))


  edge.mat <- tr$edge
  n.sample <- tr$Nnode + 1
  t.tot <- max(ape::node.depth.edgelength(tr))
  n.t <- t.tot - ape::node.depth.edgelength(tr)
  t.dat <- data.frame(
    lab.1 = edge.mat[, 1],
    lab.2 = edge.mat[, 2],
    t.end = n.t[edge.mat[, 1]],
    t.start = n.t[edge.mat[, 2]]
  )

  t.dat <- t.dat[order(t.dat$lab.2),]

  coal.t <- sort(n.t[(n.sample + 1):length(n.t)])

  n.c.event <- length(coal.t)
  tmp.s.t <- round(n.t[1:n.sample], digits = tol)

  ## This for loop is only for rounding
  for (i in 1:n.sample) {
    tmp.ind <- which(tmp.s.t == tmp.s.t[i])
    if (length(tmp.ind) > 1) {
      group.t <- min(n.t[tmp.ind])
      n.t[tmp.ind] <- group.t
      t.dat[tmp.ind, 4] <- group.t
    }
  }


  sample.t <- sort(unique(n.t[1:n.sample]))
  n.s.event <- length(sample.t)

  ## This flag should* never be triggered after the rounding
  stopifnot(n.s.event == length(unique(tmp.s.t)))


  encod.length <- n.sample + n.c.event

  # ##############################################
  # ###### This is only to get u.info to get times
  # u.t <- data.frame(
  #   t = c(coal.t, sample.t),
  #   type = c(rep("c",
  #                n.c.event), rep("s", n.s.event)),
  #   c.id = c(seq(n.c.event),
  #            rep(-9, n.s.event)),
  #   s.id = c(rep(-9, n.c.event), seq(n.s.event)),
  #   stringsAsFactors = FALSE
  # )
  #
  # u.t <- u.t[order(u.t$t, decreasing = TRUE),]
  #
  # rownames(u.t) <- paste("u", 1:(n.c.event + n.s.event), sep = ".")
  # ###### End
  # ################################################


  #### This is another u.t to get encod
  my.u.t <- data.frame(t = n.t,
                    lab = seq(1, length(n.t)),
                    type = c(rep("s", n.sample),
                             rep("c", n.c.event)))

  lab.parent <- rep(NA, nrow(my.u.t))
  for (i in 1:nrow(my.u.t)) {
    temp <- t.dat$lab.1[which(t.dat$lab.2 == i)]
    lab.parent[i] <- ifelse(length(temp) > 0, temp, NA)
  }
  my.u.t$lab.parent <- lab.parent


  my.u.t <- my.u.t[order(my.u.t$t, decreasing = TRUE),]

  lab.new <- rep(NA, nrow(my.u.t))
  lab.new[my.u.t$type == "c"] <- 2:(n.c.event + 1)
  lab.new[my.u.t$type == "s"] <- -1 * my.u.t$lab[my.u.t$type == "s"]

  my.u.t$lab.new <- lab.new


  encod <-
    rename.labels(my.u.t$lab.parent, my.u.t$lab, my.u.t$lab.new, na.replace = 1)
  encod[my.u.t$type == "s"] <- -1 * encod[my.u.t$type == "s"]

  class(encod) <- c(class(encod), "encod")


  attr(encod, "summarize_phylo_times") <- summarize_phylo_times

  times <- sort(c(summarize_phylo_times$samp_times, summarize_phylo_times$coal_times))

  attr(encod, "times") <- times

  return(encod)
}

summarize_phylo <- function (phy, ...)
{
  hgpstat <- phylodyn:::heterochronous_gp_stat(phy, ...)
  return(list(samp_times = hgpstat$samp_times, n_sampled = hgpstat$n_sampled,
              coal_times = hgpstat$coal_times))
}

# get_times_from_uinfo <- function(u.info){
#   n.coal <- max(u.info$c.id)
#   n.event.1 <- dim(u.info)[1]
#   c.match.ind.1 <- match(seq(n.coal, 1), u.info$c.id)
#   n.s.1 <- c(diff(c.match.ind.1) - 1, n.event.1 - c.match.ind.1[n.coal])
#   n.s.comb <- pmax(n.s.1, n.s.1)
#   n.a.1 <- n.s.comb - n.s.1
#   fmat.dim <- n.event.1 + sum(n.a.1)
#   comb.Fmat.1 <- matrix(0, nrow = fmat.dim, ncol = fmat.dim)
#   insert.info.1 <- phylodyn:::insert.event(e.vec = u.info$type, t.vec = u.info$t,
#                                            n.a.vec = n.a.1)
#   w.mat.1 <- phylodyn:::create.weight.mat.hetero(insert.info.1$comb.t)[-1,-1]
#   return(times_from_wt(w.mat.1))
# }

hEncod_times_fillin <- function(hEncod){

  n <- sum(hEncod > 0)

  times <- rep(NA, length(hEncod))
  coal_times <- times[hEncod > 0] <- (n):1

  m <- length(times)
  o <- is.na(times)

  count <- n

  cu <- cumsum(hEncod > 0)
  times <- n + 1 - cu
  times[hEncod < 0] <- times[hEncod < 0] - .5
  times[times == .5] <- 0

  samp_times <- times[hEncod < 0]


  n_sampled <- unname(table(samp_times))
  samp_times <- sort(unique(samp_times))
  coal_times <- sort(coal_times)

  summarize_phylo_times <- list(samp_times = samp_times,
                                n_sampled = n_sampled,
                                coal_times = coal_times
                                )

  return(summarize_phylo_times)

}


#' @param hEncod
#' @param times
#'
#'
#' @export tree_from_hEncod

tree_from_hEncod <- function(hEncod, summarize_phylo_times = NULL){


  n <- sum(hEncod > 0)

  hEncod <- preprocess_hEncod(hEncod)



  if (is.null(summarize_phylo_times)) {
    if (is.null(attr(hEncod, "summarize_phylo_times"))) {
      summarize_phylo_times <- hEncod_times_fillin(hEncod)
    }
    else{
      summarize_phylo_times <- attr(hEncod, "summarize_phylo_times")
    }
  }

  samp_times <- sort(summarize_phylo_times$samp_times, decreasing = TRUE)
  coal_times <- sort(summarize_phylo_times$coal_times, decreasing = TRUE)

  all_samp_times <- rep(summarize_phylo_times$samp_times, summarize_phylo_times$n_sampled)

  times <- sort(c(all_samp_times, coal_times), decreasing = TRUE)




  edge <- matrix(NA, nrow = 2*n, ncol = 2)
  edge.length <- rep(NA, 2*n)

  ## This is the label of the first internal node to whose parent
  ## is to be determined in the loop below
  current_internal_node <- 2 + 1
  ## Same for tip (external node)
  current_tip <- -1

  internal_times <- sort(summarize_phylo_times$coal_times, decreasing = TRUE)

  for(i in 1:nrow(edge)){

    if(hEncod[i + 1] < 0){
      edge[i, ] <- c(-hEncod[i + 1], current_tip)
      edge.length[i] <- internal_times[(-1)*hEncod[i + 1] - 1] - times[i + 1]

      current_tip <- current_tip - 1

    } else{
      edge[i, ] <- c(hEncod[i + 1], current_internal_node)
      edge.length[i] <- internal_times[hEncod[i + 1] - 1] - times[i + 1]


      current_internal_node <- current_internal_node + 1
    }
  }

  old.labels <- c((-1):(current_tip + 1), 2:(current_internal_node-1))
  new.labels <- 1:length(old.labels)

  edge[, 1] <- rename.labels(edge[, 1], old.labels, new.labels)
  edge[, 2] <- rename.labels(edge[, 2], old.labels, new.labels)


  tip.label <- paste0("t", 1:(n+1))


  tr <- list(
    edge = edge,
    edge.length = edge.length,
    Nnode = n,
    tip.label = tip.label
  )

  class(tr) <- "phylo"

  return(tr)
}

#' @param hEncod
preprocess_hEncod <- function(hEncod){

  abs_hEncod <- abs(hEncod)

  n <- sum(hEncod > 0) + 2

  sum_neg <- sum(hEncod < 0)

  for(i in 2:(n-1)){
    count <- sum(abs_hEncod == i)
    if(count == 2) {
      next
    }
    else if(count == 1) {
      hEncod <- append(hEncod, -i)
    }
    else if(count == 0) {
      hEncod <- append(hEncod, rep(-i, 2))
    }
    else {
      stop(sprintf("%d appears %d times in hEncod, must be 0, 1, or 2", i, count))
    }
  }


  return(hEncod)

}

#' @param str
#' @param ref1
#' @param ref2
#' @param na.replace
rename.labels <- function(str, ref1, ref2, na.replace = NA){
  temp <- ref2[match(str, ref1)]
  temp[is.na(temp)] <- na.replace

  temp
}


#' @param hEncod
#' @param check          Do we check if the input myencod is valid?
#'                           Default is FALSE
#'
#' @export proposal_hEncod
proposal_hEncod <- function(hEncod, check = FALSE, digits.tol = 6){

  ## Note: This chain is not irreducible.  Need to fix.

  hEncod <- round(hEncod, digits = digits.tol)

  abs_hEncod <- abs(hEncod)
  negative_entries <- which(hEncod < 0)

  # if(check) stopifnot(is.myencod(myencod))

  n <- (length(hEncod) + 1) / 2


  # if(n <= 3) return(hEncod)

  ret <- hEncod

  i <- sample.vec(seq(2, length(hEncod)), 1) # choose an index to change

  if(hEncod[i] > 0){

    temp <- hEncod
    temp[i] <- Inf
    temp <- temp[temp > 0]
    r <- which(temp == Inf)

    eligible <- 2:r

    part_to_check <- abs_hEncod[-c(i, negative_entries[negative_entries > i])]
    tab <- table(part_to_check)

    ineligible <- names(tab)[tab >= 2]

    eligible <- setdiff(eligible, as.numeric(ineligible))

    if(length(eligible) == 0) browser()

    proposed <- sample.vec(eligible, 1)

    ret[i] <- proposed



    temp <- which(hEncod == (-proposed))
    temp <- temp[length(temp)]

    if(length(temp) == 0) {
      return(ret)
    } else{

      if(is.na(temp)) browser()

      if(TRUE){
        t1 <- sum(hEncod[1:temp] > 0) + 1
        if(! abs(hEncod[i]) <= t1) browser()
      }

      ret[temp] <- -hEncod[i]




      if (length(ret) != length(hEncod))
        browser()

      # if(check) stopifnot(is.myencod(ret))
      # if (check) {
      #   if (!is.myencod(ret))
      #     browser()
      # }
      return(ret)
    }



  }
  else {

    ### This is not the best chain here.


    ## This is the maximum abs(hEncod[i]) could have been
    ri <- sum( hEncod[1:i] > 0 ) + 1

    j <- sample.vec(negative_entries, 1)
    ## This is the maximum abs(hEncod[i]) could have been
    rj <- sum( hEncod[1:j] > 0 ) + 1

    if( (abs_hEncod[j] <= ri) &
        (abs_hEncod[i] <= rj) ){
      ret[j] <- hEncod[i]
      ret[i] <- hEncod[j]
    }

    return(ret)
  }


  return(-1)
}

#' Check if a string is a valid hEncod
#'
#' Check if a string correctly encodes a heterochronous ranked tree shape
#'
#' @param hEncod
#'
#' @export is.hEncod
is.hEncod <- function(hEncod, digits.tol = 6){
  hEncod <- round(hEncod, digits = digits.tol)

  abs_hEncod <- abs(hEncod)
  sign_hEncod <- (hEncod > 0)
  negative_entries <- which(hEncod < 0)

  # if(check) stopifnot(is.myencod(myencod))

  n <- (length(hEncod) + 1) / 2

  ## First entry must be 1 (root condition)
  if(! hEncod[1] == 1) return(FALSE)

  ## There should be exactly n leaves, i.e. negative entries
  if(! sum(!sign_hEncod) == n) return(FALSE)

  ## Each internal node should appear exactly twice
  if(! all(table(abs_hEncod[-1]) == 2)) return(FALSE)


  ## No internal node should occur before it was created
  num_branch_events <- 1
  for(i in 1:length(hEncod)){

    if(abs_hEncod[i] > num_branch_events) return(FALSE)
    if(1*sign_hEncod[i] > 0) num_branch_events <- num_branch_events + 1
  }

  return(TRUE)
}


#' @param hEncod
#' @param check          Do we check if the input myencod is valid?
#'                           Default is FALSE
#'
#' @export proposal_hEncod_new
proposal_hEncod_new <- function(hEncod, check = FALSE, digits.tol = 6, max.loop = 1e3){

  ## Note: This chain is not irreducible.  Need to fix.

  hEncod <- round(hEncod, digits = digits.tol)

  abs_hEncod <- abs(hEncod)
  sign_hEncod <- (hEncod > 0)
  negative_entries <- which(hEncod < 0)

  if(check) stopifnot(is.hEncod(hEncod))

  n <- (length(hEncod) + 1) / 2


  # if(n <= 3) return(hEncod)

  ret <- hEncod

  for(loop.ind in 1:max.loop){

    # choose two indexes to change for abs_hEncod
    i.t <- sample.vec(seq(2, length(hEncod)), 1)
    j.t <- sample.vec(setdiff(seq(2, length(hEncod)), i.t), 1)

    ret <- abs_hEncod
    ret[j.t] <- abs_hEncod[i.t]
    ret[i.t] <- abs_hEncod[j.t]


    # choose two indexes to change for sign_hEncod
    # i.s <- sample.vec(seq(2, length(hEncod)), 1)
    # j.s <- sample.vec(setdiff(seq(2, length(hEncod)), i.s), 1)

    sign_ret <- sign_hEncod
    # sign_ret[j.s] <- sign_hEncod[i.s]
    # sign_ret[i.s] <- sign_hEncod[j.s]

    sign_ret <- 2*sign_ret - 1

    ret <- ret * sign_ret

    if(is.hEncod(ret)) return(ret)
  }

  return(hEncod)
}
