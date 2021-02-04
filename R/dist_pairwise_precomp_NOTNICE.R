dist_pairwise_precomp <- function (tmp.info.1, tmp.fmat.1, tmp.info.2, tmp.fmat.2, dist.method = "l1", weighted = FALSE)
{
  if (!(dist.method == "l1" | dist.method == "l2")) {
    stop("dist.method should be either \"l1\" or \"l2\"")
  }
  # if (tr.1$Nnode != tr.2$Nnode) {
  #   stop("Two trees must have the same number of taxa.")
  # }
  # tr.dat.1 <- gen.tr.data(tr.1)
  # tr.dat.2 <- gen.tr.data(tr.2)
  # tmp.info.1 <- tr.dat.1$u.info
  # tmp.info.2 <- tr.dat.2$u.info
  # tmp.fmat.1 <- tr.dat.1$Fmat
  # tmp.fmat.2 <- tr.dat.2$Fmat
  n.coal <- max(tmp.info.1$c.id)
  stopifnot(n.coal == max(tmp.info.2$c.id))
  n.event.1 <- dim(tmp.info.1)[1]
  n.event.2 <- dim(tmp.info.2)[1]
  c.match.ind.1 <- match(seq(n.coal, 1), tmp.info.1$c.id)
  c.match.ind.2 <- match(seq(n.coal, 1), tmp.info.2$c.id)

  n.s.1 <- c(diff(c.match.ind.1) - 1, n.event.1 - c.match.ind.1[n.coal])
  n.s.2 <- c(diff(c.match.ind.2) - 1, n.event.2 - c.match.ind.2[n.coal])
  n.s.comb <- pmax(n.s.1, n.s.2)
  n.a.1 <- n.s.comb - n.s.1
  n.a.2 <- n.s.comb - n.s.2
  stopifnot((n.event.1 + sum(n.a.1)) == (n.event.2 + sum(n.a.2)))
  fmat.dim <- n.event.1 + sum(n.a.1)
  comb.Fmat.1 <- matrix(0, nrow = fmat.dim, ncol = fmat.dim)
  comb.Fmat.2 <- matrix(0, nrow = fmat.dim, ncol = fmat.dim)
  insert.info.1 <- phylodyn:::insert.event(e.vec = tmp.info.1$type, t.vec = tmp.info.1$t,
                                n.a.vec = n.a.1)
  insert.info.2 <- phylodyn:::insert.event(e.vec = tmp.info.2$type, t.vec = tmp.info.2$t,
                                n.a.vec = n.a.2)
  a.rind.1 <- sort(insert.info.1$fake.ind, decreasing = TRUE)
  a.cind.1 <- insert.info.1$fake.ind + 1
  for (rind in dim(tmp.fmat.1)[1]:2) {
    for (cind in rind:2) {
      rind.new <- insert.info.1$ind.map[rind, 2]
      cind.new <- insert.info.1$ind.map[cind - 1, 2] +
        1
      comb.Fmat.1[rind.new, cind.new] <- tmp.fmat.1[rind,
                                                    cind]
    }
  }
  for (rind in a.rind.1) {
    comb.Fmat.1[rind, 2:rind] <- comb.Fmat.1[rind + 1, 2:rind]
  }
  for (cind in a.cind.1) {
    comb.Fmat.1[cind:fmat.dim, cind] <- comb.Fmat.1[cind:fmat.dim,
                                                    cind - 1]
  }
  stopifnot(all(comb.Fmat.1[upper.tri(comb.Fmat.1)] == 0))
  a.rind.2 <- sort(insert.info.2$fake.ind, decreasing = TRUE)
  a.cind.2 <- insert.info.2$fake.ind + 1
  for (rind in dim(tmp.fmat.2)[1]:2) {
    for (cind in rind:2) {
      rind.new <- insert.info.2$ind.map[rind, 2]
      cind.new <- insert.info.2$ind.map[cind - 1, 2] +
        1
      comb.Fmat.2[rind.new, cind.new] <- tmp.fmat.2[rind,
                                                    cind]
    }
  }
  for (rind in a.rind.2) {
    comb.Fmat.2[rind, 2:rind] <- comb.Fmat.2[rind + 1, 2:rind]
  }
  for (cind in a.cind.2) {
    comb.Fmat.2[cind:fmat.dim, cind] <- comb.Fmat.2[cind:fmat.dim,
                                                    cind - 1]
  }
  stopifnot(all(comb.Fmat.2[upper.tri(comb.Fmat.2)] == 0))
  if (weighted) {
    w.mat.1 <- phylodyn:::create.weight.mat.hetero(insert.info.1$comb.t)
    w.mat.2 <- phylodyn:::create.weight.mat.hetero(insert.info.2$comb.t)
    if (dist.method == "l1") {
      dist <- sum(abs(comb.Fmat.1 * w.mat.1 - comb.Fmat.2 *
                        w.mat.2))
    }
    else if (dist.method == "l2") {
      dist <- sqrt(sum((comb.Fmat.1 * w.mat.1 - comb.Fmat.2 *
                          w.mat.2)^2))
    }
    else {
      stop("unsupported dist.method")
    }
  }
  else {
    if (dist.method == "l1") {
      dist <- sum(abs(comb.Fmat.1 - comb.Fmat.2))
    }
    else if (dist.method == "l2") {
      dist <- sqrt(sum((comb.Fmat.1 - comb.Fmat.2)^2))
    }
    else {
      stop("unsupported dist.method")
    }
  }
  return(dist)
}
