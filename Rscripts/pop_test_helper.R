gen.tr.data <- function(tr, tol=13) {
      ## ================================================================
      # Get depth of each coalescent node and relabel
      # The time starts from zero at the first sampling event (u_{n+m-1}).
      # This is a revised version to correct for the sampling time numerical issue.
      # Instead of rounding up all edge and branch lengths, it clusters
      # sampling events occuring within tol into one sampling event.
      #
      ## Input:
      #   tr: phylo object tree.
      #   tol: numerical precision to round up time at each node to avoid
      #           the nodes sampled at the same time is treated as sampled
      #           at different times due to numerical error.
      #
      ## Output:
      #   list object containing the following
      #   Fmat: F-matrix of the input tr
      #   u.info: data.frame with time, event type, and other info of the tree
      #
      # Note that the tip labels are always labled with 1:n.tip and
      # the internal nodes are (n.tip+1):(2*n.tip-1) in phylo.
      
      if (class(tr) != 'phylo') {
            stop('The input tree must be a phylo object')
      }
      
      edge.mat <- tr$edge
      n.sample <- tr$Nnode + 1
      
      t.tot <- max(ape::node.depth.edgelength(tr))
      n.t <-  t.tot - ape::node.depth.edgelength(tr)
      
      edge.mat <- tr$edge
      t.dat <- data.frame(lab.1=edge.mat[,1], lab.2=edge.mat[,2],
                          t.end=n.t[edge.mat[,1]], t.start=n.t[edge.mat[,2]])
      t.dat <- t.dat[order(t.dat$lab.2), ]
      
      # coalescent times
      coal.t <- sort(n.t[(n.sample+1) : length(n.t)])
      n.c.event <- length(coal.t) #number of coalescent event
      
      if (any(diff(coal.t) == 0)) {
            stop('more than one coalescent event at a given time.')
      }
      
      # sampling times
      # correct for numerical issue in sampling time
      tmp.s.t <- round(n.t[1 : n.sample], digits=tol)
      for (i in 1:n.sample) {
            tmp.ind <- which(tmp.s.t == tmp.s.t[i])
            if (length(tmp.ind) > 1) {
                  group.t <- min(n.t[tmp.ind]) # replace with min of grouped sampling time
                  n.t[tmp.ind] <- group.t
                  t.dat[tmp.ind, 4] <- group.t
            }
      }
      sample.t <- sort(unique(n.t[1 : n.sample]))
      n.s.event <- length(sample.t) # number of sampling events
      stopifnot(n.s.event == length(unique(tmp.s.t)))
      
      
      # combined time points for F-matrix
      u.t <- data.frame(t=c(coal.t, sample.t),
                        type=c(rep('c', n.c.event), rep('s', n.s.event)),
                        c.id=c(seq(n.c.event), rep(-9, n.s.event)),
                        s.id=c(rep(-9, n.c.event), seq(n.s.event)),
                        stringsAsFactors=FALSE)
      u.t <- u.t[order(u.t$t, decreasing=TRUE),]
      rownames(u.t) <- paste('u', 1:(n.c.event+n.s.event), sep='.')
      
      ## Construct F matrix
      # n.col = n.row = (# sampling event) + (# of coalescent event) - 1
      # F(i,j) = # of lineages that exist and do not coalesce in (u_j, u_{j-1}).
      
      f.dim <- n.s.event + n.c.event
      Fmat <- matrix(0, nrow=f.dim, ncol=f.dim)
      
      for (i in 2:f.dim) {
            for (j in 2:i) {
                  u.start <- u.t$t[i]
                  u.end <- u.t$t[j-1]
                  Fmat[i,j] <- sum((t.dat$t.end >= u.end) & (t.dat$t.start <= u.start))
            }
      }
      
      return(list(Fmat=Fmat, u.info=u.t))
}


create.wmat <- function(u.t) {
      # This function creates a matrix, W, for the branch length weight
      # to be multiplied (element-wise) to the F-matrix in distance metric
      # computation.
      # Input:
      #   u.t: vector of times for each event
      #        (in decreasing order with current sampling time being 0)
      # Output:
      #    For a tree with n taxa, W is a n x n matrix with entries:
      #   upper tri: all zeroes
      #   first col: all zeroes
      #   W_{ij}: u[j-1] - u[i] (2 <= i,j <= n)
      
      n.dim <- length(u.t)
      col.vec <- c(0, u.t[-n.dim])
      row.vec <- c(0, u.t[-1])
      
      w.mat <- matrix(rep(col.vec, each=n.dim), nrow=n.dim)
      w.mat <- w.mat - matrix(rep(row.vec, times=n.dim), nrow=n.dim)
      w.mat[upper.tri(w.mat)] <- 0
      w.mat[,1] <- 0
      
      return(w.mat)
}

create.wmat2 <- function(u.t) {
   # This function creates a matrix, W, for the branch length weight
   # to be multiplied (element-wise) to the F-matrix in distance metric
   # computation.
   # Input:
   #   u.t: vector of times for each event
   #        (in decreasing order with current sampling time being 0)
   # Output:
   #    For a tree with n taxa, W is a n x n matrix with entries:
   #   upper tri: all zeroes
   #   first col: all zeroes
   #   W_{ij}: u[j-1] - u[i] (2 <= i,j <= n)
   
   n.dim <- length(u.t)
   col.vec <- -c(0,diff(u.t))
   
   w.mat <- matrix(rep(col.vec, each=n.dim), nrow=n.dim)
   w.mat[upper.tri(w.mat)] <- 0
   w.mat[,1] <- 0
   
   return(w.mat)
}



create.dmat.same.s <- function(Fvec.mat.1, Wvec.mat.1,
                               Fvec.mat.2, Wvec.mat.2,
                               dist.method='l1', weighted=F) {
      ## This function computes distance matrix between heterochronous trees 
      #  with the same number of sampling event. Since all trees have the same 
      #  number of sampling event and the number of taxa, all trees share the
      #  same F-matrix dimension. No time-subdivision step is necessary.
      #
      ## This function creates a distance matrix D, where the entry D[i,j]
      #  represents the distance between i-th and j-th trees.
      #
      ## Input: 
      #       Fvec.mat.1, Wvec.mat.1: Matrices of F-matrix & W-matrix vectors for group 1
      #       Fvec.mat.2, Wvec.mat.2: Matrices of F-matrix & W-matrix vectors for group 2
      #       n.tr: number of trees in each group
      #       dist.method: 'l1' = sum_{ij}(abs(F1[i,j] - F2[i,j]))
      #                    'l2' = sqrt(sum_{ij}(F1[i,j] - F2[i,j])^2)
      #       weighted: T/F, weighted by time interval when computing dist.pairwise
      ## Output:
      #       Dmat: distance matrix with entries Dmax[i,j] denotes distance
      #             between tree i and tree j (represented as F-matrix). 
      #             By definition, this is a symmetric matrix with zero diagonal.
      #
      # Fvec.mat, Wvec.mat: 
      #   col: trees
      #   row: Fvec or Wvec
      
      # ===== group.1 =====
      if (weighted) {
            mat.1 <- Fvec.mat.1 * Wvec.mat.1
      } else {
            mat.1 <- Fvec.mat.1
      }
      rm(Fvec.mat.1, Wvec.mat.1)
      n.tr.1 <- dim(mat.1)[2]
      
      # ===== group.2 =====
      if (weighted) {
            mat.2 <- Fvec.mat.2 * Wvec.mat.2
      } else {
            mat.2 <- Fvec.mat.2
      }
      rm(Fvec.mat.2, Wvec.mat.2)
      n.tr.2 <- dim(mat.2)[2]
      
      # ===== construct distance matrix =====
      # row: group.1, col: group.2
      Dmat <- matrix(NA, nrow=n.tr.1, ncol=n.tr.2)
      
      for (i in 1:n.tr.1) {
            if (dist.method == 'l1') {
                  # L1,1 norm
                  dist <- colSums(abs(mat.2 - mat.1[ , i]))
            } else if (dist.method == 'l2') {
                  # Frobenius norm
                  dist <- sqrt(colSums((mat.2 - mat.1[ , i])^2))
            } else {
                  stop('unsupported dist.method')
            }
            
            Dmat[i, ] <- dist
            # rm(dist)
      }
      
      stopifnot(!any(is.na(Dmat)))
      
      return(Dmat)
}


gen.lower.tri.ind <- function(n) {
      ## This function generateds indices of lower triangle of a matrix of
      #  size n. 
      ## Input: 
      #   n: dimension of a matrix
      z <- sequence(n)
      return(cbind(
            row = unlist(lapply(2:n, function(x) x:n), use.names = FALSE),
            col = rep(z[-length(z)], times = rev(tail(z, -1))-1)))
}


# ============================================================
# construct pairwise comparison vec
# ============================================================

patch.dmat <- function(n.model, n.sim, save.dir,
                       dist.method, weighted) {
      comp.ind <- rbind(gen.lower.tri.ind(n.model),
                        matrix(rep(1:n.model, each=2), ncol=2, byrow=T))
      comp.ind <- comp.ind[order(comp.ind[,1], comp.ind[,2]), ]
      n.comp <- dim(comp.ind)[1]
      
      #cl <- makeCluster(detectCores(), outfile='log.txt')
      cl <- parallel::makeCluster(12, setup_strategy = "sequential")
      registerDoParallel(cl)
      
      dmat.list <- foreach(j = 1:n.comp, .export=c('create.dmat.same.s')) %dopar% {
            ind.1 <- comp.ind[j, 1]
            ind.2 <- comp.ind[j, 2]
            
            # group 1 
            load(paste(save.dir, ind.1, '_FWvec_mat.RData', sep=''))
            Fvec.mat.1 <- Fvec.mat
            Wvec.mat.1 <- Wvec.mat
            rm(Fvec.mat, Wvec.mat)
            
            # group 2
            load(paste(save.dir, ind.2, '_FWvec_mat.RData', sep=''))
            Fvec.mat.2 <- Fvec.mat
            Wvec.mat.2 <- Wvec.mat
            rm(Fvec.mat, Wvec.mat)
            
            # compute dmat
            create.dmat.same.s(Fvec.mat.1=Fvec.mat.1, Wvec.mat.1=Wvec.mat.1,
                               Fvec.mat.2=Fvec.mat.2, Wvec.mat.2=Wvec.mat.2,
                               dist.method=dist.method, weighted=weighted)
      }
      
      stopCluster(cl)
      
      dmat.dim <- n.model * n.sim
      Dmat <- matrix(NA, nrow=dmat.dim, ncol=dmat.dim)
      
      for (i in 1:n.comp) {
            ind.1 <- comp.ind[i, 1]
            ind.2 <- comp.ind[i, 2]
            range.1 <- (1 + (ind.1-1)*n.sim):(ind.1*n.sim)
            range.2 <- (1 + (ind.2-1)*n.sim):(ind.2*n.sim)
            
            Dmat[range.1, range.2] <- dmat.list[[i]]
      }
      
      tmp.dist <- Dmat[lower.tri(Dmat)]
      stopifnot(!any(is.na(tmp.dist)))
      Dmat <- t(Dmat)
      Dmat[lower.tri(Dmat)] <- tmp.dist
      
      stopifnot(!any(is.na(Dmat)))
      
      return(Dmat)
}


find.medoid <- function(dmat) {
      # This function find an index of a medoid tree. The indices of trees follow
      # the same ordering as in the input distance matrix dmat. 
      # Here, which is used instead of which.min in case there're more than one
      # index corresponding to the minimum value. This might have very small
      # probability if we include branch length. 
      total.dist <- apply(dmat, 1, sum)
      return(which(total.dist == min(total.dist)))
}

find.frechet.mean <- function(dmat) {
      # This function find a sample frechet mean of the given sampled trees. 
      total.dist <- apply(dmat^2, 1, sum)
      return(which(total.dist == min(total.dist)))
}


addalpha <- function(colors, alpha=1.0) {
      r <- col2rgb(colors, alpha=T)
      # Apply alpha
      r[4,] <- alpha*255
      r <- r/255.0
      return(rgb(r[1,], r[2,], r[3,], r[4,]))
}


plot.mds <- function(n.model, n.sim, mds.points, m.pts, plt.legend=NA,
                     col.set='Dark2',
                     plt.legend.loc='bottomright') {
      
      if (n.model <= 2) {
            cols <- addalpha(brewer.pal(5, col.set), 0.5)
            cols.sol <- addalpha(brewer.pal(5, col.set), 1.0)
      } else if (n.model <= 8) {
            cols <- addalpha(brewer.pal(n.model, col.set), 0.5)
            cols.sol <- addalpha(brewer.pal(n.model, col.set), 1.0)
      } else if (n.model > 9) {
            qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
            col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
            cols.tmp <- sample(col_vector, n.model)
            cols <- addalpha(cols.tmp, 0.5)
            cols.sol <- addalpha(cols.tmp, 1.0)
      }
      
      par(pty='s')
      plot(mds.points[1:n.sim,],
           col=cols[1], pch=16, xlab='MDS axis 1', ylab='MDS axis 2',
           main=paste('n.sim =', n.sim, 'n.tip =', n.tip),
           xlim=c(min(mds.points[,1]), max(mds.points[,1])),
           ylim=c(min(mds.points[,2]), max(mds.points[,2])))
      for (i in 2:n.model) {
            plt.ind <- ((i-1)*n.sim + 1):(i*n.sim)
            points(mds.points[plt.ind, ], col=cols[i], pch=16)
      }
      
      # plot medoid points
      for (i in 1:n.model) {
            points(mds.points[m.pts[[i]], 1], mds.points[m.pts[[i]], 2],
                   bg=cols.sol[i], pch=24, cex=2, col="black", lwd = 2)
      }
      
      if (n.model <= 8) {
            legend(plt.legend.loc, legend=plt.legend, col=cols.sol, pch=17)
      }
      
      return(cols)
}
