library(parallel)
library(doParallel)
library(foreach)
library(RColorBrewer)

rm(list=ls())

base.dir <- '~/Desktop/TreeMetrics/R/for_julia/'
setwd(base.dir)
data.dir <- './data/'
save.dir <- data.dir
source('./pop_test_helper.R')

# load('~/Desktop/TreeMetrics/R/for_paper/5_1_2_demographic/data/sim_trees_nsim_1000_ntip_100.RData')
# all.tr <- unlist(tr.list, recursive=F)
# rm(tr.info.list)

treetypes <- c('uniform', 'exponential', 'logistic')

n.sim <- 1000
n.model <- 3
n.tip <- 100

# Correct the numerical issue at the tips
corr.flag <- F
if (corr.flag) {
      all.tr.info <- lapply(all.tr, gen.tr.data, tol=8)
      stopifnot(all(sapply(1:length(all.tr.info), 
                           function(i) dim(all.tr.info[[i]][[1]])[1]) == 100))
      tr.info.list <- list()
      tr.info.list[[1]] <- all.tr.info[1:1000]
      tr.info.list[[2]] <- all.tr.info[1001:2000]
      tr.info.list[[3]] <- all.tr.info[2001:3000]
      save(tr.list, all.tr.info, tr.info.list,
           file=paste(data.dir, 'pophist_corr.RData', sep=''))
} else {
      load(paste(data.dir, 'pophist_corr.RData', sep=''))
}

# =============================================================
# convert Fmat in tr.info.list to Fvec and contrust Fvec matrix
# and Wmat to Wvec and construct Wvec matrix
# =============================================================
## For each tree, convert F-matrix into F-vector and Weight matrix into
#  W-vector by taking relevant entries (non-upper-diagonal) only.
#  Each col: each tree

fmat.dim <- 100
create.FWvec <- T

if (create.FWvec) {
      Fvec.length <- fmat.dim * (fmat.dim - 1) /2
      
      for (i in 1:n.model) {
            gene <- treetypes[i]
            print(paste('processing model =', i, 'out of', n.model))
            Fvec.mat <- matrix(NA, nrow=Fvec.length, ncol=n.sim)
            Wvec.mat <- matrix(NA, nrow=Fvec.length, ncol=n.sim)
            
            for (j in 1:n.sim) {
                  tmp.fmat <- tr.info.list[[i]][[j]]$Fmat
                  tmp.wmat <- create.wmat2(tr.info.list[[i]][[j]]$u.info$t)
                  stopifnot(all(tmp.fmat[1, ] == 0) & all(tmp.fmat[ , 1] == 0))
                  stopifnot(all(tmp.wmat[1, ] == 0) & all(tmp.wmat[ , 1] == 0))
                  tmp.fmat <- tmp.fmat[-1, -1]
                  tmp.wmat <- tmp.wmat[-1, -1]
                  
                  Fvec.mat[ , j] <- tmp.fmat[!upper.tri(tmp.fmat)]
                  Wvec.mat[ , j] <- tmp.wmat[!upper.tri(tmp.wmat)]
                  rm(tmp.fmat, tmp.wmat)
            }
            
            stopifnot(!any(is.na(Fvec.mat)))
            stopifnot(!any(is.na(Wvec.mat)))
            
            vec.mat.save.f <- paste(save.dir, i, '_FWvec_mat.RData', sep='')
            save(gene, Fvec.mat, Wvec.mat, file=vec.mat.save.f)
            
            print(format(object.size(Fvec.mat), units='Mb'))
            print(format(object.size(Wvec.mat), units='Mb'))
            
            rm(gene, Fvec.mat, Wvec.mat)
      }
}

rm(tr.info.list)


# ============================================================
# MDS
# ============================================================
mds.analysis <- function(dmat, save.f.pre) {
      # dmat: distance matrix
      # save.f.pre: prefix for same names for the figures and data save
      
      ## MDS
      mds <- cmdscale(dmat, eig=T, k=3)
      
      ## Find medoid points for each model and plot
      m.pts <- list()
      for (i in 1:n.model) {
            tr.ind <- ((i-1)*n.sim + 1):(i*n.sim)
            m.pts[[i]] <- find.frechet.mean(dmat[tr.ind, tr.ind]) + (i-1)*n.sim
      }
      
      png(paste(save.f.pre, 'mds.png', sep=''), width=800, height=800)
      par(mfrow=c(1,1), pty="s", mar=c(4,4,2,2), cex=1.5)
      plot.mds(n.model, n.sim, mds$points[,1:2], m.pts,  
               plt.legend=treetypes, plt.legend.loc='topright')
      dev.off()
      
      ## First 10 Eigenvalues and Goodness of Fit
      png(paste(save.f.pre, 'mds_gof.png', sep=''), width=1000, height=500)
      par(mfrow=c(1,2), pty="s", mar=c(4,4,2,2), cex=1.5)
      plot(mds$eig[1:10]/sum(mds$eig), type="h", lwd=3, las=1, col='blue',
           xlab="Number of dimensions", ylab="Eigenvalues")
      plot(mds$eig/sum(mds$eig), type="l", lwd=3, las=1, col='blue', 
           xlab="Number of dimensions", ylab="Eigenvalues")
      dev.off()
      print(mds$GOF)
      
      ## Number of negative eigenvalues
      neg.eig.ind <- which(mds$eig < 0)
      neg.eig <- mds$eig[neg.eig.ind]
      print(paste('number of negative eigenvalues:', length(neg.eig.ind), 
                  'of', length(mds$eig)))
      print(paste('(scaled) max', max(neg.eig)/sum(mds$eig), 
                  'and (scaled) min', min(neg.eig)/sum(mds$eig)))    
      
      ## Plot medoid trees.
      png(paste(save.f.pre, 'medoid_tr.png', sep=''), width=1600, height=1000)
      par(mfrow=c(2,4), oma=c(1,1,1,1), mar=c(2,2,2,2), cex=1.5)
      for(i in 1:n.model) {
            m.tr.ind <- m.pts[[i]][1] - (i-1)*n.sim
            plot.phylo(ladderize(tr.list[[i]][[m.tr.ind]]), show.tip.label=F,
                       main=paste('model =', treetypes[i]))
            axisPhylo()
      }
      mtext("Medoid Tree", side=3, line=-21.5, outer=T, font=2, cex=2)
      dev.off()
      
      save(dmat, mds, m.pts, file=paste(save.f.pre, 'mds_data.RData', sep=''))
}


# ============================================================
# L1 metric weighted
# ============================================================
print('computing weighted L1 distance same.s')
start.time <- Sys.time()

dmat.L1.w <- patch.dmat(n.model=n.model, n.sim=n.sim, save.dir=save.dir,
                        dist.method='l1', weighted=T)

end.time <- Sys.time()
print(end.time - start.time)
print(end.time)

save(dmat.L1.w, file='L1_wtd_dmat_sames.RData')
mds.analysis(dmat.L1.w, save.f.pre=paste('weighted_L1_sames_nsim', n.sim, 
                                         'ntip', n.tip, '', sep='_'))

# ============================================================
# L2 metric weighted
# ============================================================

print('computing weighted L2 distance same.s')
start.time <- Sys.time()
dmat.L2.w <- patch.dmat(n.model=n.model, n.sim=n.sim, save.dir=save.dir,
                        dist.method='l2', weighted=T)
end.time <- Sys.time()
print(end.time - start.time)
print(end.time)

save(dmat.L2.w, file='L2_wtd_dmat_sames.RData')
mds.analysis(dmat.L2.w, save.f.pre=paste('weighted_L2_sames_nsim', n.sim, 
                                         'ntip', n.tip, '', sep='_'))









