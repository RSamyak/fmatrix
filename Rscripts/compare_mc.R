
library(fmatrix)



B <- 500
R <- 50

n <- 8

set.seed(1)
start <- gen_Fmat(rcoal(n))

temp <- array(data = list(start), dim = c(2, B, R))



for(i in 1:B){
  for(j in 1:R){
    if(j == 1) {
      temp[[1, i, j]] <- temp[[2, i, j]] <- start
      current_matching <- matching(mytree_from_F(start))
      current_myencod <- my_encod(start)
      next
    }

    current_matching <- proposal_matching(current_matching)
    temp[[1, i, j]] <- gen_Fmat(fmatrix:::tree_from_matching(current_matching))
    current_myencod <- proposal_myencod(current_myencod)
    temp[[2, i, j]] <- Fmat_from_myencod(current_myencod)

  }
}

which1 <- function(u, list = F.list8){
  which(sapply(list, identical, u))
}

temp3a <- temp3b <- matrix(NA, B, R)

for(i in 1:B){
  for(j in 1:R){
    temp3a[i, j] <- which1(temp[[1, i, j]])
    temp3b[i, j] <- which1(temp[[2, i, j]])
  }
}


tempplot <- function(j){
  plot.new()
  par(mfrow = c(1, 2))
  hist(temp3a[, j], xlim = c(1, 271), breaks = 30, ylim = c(0, B/4),
       main = "Adjacent Transpositions MC",
       xlab = NULL,
       sub = sprintf("After %d steps", j))

  abline(h = B/30)
  hist(temp3b[, j], xlim = c(1, 271), breaks = 30, ylim = c(0, B/4),
       main = "Encodings MC",
       xlab = NULL,
       sub = sprintf("After %d steps", j))
  abline(h = B/30)
}

if(!dir.exists("MCpng")) dir.create("MCpng")
for(ii in 1:R){
  png(filename = file.path("MCpng", sprintf("png%d.png", ii)), width = 800, height = 400)
  tempplot(ii)
  dev.off()
}


# "convert -delay 30 *.png out.gif"
