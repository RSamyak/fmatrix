
set.seed(1)

Fmat <- rFmat(8)

plotF(Fmat, node.labels = list(text = c("(1)", "(2)" , "(3)", "(3)", "(4)", "(5)", "(7)"), bg = "yellow"))

Fmat

fmatrix:::my_encod(Fmat)

View(proposal_matching)


# for(i in 2:length(temp)){
#   if(i%%2 == 1 & temp[i] < temp[i-1]) {cat(i); break}
#   if(i%%2 == 0 & temp[i] > temp[i-1]) {cat(i); break}
# }
# cat("done")


str <- c(2, 3, 3, 2, 4, 5, 6, 7)

str2 <- c(2, 3, 3, 2, 5, 5, 6, 7)

temp <- Fmat_from_myencod(str)
is.Fmat(temp)

temp2 <- Fmat_from_myencod(str2)
is.Fmat(temp2)

par(mfrow = c(1,2))
plotF(temp, node.labels = T)
plotF(temp2, node.labels = T)

encod.list <- lapply(F.list5, my_encod)
plotF.list(F.list5, node.labels = T)
View(plotF.list)
