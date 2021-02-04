### Bug fixed; the bug was in rFmat and is.Fmat, not in tree_from_F.

path <- "./R/"

for(file in list.files(path)) source(file.path(path, file))

n <- 25
set.seed(1)
temp <- rFmat(n)

is.Fmat(temp)
View(temp)

debug(tree_from_F)
tree_from_F(rbind(rep(0,n),cbind(rep(0,n-1),temp)),seq(1,n-1))


