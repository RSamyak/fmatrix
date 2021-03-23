n <- 10

library(fmatrix)

outF <- gurobi_mean(n)

M <- kingman_m(n)

outD <- fmatrix:::Dmat_from_Fmat(outF)



MD <- fmatrix:::Dmat_from_Fmat(M)


round(MD, 2)

# 1*(MD - floor(MD) == .5)
testD <- round(MD, 0)

testD - outD
outD


sum(outF)
sum(M)


# m_ij = j*(j+1)/i



