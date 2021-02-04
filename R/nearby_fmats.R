#' @param qmat
#' @param rd
#' @param tol
#' @param MAX.LENGTH
#'
#' @export
nearby_Fmats <- function(qmat, rd = .25, tol = 1e-5, MAX.LENGTH = 300){
  if(is.Fmat(qmat)) return(qmat)

  n <- ncol(qmat) + 1

  if(! nrow(qmat) == (n-1)) stop("Number of rows and columns do not match")

  Fmat <- matrix(NA, n-1, n-1)

  if(! all(abs(qmat[upper.tri(qmat)]-0) < tol)) stop("Upper triangular part of qmat is not zero")

  Fmat[upper.tri(Fmat)] <- 0

  if(! all(abs(diag(qmat) - 2:n) < tol) ) stop("Diagonal of qmat improper")

  diag(Fmat) <- 2:n

  for(i1 in 1:(n-2)){
    if(abs(qmat[i1+1,i1] - i1) > tol) stop("Subdiagonal of qmat improper")
    Fmat[i1+1, i1] <- i1
  }

  ret.list <- list()
  ind.list <- list()
  ret.list <- append(ret.list, list(Fmat))
  ind.list <- append(ind.list, list(c(3,1)))

  if(n >= 4){

    ilist <- 0

    while(ilist < length(ret.list)){

      ilist <- ilist + 1
      temp <- ret.list[[ilist]]

      skiptill <- ind.list[[ilist]]

      for(i1 in 3:(n-1)){
        if(i1 < skiptill[1]) next

        if(i1 == skiptill[1] & skiptill[2] > 1){}
        else if(!all(qmat[i1, 1] %in% c(temp[i1-1,1], temp[i1-1,1] - 1),
                     qmat[i1, 1] >= 0
        )) {

          d.same <- abs(qmat[i1,1] - temp[i1-1,1])
          d.drop <- abs(qmat[i1,1] - temp[i1-1,1] + 1)

          if(d.same < tol){
            temp[i1, 1] <- temp[i1-1, 1]
          }
          else if(d.drop < tol){
            temp[i1, 1] <- max(temp[i1-1, 1] - 1, 0)
          }
          else if(d.same < rd){
            temp[i1, 1] <- temp[i1-1, 1]
          }
          else if(d.drop < rd){
            temp[i1, 1] <- max(temp[i1-1, 1] - 1, 0)
          }
          else {
            temp[i1, 1] <- temp[i1-1, 1]
            if(temp[i1-1, 1] - 1 >= 0){
              tempnew <- temp
              tempnew[i1, 1] <- max(temp[i1-1, 1] - 1, 0)
              if(length(ret.list) < MAX.LENGTH){
                ret.list <- append(ret.list, list(tempnew))
                ind.list <- append(ind.list, list(c(i1, 1+1)))
              }
            }

          }

        } else {
          temp[i1, 1] <- qmat[i1, 1]
        }

        if(i1 >= 4){
          for(i2 in 2:(i1 - 2)){

            if(i1 == skiptill[1] & i2 < skiptill[2]) next

            if(!all(qmat[i1, i2] %in% c(temp[i1-1, i2], temp[i1-1, i2] - 1),
                    qmat[i1, i2] >= temp[i1-1, i2],
                    temp[i1-1, i2] - qmat[i1, i2] >= temp[i1-1, i2-1] - temp[i1, i2-1]
            )) {

              d.same <- abs(qmat[i1,i2] - temp[i1-1,i2])
              d.drop <- abs(qmat[i1,i2] - temp[i1-1,i2] + 1)

              if(temp[i1-1, i2] - 1 < temp[i1, i2-1]) d.drop <- Inf

              last.drop <- abs(temp[i1-1, i2-1] - temp[i1, i2-1] - 1)


              if(last.drop < tol){
                temp[i1, i2] <- max(temp[i1-1, i2] - 1, 0)
              }

              else if(d.same < tol){
                temp[i1, i2] <- temp[i1-1, i2]
              }
              else if(d.drop < tol){
                temp[i1, i2] <- max(temp[i1-1, i2] - 1, 0)
              }
              else if(d.same < rd){
                temp[i1, i2] <- temp[i1-1, i2]
              }
              else if(d.drop < rd){
                temp[i1, i2] <- max(temp[i1-1, i2] - 1, 0)
              }
              else {

                temp[i1, i2] <- temp[i1-1, i2]

                if(temp[i1-1, i2] - 1 >= max(temp[i1, i2-1], 0)){
                  tempnew <- temp
                  tempnew[i1, i2] <- max(temp[i1-1, i2] - 1, 0)
                  if(length(ret.list) < MAX.LENGTH){
                    ret.list <- append(ret.list, list(tempnew))
                    ind.list <- append(ind.list, list(c(i1, i2+1)))
                  }
                }
              }

            } else {
              temp[i1, i2] <- qmat[i1, i2]
            }
          }
        }
      }

      ret.list[[ilist]] <- temp

    }

  }

  return(ret.list)
}

# ## Candidate Fmat
# nearby_Fmat <- function(qmat, rd = .25, tol = 1e-5, ret.flag = F){
#   if(is.Fmat(qmat)) return(qmat)
#   n <- ncol(qmat) + 1
#
#   if(! nrow(qmat) == (n-1)) stop("Number of rows and columns do not match")
#
#   Fmat <- matrix(NA, n-1, n-1)
#   flag <- matrix(NA, n-1, n-1)
#
#   if(! all(abs(qmat[upper.tri(qmat)]-0) < tol)) stop("Upper triangular part of qmat is not zero")
#
#   Fmat[upper.tri(Fmat)] <- 0
#
#   if(! all(abs(diag(qmat) - 2:n) < tol) ) stop("Diagonal of qmat improper")
#
#   diag(Fmat) <- 2:n
#
#   for(i1 in 1:(n-2)){
#     if(abs(qmat[i1+1,i1] - i1) > tol) stop("Subdiagonal of qmat improper")
#     Fmat[i1+1, i1] <- i1
#   }
#
#   if(n >= 4){
#
#
#
#     for(i1 in 3:(n-1)){
#       if(!all(qmat[i1, 1] %in% c(Fmat[i1-1,1], Fmat[i1-1,1] - 1),
#               qmat[i1, 1] >= 0
#       )) {
#
#         d.same <- abs(qmat[i1,1] - Fmat[i1-1,1])
#         d.drop <- abs(qmat[i1,1] - Fmat[i1-1,1] + 1)
#
#         if(d.same < tol){
#           Fmat[i1, 1] <- Fmat[i1-1, 1]
#           flag[i1, 1] <- 1
#         }
#         else if(d.drop < tol){
#           Fmat[i1, 1] <- max(Fmat[i1-1, 1] - 1, 0)
#           flag[i1, 1] <- 1
#         }
#         else if(d.same < rd){
#           Fmat[i1, 1] <- Fmat[i1-1, 1]
#           flag[i1, 1] <- 0
#         }
#         else if(d.drop < rd){
#           Fmat[i1, 1] <- max(Fmat[i1-1, 1] - 1, 0)
#           flag[i1, 1] <- 0
#         }
#         else if(d.same < d.drop) {
#           Fmat[i1, 1] <- Fmat[i1-1, 1]
#           flag[i1, 1] <- -1
#         }
#         else {
#           Fmat[i1, 1] <- max(Fmat[i1-1, 1] - 1, 0)
#           flag[i1, 1] <- -1
#         }
#
#       } else {
#         flag[i1, 1] <- 1
#         Fmat[i1, 1] <- qmat[i1, 1]
#       }
#
#       if(i1 >=4){
#         for(i2 in 2:(i1 - 2)){
#           if(!all(qmat[i1, i2] %in% c(qmat[i1-1, i2], qmat[i1-1, i2] -1),
#                   qmat[i1-1, i2] - qmat[i1, i2] >= qmat[i1-1, i2-1] - qmat[i1, i2-1]
#           )) {
#
#             d.same <- abs(qmat[i1,i2] - Fmat[i1-1,i2])
#             d.drop <- abs(qmat[i1,1] - Fmat[i1-1,i2] + 1)
#
#             if(abs(Fmat[i1-1, i2-1] - Fmat[i1, i2-1] - 1) < tol){
#               Fmat[i1, i2] <- max(Fmat[i1-1, i2] - 1, 0)
#               flag[i1, i2] <- flag[i1, i2-1]
#             }
#
#             else if(d.same < tol){
#               Fmat[i1, i2] <- Fmat[i1-1, i2]
#               flag[i1, i2] <- 1
#             }
#             else if(d.drop < tol){
#               Fmat[i1, i2] <- max(Fmat[i1-1, i2] - 1, 0)
#               flag[i1, i2] <- 1
#             }
#             else if(d.same < rd){
#               Fmat[i1, i2] <- Fmat[i1-1, i2]
#               flag[i1, i2] <- 0
#             }
#             else if(d.drop < rd){
#               Fmat[i1, i2] <- max(Fmat[i1-1, i2] - 1, 0)
#               flag[i1, i2] <- 0
#             }
#             else if(d.same < d.drop) {
#               Fmat[i1, i2] <- Fmat[i1-1, i2]
#               flag[i1, i2] <- -1
#             }
#             else {
#               Fmat[i1, i2] <- max(Fmat[i1-1, i2] - 1, 0)
#               flag[i1, i2] <- -1
#             }
#
#           } else {
#             flag[i1, i2] <- 1
#             Fmat[i1, i2] <- qmat[i1, i2]
#           }
#         }
#       }
#     }
#   }
#
#   if(ret.flag){
#     return(list(Fmat, flag))
#   } else return(Fmat)
# }
