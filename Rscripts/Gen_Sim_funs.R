##This function generates a form of F matrices, I tried it for n<=10, it takes long for larger values
construct_trees<-function(n){
  k<-(n-2)*(n-1)/2
  k2<-(n-3)*(n-2)/2
  initial<-rep(1,k)
  first.column<-seq(1,k2)*(seq(1,k2)-1)/2 + 1
  first.column<-first.column[first.column<k]
  diagonal<-(seq(1,k2+1)*(seq(1,k2+1)-1)/2)[-1]
  diagonal<-diagonal[diagonal<=k]
  initial[diagonal]<-seq(1,length(first.column))
  initial<-rbind(initial,initial)
  initial[2,2]<-0
  #x_index, column, row
  F.matrixRel<-cbind(first.column,rep(1,length(first.column)),seq(1,length(first.column)))
  j<-2
  F.matrixRel<-rbind(F.matrixRel,cbind(first.column[j]+seq(1,j-1),2:j,rep(j,j-1)))
  for (j in 3:length(first.column)){ #row in F matrix
    F.matrixRel<-rbind(F.matrixRel,cbind(first.column[j]+seq(1,j-1),2:j,rep(j,j-1)))
    #column in F matrix
    for (cc in 1:(j-1)){
      xval<-F.matrixRel[F.matrixRel[,2]==cc & F.matrixRel[,3]==j,1]
      nrows.to<-nrow(initial)
      for (i in 1:nrows.to){
        if (cc==1){ #it is a first column element
          prevxval<-F.matrixRel[F.matrixRel[,2]==cc & F.matrixRel[,3]==(j-1),1]
          cond1<-max(c(0,initial[i,prevxval]-1))
          cond2<-initial[i,prevxval]
          newval<-seq(cond1,cond2)
          initial[i,xval]<-newval[1]
          if (length(newval)>1){
            for (z in 2:length(newval)){
              temp<-initial[i,]
              temp[xval]<-newval[z]
              initial<-rbind(initial,temp)
            }
          }
        }else{
          prevxval<-F.matrixRel[F.matrixRel[,2]==(cc-1) & F.matrixRel[,3]==j,1]
          prevxval2<-F.matrixRel[F.matrixRel[,2]==cc & F.matrixRel[,3]==(j-1),1]
          prevxval3<-F.matrixRel[F.matrixRel[,2]==(cc-1) & F.matrixRel[,3]==(j-1),1]
          cond1<-max(c(0, initial[i,prevxval], initial[i,prevxval2]-1, initial[i,prevxval]+initial[i,prevxval2]-initial[i,prevxval3]-1))
          cond2<-min(initial[i,prevxval]+initial[i,prevxval2]-initial[i,prevxval3],initial[i,prevxval2])
          newval<-seq(cond1,cond2)
          initial[i,xval]<-newval[1]
          if (length(newval)>1){
            for (z in 2:length(newval)){
              temp<-initial[i,]
              temp[xval]<-newval[z]
              initial<-rbind(initial,temp)
            }
          }
        }
        
      }
      
    }
  }
  return(list(F.matrixRel=F.matrixRel,res=initial)) 
}


#This functions generates a tree in phylo format to plot from an F-matrix
tree_from_F <- function(matF, coal_times){
  #generate an ape tree (Jaime's code)
  #F is the actual form used in the code that differs from the paper's notation
  
  n= dim(matF)[1]
  edge=matrix(rep(0,6*n-6),ncol=3)
  edge[,2]= (1:(2*n-2))
  vintages = c()
  times=c(rep(0,n),coal_times)
  for (j in n:2){
    new_node = 2*n-j+1
    next_leaf = intersect(which(edge[,1]==0),1:n)[1]
    F_difference = rev(matF[,j]-matF[,j-1])
    
    if (F_difference[1]==2){
      edge[next_leaf:(next_leaf+1),1]=new_node
      vintages=c(vintages, new_node)
    }
    else if (F_difference[1]==1){
      selected_vintage = which(F_difference == 2)[1]+n-1
      edge[selected_vintage,1]=new_node
      edge[next_leaf,1]=new_node
      vintages = c(vintages[vintages != selected_vintage],new_node)
    }
    else {
      selected_vintage1 =which(F_difference == 1)[1]+n-1
      selected_vintage2 =which(F_difference == 2)[1]+n-1
      edge[selected_vintage1,1]=new_node
      edge[selected_vintage2,1]=new_node
      vintages = vintages[vintages!=selected_vintage1]
      vintages = vintages[vintages!=selected_vintage2]
      vintages<-c(vintages,new_node)
    }
  }
  #edge[5:8,2]=c(6,7,8,5)
  edge[1:n,]=edge[order(edge[1:n,2]),]
  #edge=edge[order(edge[,1]),]
  
  for (j in 1:(2*n-2)) {
    #I don't understand this
    edge[j,3]=times[edge[j,1]]-times[edge[j,2]]
  }
  edge[,1]=3*n-edge[,1]
  edge[-(1:n),2]=3*n-edge[-(1:n),2]
  
  final_tree=rcoal(n,br=coal_times)
  final_tree$edge=edge[,-3]
  final_tree$edge.length=edge[,3]
  final_tree$Nnode=n-1
  class(final_tree) <- "phylo"
  final_tree <- reorder(final_tree,"postorder")
  final_tree$edge[final_tree$edge[, 2] <= n, 2] <- 1:n
  return(final_tree)
}

mytree_from_F <- function(Fmat){
  n <- ncol(Fmat) + 1
  tree_from_F(rbind(rep(0,n),cbind(rep(0,n-1),Fmat)),seq(1,n-1))
}

#This function converts the output of the construct_trees function to 
Fmat_from_construct <- function(n, output, relpos){
  # Expected syntax:
  # F.list<-list()
  # for (j in 1:nrow(res$res)){
  #   F.list[[j]]<- Fmat_from_construct(n, res$res[j,], res$F.matrixRel)
  # }
  out <- matrix(0,nrow=n-1,ncol=n-1)
  diag(out) <- seq(2,n)
  for (i in 1:length(output)){
    coord <- relpos[relpos[,1]==i,]
    out[(coord[3]+1),coord[2]] <- output[i]
  }
  return(out)
}

#Simple function to return number of cherries
num_cherries <- function(Fmat){
  n <- dim(Fmat)[1]
  c<-sum(diff(Fmat[n,])==2)
  return(c)
}

##Brute force frechet mean
brute.mean<-function(sampleF,totalF.list, dist=c("l1","l2")){
  dist <- match.arg(dist)
  dist.matrix<-matrix(0,nrow=length(totalF.list),ncol=length(sampleF))
  if(dist == "l1") for (j in 1:length(totalF.list)){
    for (i in 1:length(sampleF)){
      dist.matrix[j,i]<-sum(abs(totalF.list[[j]]-sampleF[[i]]))
    }
  }
  if(dist == "l2") for (j in 1:length(totalF.list)){
    for (i in 1:length(sampleF)){
      dist.matrix[j,i]<-sum((totalF.list[[j]]-sampleF[[i]])**2)
    }
  }
  total<-apply(dist.matrix,1,sum)
  return(which.min(total))
}

## Wrapper to plot tree given F-matrix
plotF <- function(Fmat){
  n <- ncol(Fmat) + 1
  tree<-tree_from_F(rbind(rep(0,n),cbind(rep(0,n-1),Fmat)),seq(1,n-1))
  ape::plot.phylo(ape::ladderize(tree), direction = "downwards")
}

plotF.list <- function(F.list){
  n <- length(F.list)
  t <- ceiling(sqrt(n))
  par(mfrow=c(t,t))
  for(i in 1:n){
    plotF(F.list[[i]])
  }
  # dev.off()
}

## Simple function to calculate probability of F-matrix under Kingman Model
kingman_likelihood <- function(Fmat){
  n <- ncol(Fmat) + 1
  c <- num_cherries(Fmat)
  fact <- (2^(n-1))/factorial(n-1)
  return(fact*2^(-c))
}


# Functions to calculate likelihood of F-matrix under Blum-Francois' Beta-splitting model (only beta_BF_likelihood to be exported, rest are internal)
beta_BF_expression <- function(num, denom){
  k <- length(num)
  if(!length(denom)==k) stop("same length needed")
  
  t <- function(x){
    out <- 1
    for(ind in 1:k){
      out <- out * (x+num[ind])/(2*x + denom[ind])
    }
    return(out)
  }
  return(t)
}

find_drop <- function(upperrow, lowerrow){
  check <- lowerrow == (upperrow - 1)
  out <- which(check)[1]
  return(out + 1)
}

my_encod <- function(Fmat){
  n <- ncol(Fmat) + 1
  nn <- n - 1
  
  out <- c()
  
  for(ind1 in 1:(nn-1)){
    out <- append(out, find_drop(Fmat[ind1, ], Fmat[ind1 + 1, ]))
  }
  return(out)
}

child_seq <- function(myencod, node){
  nn <- length(myencod) + 1
  
  left <- 2
  
  childL <- c()
  childR <- c()
  
  for (ind1 in (node-1):(nn-1)){
    
    if(myencod[ind1] == node){
      if(left == 2) childR <- append(childR, ind1 + 2)
      if(left == 1) childL <- append(childL, ind1 + 2)
      
      left <- left - 1
    } 
    else if(myencod[ind1] %in% childL){
      childL <- append(childL, ind1 + 2)
    } else if(myencod[ind1] %in% childR){
      childR <- append(childR, ind1 + 2)
    }
  }
  return(list(childL=childL
              , childR=childR
  ))
}


beta_BF_numden_from_encod <- function(myencod, beta = 0){
  
  n <- length(myencod) + 2
  
  num <- c()
  den <- c()
  
  for(ind1 in 4:n){ # ind1 is size of tree
    
    for(ind2 in (ind1-1):2 ){ # ind2 looks for ancestry
      
      childs <- child_seq(myencod[1:(ind1-2)], ind2)
      childL <- childs$childL
      childR <- childs$childR
      
      if(ind1 %in% childL){
        l <- length(childL) 
        r <- length(childR) + 1
        if(!(l==1 & r == 1)){
          num <- append(num, l)
          den <- append(den, l+r)
        }
      } else if(ind1 %in% childR){
        l <- length(childL) + 1
        r <- length(childR) 
        if(!(l==1 & r == 1)){
          num <- append(num, r)
          den <- append(den, l+r)
        }
      }
      
    }
    
  }
  return(list(num=num,den=den))
}


#### Note this is Blum and Francois (2006) Beta-splitting model, not Aldous; as described in Sainudiin, Veber (2016)
#### Only export this function, the rest are called
beta_BF_likelihood <- function(Fmat, betas=0){
  n <- ncol(Fmat) + 1
  
  if(n <= 3) return(rep(1,length(betas)))
  
  numden <- beta_BF_numden_from_encod(my_encod(Fmat))
  num <- numden$num
  den <- numden$den
  
  func <- beta_BF_expression(num, den)
  return(sapply(betas, FUN = func))
}



# # Aldous Beta-splitting model - no explicit likelihood, what to do?
# beta_aldous_likelihood <- function(Fmat, betas=0){
#   n <- ncol(Fmat) + 1
#   
#   if(n <= 3) return(rep(1,length(betas)))
#   
#   
# }


## Naive weighted mean function for F-matrices, output need not be F-matrix, but rounded F-matrix is close?
meanF <- function(F.list, weights=rep(1,length(F.list))){
  n <- ncol(F.list[[1]]) + 1
  
  out <- matrix(0, n-1, n-1)
  for (i in 1:length(F.list)){
    out <- out + weights[i]*F.list[[i]]
  }
  return(out/sum(weights))
}

## Function to test if a given matrix is a valid F-matrix, returns a negative value if not
is.Fmat.debug <- function(Fmat, tol=0){
  n <- ncol(Fmat) + 1
  if(! nrow(Fmat) == (n-1)) return(-0.1)
  if(! all( abs(Fmat[upper.tri(Fmat)] - 0) <= tol)) return(-2)
  if(! all( abs(diag(Fmat) - 2:n) <= tol)) return(-3)
  
  for(i1 in 1:(n-2)){
    if(!abs(Fmat[i1+1,i1] - i1) <= tol) return(-4)
  }
  
  if(n >=4 ){
    for(i1 in 3:(n-1)){
      if(!all(any(abs(Fmat[i1, 1] - Fmat[i1-1,1]) <= tol,
                  abs(Fmat[i1, 1] - Fmat[i1-1,1] + 1) <= tol
      ),
      Fmat[i1, 1] >= -tol
      )) return(-5)
      
      if(i1 >=4){
        for(i2 in 2:(i1 - 2)){
          if(!all(any(abs(Fmat[i1, i2] - Fmat[i1-1,i2]) <= tol,
                      abs(Fmat[i1, i2] - Fmat[i1-1,i2] + 1) <= tol
          ),
          Fmat[i1-1, i2] - Fmat[i1, i2] >= Fmat[i1-1, i2-1] - Fmat[i1, i2-1] - tol
          )) return(-6)
        }
      }
    }
  }
  
  return(T)
}

is.Fmat <- function(Fmat){
  return(is.Fmat.debug(Fmat)>0)
}

## Closest Fmat, given full F.list
closest_Fmat <- function(qmat, F.list){
  if(is.Fmat(qmat)) return(qmat) ## Note this is assuming qmat is in F.list
  return(F.list[[brute.mean(list(qmat),F.list)]])
}

where.Fmat.debug <- function(Fmat, tol=0){
  n <- ncol(Fmat) + 1
  if(! nrow(Fmat) == (n-1)) return(-0.1)
  if(! all( abs(Fmat[upper.tri(Fmat)] - 0) <= tol)) return(-2)
  if(! all( abs(diag(Fmat) - 2:n) <= tol)) return(-3)
  
  for(i1 in 1:(n-2)){
    if(!abs(Fmat[i1+1,i1] - i1) <= tol) return(-4)
  }
  
  if(n >= 4){
    
    flag <- matrix(NA, n-1, n-1)
    
    for(i1 in 3:(n-1)){
      if(!all(any(abs(Fmat[i1, 1] - Fmat[i1-1,1]) <= tol,
                  abs(Fmat[i1, 1] - Fmat[i1-1,1] + 1) <= tol
      ),
      Fmat[i1, 1] >= -tol
      )) {
        flag[i1, 1] <- F
      } else flag[i1, 1] <- T
      
      if(i1 >=4){
        for(i2 in 2:(i1 - 2)){
          if(!all(any(abs(Fmat[i1, i2] - Fmat[i1-1,i2]) <= tol,
                      abs(Fmat[i1, i2] - Fmat[i1-1,i2] + 1) <= tol
          ),
          Fmat[i1-1, i2] - Fmat[i1, i2] >= Fmat[i1-1, i2-1] - Fmat[i1, i2-1] - tol
          )) {
            flag[i1, i2] <- F
          } else flag[i1, i2] <- T
        }
      }
    }
  }
  
  return(flag)
}



## Candidate Fmat
nearby_Fmat <- function(qmat, rd = .25, tol = 1e-5, ret.flag = F){
  if(is.Fmat(qmat)) return(qmat)
  n <- ncol(qmat) + 1
  
  if(! nrow(qmat) == (n-1)) return(-0.1)
  
  Fmat <- matrix(NA, n-1, n-1)
  flag <- matrix(NA, n-1, n-1)
  
  if(! all(abs(qmat[upper.tri(qmat)]-0) < tol)) return(-2)
  
  Fmat[upper.tri(Fmat)] <- 0
  
  if(! all(abs(diag(qmat) - 2:n) < tol) ) return(-3)
  
  diag(Fmat) <- 2:n
  
  for(i1 in 1:(n-2)){
    if(abs(qmat[i1+1,i1] - i1) > tol) return(-4)
    Fmat[i1+1, i1] <- i1
  }
  
  if(n >= 4){
    
    
    
    for(i1 in 3:(n-1)){
      if(!all(qmat[i1, 1] %in% c(Fmat[i1-1,1], Fmat[i1-1,1] - 1),
              qmat[i1, 1] >= 0
      )) {
        
        d.same <- abs(qmat[i1,1] - Fmat[i1-1,1])
        d.drop <- abs(qmat[i1,1] - Fmat[i1-1,1] + 1)
        
        if(d.same < tol){
          Fmat[i1, 1] <- Fmat[i1-1, 1]
          flag[i1, 1] <- 1
        }
        else if(d.drop < tol){
          Fmat[i1, 1] <- Fmat[i1-1, 1] - 1
          flag[i1, 1] <- 1
        }
        else if(d.same < rd){
          Fmat[i1, 1] <- Fmat[i1-1, 1]
          flag[i1, 1] <- 0
        }
        else if(d.drop < rd){
          Fmat[i1, 1] <- Fmat[i1-1, 1] - 1
          flag[i1, 1] <- 0
        } 
        else if(d.same < d.drop) {
          Fmat[i1, 1] <- Fmat[i1-1, 1]
          flag[i1, 1] <- -1
        } 
        else {
          Fmat[i1, 1] <- Fmat[i1-1, 1] - 1
          flag[i1, 1] <- -1
        }
        
      } else {
        flag[i1, 1] <- 1
        Fmat[i1, 1] <- qmat[i1, 1]
      }
      
      if(i1 >=4){
        for(i2 in 2:(i1 - 2)){
          if(!all(qmat[i1, i2] %in% c(qmat[i1-1, i2], qmat[i1-1, i2] -1),
                  qmat[i1-1, i2] - qmat[i1, i2] >= qmat[i1-1, i2-1] - qmat[i1, i2-1]
          )) {
            
            d.same <- abs(qmat[i1,i2] - Fmat[i1-1,i2])
            d.drop <- abs(qmat[i1,1] - Fmat[i1-1,i2] + 1)
            
            if(abs(Fmat[i1-1, i2-1] - Fmat[i1, i2-1] - 1) < tol){
              Fmat[i1, i2] <- Fmat[i1-1, i2] - 1
              flag[i1, i2] <- flag[i1, i2-1]
            }
            
            else if(d.same < tol){
              Fmat[i1, i2] <- Fmat[i1-1, i2]
              flag[i1, i2] <- 1
            }
            else if(d.drop < tol){
              Fmat[i1, i2] <- Fmat[i1-1, i2] - 1
              flag[i1, i2] <- 1
            }
            else if(d.same < rd){
              Fmat[i1, i2] <- Fmat[i1-1, i2]
              flag[i1, i2] <- 0
            }
            else if(d.drop < rd){
              Fmat[i1, i2] <- Fmat[i1-1, i2] - 1
              flag[i1, i2] <- 0
            } 
            else if(d.same < d.drop) {
              Fmat[i1, i2] <- Fmat[i1-1, i2]
              flag[i1, i2] <- -1
            } 
            else {
              Fmat[i1, i2] <- Fmat[i1-1, i2] - 1
              flag[i1, i2] <- -1
            }
            
          } else {
            flag[i1, i2] <- 1
            Fmat[i1, i2] <- qmat[i1, i2]
          }
        }
      }
    }
  }
  
  if(ret.flag){
    return(list(Fmat, flag))
  } else return(Fmat)
}


nearby_Fmats <- function(qmat, rd = .25, tol = 1e-5, MAX.LENGTH = 300){
  if(is.Fmat(qmat)) return(qmat)
  
  n <- ncol(qmat) + 1
  
  if(! nrow(qmat) == (n-1)) return(-0.1)
  
  Fmat <- matrix(NA, n-1, n-1)
  
  if(! all(abs(qmat[upper.tri(qmat)]-0) < tol)) return(-2)
  
  Fmat[upper.tri(Fmat)] <- 0
  
  if(! all(abs(diag(qmat) - 2:n) < tol) ) return(-3)
  
  diag(Fmat) <- 2:n
  
  for(i1 in 1:(n-2)){
    if(abs(qmat[i1+1,i1] - i1) > tol) return(-4)
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
            temp[i1, 1] <- temp[i1-1, 1] - 1
          }
          else if(d.same < rd){
            temp[i1, 1] <- temp[i1-1, 1]
          }
          else if(d.drop < rd){
            temp[i1, 1] <- temp[i1-1, 1] - 1
          } 
          else {
            temp[i1, 1] <- temp[i1-1, 1]
            if(temp[i1-1, 1] - 1 >= 0){
              tempnew <- temp
              tempnew[i1, 1] <- temp[i1-1, 1] - 1
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
                    temp[i1-1, i2] - qmat[i1, i2] >= temp[i1-1, i2-1] - temp[i1, i2-1]
            )) {
              
              d.same <- abs(qmat[i1,i2] - temp[i1-1,i2])
              d.drop <- abs(qmat[i1,i2] - temp[i1-1,i2] + 1)
              
              last.drop <- abs(temp[i1-1, i2-1] - temp[i1, i2-1] - 1)
              
              
              if(last.drop < tol){
                temp[i1, i2] <- temp[i1-1, i2] - 1
              }
              
              else if(d.same < tol){
                temp[i1, i2] <- temp[i1-1, i2]
              }
              else if(d.drop < tol){
                temp[i1, i2] <- temp[i1-1, i2] - 1
              }
              else if(d.same < rd){
                temp[i1, i2] <- temp[i1-1, i2]
              }
              else if(d.drop < rd){
                temp[i1, i2] <- temp[i1-1, i2] - 1
              } 
              else {
                
                temp[i1, i2] <- temp[i1-1, i2]
                
                if(temp[i1-1, i2] - 1 >= 0){
                  tempnew <- temp
                  tempnew[i1, i2] <- temp[i1-1, i2] - 1
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


#### Fmat - Dmat and back
Fmat_from_Dmat <- function(Dmat){
  n <- ncol(Dmat) + 1
  
  Fmat <- matrix(0, n-1, n-1)
  Fmat[, 1] <- Dmat[, 1]
  
  if(n>=2) {
    for(i in 2:(n-1)){
      if(i == n-1) {
        Fmat[n-1, n-1] <- sum(Dmat[n-1, 1:(n-1)])
        break
      }
      Fmat[i:(n-1), i] <- apply(Dmat[i:(n-1), 1:i], 1, sum)
    }
    
  }
  return(Fmat)
}

Dmat_from_Fmat <- function(Fmat){
  n <- ncol(Fmat) + 1
  
  Dmat <- matrix(0, n-1, n-1)
  Dmat[, 1] <- Fmat[, 1]
  
  if(n>=2) for(i in 2:(n-1)){
    Dmat[i:(n-1), i] <- Fmat[i:(n-1), i] - Fmat[i:(n-1), i-1]
  }
  return(Dmat)
}

#### gen_Fmat for forgetting edge lengths
forget_lengths <- function(tree){
  m <- max(tree$edge)
  n <- m - tree$Nnode
  
  o <- which(tree$edge[,2] <= n)
  tree$edge.length[o] <- m+1-tree$edge[o,1]
  tree$edge.length[-o] <- tree$edge[-o,2] - tree$edge[-o,1]
  
  return(tree)
}

gen_Fmat_unweighted <- function(tree){
  gen_Fmat(forget_lengths(tree))
}


#### Converting phylo tree to and from perfect matchings
matching<-function(tree){
  o <- order(tree$edge[,1], decreasing = T)
  return(tree$edge[o, 2])
}

tree_from_matching <- function(config){
  
  n<-length(config)/2 + 1
  
  tree<-rcoal(n,br=1)
  
  parents<-rep(seq(n+1,2*n-1),each=2)
  
  tree$edge[,1]<-rev(parents)
  tree$edge[,2]<-config
  
  m <- max(tree$edge)
  n <- m - tree$Nnode
  
  o <- which(tree$edge[,2] <= n)
  tree$edge.length[o] <- m+1-tree$edge[o,1]
  tree$edge.length[-o] <- tree$edge[-o,2] - tree$edge[-o,1]
  return(tree)
}

proposal <- function(config){
  n <- length(config)/2 + 1
  
  ret <- config
  
  i <- sample(seq(1,n-2),1) # choose a pair to change
  
  selected <- config[2*i + 1:2]
  
  eligible <- selected[selected > (2*n - i) | selected <= n] # either nodes which are lower or leaves
  
  if(length(eligible)==1) sampleR <- eligible
  else sampleR <- sample(eligible, size=1)
  
  sampleL <- sample(config[2*i - 1:0], size=1)

  place1 <- which(ret == sampleL)
  place2 <- which(ret == sampleR)
  
  ret[place1] <- sampleR
  ret[place2] <- sampleL
  
  return(ret)
}

#### Better storage for distance than list?
matrix_list <- function(F.list){
  sapply(F.list, function(u){
   as.vector(u[lower.tri(u)]) 
  })
}

reduce_to_vector <- function(Fmat){
  as.vector(Fmat[lower.tri(Fmat)])
}
