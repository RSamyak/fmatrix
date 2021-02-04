##This function generates a form of F matrices, I tried it for n<=10, it takes long for larger values
#you need library(ape)
library(ape)
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

####Example
n<-6
res<-construct_trees(n)
nrow(res$res)

#Indeed for n=5, there are 5 trees
#we construct the Fmatrices

F.list<-list()
for (j in 1:nrow(res$res)){
  F.list[[j]]<-matrix(0,nrow=n-1,ncol=n-1)
  diag(F.list[[j]])<-seq(2,n)
  for (i in 1:ncol(res$res)){
    coord<-res$F.matrixRel[res$F.matrixRel[,1]==i,]
    F.list[[j]][(coord[3]+1),coord[2]]<-res$res[j,i]
  }
}

Fmat<-F.list

#FmatProb<-rep(0,length(Fmat))
#fact<-(2^(n-1))/factorial(n-1)
#for (j in 1:length(FmatProb)){
#  c<-sum(diff(Fmat[[j]][9,])==2)
#  FmatProb[j]<-fact*2^(-c)
  
#}

#What is the maximum distance?
n.tr <- length(Fmat)
Dmat.l1 <- matrix(NA, nrow=n.tr, ncol=n.tr)
diag(Dmat.l1) <- 0
comp.ind <- t(combn(n.tr, 2))

tmp.dist <- rep(NA, dim(comp.ind)[1])
for (i in 1:dim(comp.ind)[1]) {
  tmp.ind <- comp.ind[i,]
  tmp.dist[i] <- sum(abs(Fmat[[tmp.ind[1]]] - Fmat[[tmp.ind[2]]]))
}

Dmat.l1[comp.ind] <- tmp.dist
tmp <- Dmat.l1[upper.tri(Dmat.l1)]
Dmat.l1 <- t(Dmat.l1)
Dmat.l1[upper.tri(Dmat.l1)] <- tmp

#wDmat.l1<-Dmat.l1*matrix(FmatProb,nrow=nrow(Dmat.l1),ncol=ncol(Dmat.l1),byrow=TRUE)
max(Dmat.l1)
mean<-apply(Dmat.l1,1,sum)
#meanw<-apply(wDmat.l1,1,sum)
plot(table(Dmat.l1[Dmat.l1>0])/2,ylab="frequency",main="Histogram of pairwise d1, n=10")
for (j in 1:ncol(Dmat.l1)){
  if (sum(Dmat.l1[,j]==max(Dmat.l1))>0){print(j)}
}

for (j in 1:length(mean)){
  if (mean[j]==min(mean)){print(j)}
}

#for (j in 1:length(mean)){
 # if (meanw[j]==min(meanw)){print(j)}
#}

##Distance to the mean
dist_to_mean<-rep(0,length(Fmat))
for (j in 1:length(Fmat)){
dist_to_mean[j]<-sum(abs(Fmat[[j]] - Fmat[[8]]))
}

#Kingman distance to the mean
#dist_to_meanK<-rep(0,length(Fmat))
#for (j in 1:length(Fmat)){
#  dist_to_meanK[j]<-sum(abs(Fmat[[j]] - Fmat[[2017]]))
#}

#dist_to_meanK<-cbind(dist_to_meanK,FmatProb)

#tableW<-aggregate(dist_to_meanK[,2],by=list(Distance=dist_to_meanK[,1]),FUN=sum)

plot(table(dist_to_mean),ylab="frequency",main="Histogram of d1 to the mean, n=10",xlab="d(T,Tmean)")
plot(density(dist_to_mean),xlab="distance to mean", ylab="Probability",main="Uniform vs Kingman",ylim=c(0,.15))
points(tableW,type="l",col="red")

tree<-tree_from_F(rbind(rep(0,n),cbind(rep(0,n-1),F.list[[8]])),seq(1,n-1))
plot(ladderize(tree),direction="downwards",show.tip.label = FALSE,edge.width =2)
n_trees<-nrow(res$res)
par(mfrow=c(ceiling(n_trees/3),3))
for (j in 1:n_trees){
tree<-tree_from_F(rbind(rep(0,n),cbind(rep(0,n-1),F.list[[j]])),seq(1,n-1))
plot(ladderize(tree),direction="downwards")
}

par(mfrow=c(1,2))
j<-2
tree<-tree_from_F(rbind(rep(0,n),cbind(rep(0,n-1),F.list[[j]])),seq(1,n-1))
plot(ladderize(tree),direction="downwards")
j<-7936
tree<-tree_from_F(rbind(rep(0,n),cbind(rep(0,n-1),F.list[[j]])),seq(1,n-1))
plot(ladderize(tree),direction="downwards")

##For playing with Kingman library(ape)
treeKingman<-rcoal(5)
plot(treeKingman,direction="downwards")

#####################################################################
##Brute force frechet mean
brute.mean<-function(sampleF,totalF.list){
  dist.matrix<-matrix(0,nrow=length(totalF.list),ncol=length(sampleF))
  for (j in 1:length(totalF.list)){
    for (i in 1:length(sampleF)){
      dist.matrix[j,i]<-sum(abs(totalF.list[[j]]-sampleF[[i]]))
    }
  }
  total<-apply(dist.matrix,1,sum)
  return(seq(1,length(totalF.list))[total==min(total)])
}

chosen<-sort(sample(seq(1,length(F.list)),5))
sampleF<-F.list[chosen]
totalF.list<-F.list
totalF.list[[brute.mean(sampleF,totalF.list)]]

list1<-F.list
keep<-rep(0,length(list1))
for (j in 1:length(list1)){
  if (list1[[j]][7,5]==4){
    keep[j]<-1
  }
}
list1<-list1[keep==1]
list2<-F.list[keep==0]
matdist<-matrix(0,nrow=length(list1),ncol=length(list2))
for (j in 1:length(list1)){
  for (i in 1:length(list2)){
    matdist[j,i]<-sum(abs(list1[[j]]-list2[[i]]))
  }
}
max(apply(matdist,1,max))

Frechet_mean<-function(n,sampleF){
  #work in progress, this is the idea of the naive algorithm
  k<-(n-2)*(n-1)/2
  k2<-(n-3)*(n-2)/2
  initial<-rep(1,k)
  first.column<-seq(1,k2)*(seq(1,k2)-1)/2 + 1
  first.column<-first.column[first.column<k]
  diagonal<-(seq(1,k2+1)*(seq(1,k2+1)-1)/2)[-1]
  diagonal<-diagonal[diagonal<=k]
  initial[diagonal]<-seq(1,length(first.column))
  val<-0
  for (j in 1:length(sampleF)){
    val<-val+sampleF[[j]][3,1]
  }
  initial[2]<-round(val/length(sampleF),0)
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






for (j in 1:5){
  print(sampleF[[j]][5,2])
}

potential<-totalF.list


potential2<-potential
keep<-rep(0,length(potential))

for (j in 1:length(potential)){
  if (potential[[j]][3,1]==1 & potential[[j]][4,1]==0 & potential[[j]][4,2]==1 & potential[[j]][5,1]==0 & potential[[j]][5,2]==1& potential[[j]][5,3]==3){
    keep[j]<-1
  }
}
  potential<-potential[keep==1]
  length(potential)

######################################################################
#Figure in 3D
s3d<-scatterplot3d(F21=c(0,0,1,1,1),F31=c(0,0,0,1,1),F32=c(1,2,1,1,2))
######################################################################
###Algorithm for maximum d1 distance
maxdistance<-function(n){
  if (n<4){res<-0}
  if (n==4){res<-1}
  if (n>4){
    res<-1
    for (i in 5:n){
      if (i%%2==0){
        res<-res+0.25*((i-2)^{2})  
      }else{
        res<-res+0.25*(i^2 -4*i+3)
      }
      
    }
  }
  return(res)
}

maxdistance(10)

maxd<-rep(0,10000)
for (j in 4:10003){
  maxd[j-3]<-maxdistance(j)
}
par(mfrow=c(1,2))

plot(seq(4,33),maxd[1:30],type="l",xlab="n", ylab="Maximum d1 distance")
points(seq(4,33),((1/12)*seq(4,33)^3),col="red",type="l")

plot(seq(4,10003),maxd,type="l",xlab="n", ylab="")
points(seq(4,10003),(1/12)*seq(4,10003)^3,col="red",type="l")


##let's say you have a kingman tree
tree<-rcoal(10)
plot(tree,direction="downwards")
nodelabels()
tiplabels()
tree$edge #will give you the matchings.

##Julia: Given a function N(t)= with population size, simulate a coal.time
coalsim(0,2,traj=)

