library("ape")
library("devtools")
#If you need to install phylodyn:
#install_github("JuliaPalacios/phylodyn")
library("phylodyn")

n<-5

rtree<-rcoal(n,br=rep(1,n-1))
#you can just plot it
plot(rtree,direction="downwards")
nodelabels()
tiplabels()

#This function takes a tree as input and creates a vector of perfect matchings
matching<-function(tree){
  o <- order(tree$edge[,1], decreasing = T)
  return(tree$edge[o, 2])
}

pairs<-matching(rtree)

####Markov Chain
Niter<-50

config<-matrix(0,nrow=Niter,ncol=length(pairs))

config[1,]<-pairs

for (j in 2:Niter){
  pair<-sample(seq(1,n-2),1)
  print(pair)
  eligible<-config[j-1,c(2*(pair+1)-1,2*(pair+1))] [config[j-1,c(2*(pair+1)-1,2*(pair+1))]<(pair+n)]
  if (length(eligible)==1){sampleR<-eligible}else{sampleR<-sample(eligible,size=1)}
  sampleR
  sampleL<-sample(config[j-1,c(2*(pair)-1,2*(pair))],size=1)
  sampleL
  config[j,]<-config[j-1,]
  val1<-seq(1,ncol(config))[config[j,]==sampleL]
  val2<-seq(1,ncol(config))[config[j,]==sampleR]
  config[j,val1]<-sampleR
  config[j,val2]<-sampleL
}



#this functions takes a perfect matching (config1) and returns a tree with the same characteristics of tree
construct_tree<-function(config){
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

trees_chain<-list()
trees_chain[[1]]<-rtree
for (j in 2:Niter){
  trees_chain[[j]]<-construct_tree(config[j,])
}

#check the first 4 trees of the chain
pdf("abc")
par(mfrow=c(5,5))
for(j in 1:25){
plot(trees_chain[[j]],direction="downwards")
}

#We can now count how many times we visit each "ranked tree shape"


Fmat_list <- list()
for (i in 1:length(trees_chain)) {
  Fmat_list[[i]] <- gen_Fmat(trees_chain[[i]])
}
print(Fmat_list)

uniqueList<-unique(Fmat_list)
frequencyF<-rep(0,length(uniqueList))
for (j in 1:length(frequencyF)){
  for (k in 1:length(Fmat_list)){
    frequencyF[j]<-frequencyF[j]+ifelse(sum(abs(Fmat_list[[k]]-uniqueList[[j]]))==0,1,0)
  }
}

