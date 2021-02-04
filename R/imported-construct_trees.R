#' @param n
#'
#'
#' @export
construct_trees<-function(n){
  ##This function generates a form of F matrices, I tried it for n<=10, it takes long for larger values

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
