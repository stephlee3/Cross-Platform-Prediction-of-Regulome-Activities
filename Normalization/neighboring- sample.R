#### This is the code for neighboring-sample normalization

#### select top 5 neighboring samples
cormax=0
id.micro=NULL
celltypename=NULL
cor_val=NULL
exon.celltype=NULL
num=0

for (i in 1:57){
  flag=rep(0,ncol(micro.all))
  for (k in 1:5){
    num=num+1
    cormax=0
    for(j in 1:(ncol(micro.all))){
      if ((cor(as.numeric(exon.match[,i]),as.numeric(micro.all[,j]))>cormax) & (flag[j]== 0)){
        cormax=cor(as.numeric(exon.match[,i]),as.numeric(micro.all[,j]))
        cor_val[num]=cormax
        celltypename[num]=colnames(micro.all)[j]
        id.micro[num]=j
        exon.celltype[num]=colnames(exon.match)[i]
      }
    }
    flag[id.micro[num]]=1
  }
}

result.top5=as.data.frame(cbind(celltypename,cor_val,id.micro,exon.celltype)) # show the results, including cell type name, the correlation, and corresponding exon cell type
micro.uni=micro.all[,unique(id.micro)] # neighboring samples 


#### normalization
exon.match.m=as.matrix(exon.match)
neighboring_sample_normalize=function(affydata){
  affy.trans=matrix(0,nrow=nrow(exon.match),ncol=ncol(affydata))
  exon.mean=NULL
  exon.sd=NULL
  micro.mean=NULL
  micro.sd=NULL
  b0=NULL
  a0=NULL
  for (i in 1:nrow(exon.match)){
    exon.mean[i]=mean(exon.match.m[i,])
    exon.sd[i]=sd(exon.match.m[i,])
    micro.mean[i]=mean(affydata[i,])
    micro.sd[i]=sd(affydata[i,])
    b0[i]=micro.sd[i]/exon.sd[i]
    a0[i]=micro.mean[i]-b0[i]*exon.mean[i]
    affy.trans[i,]=(affydata[i,]-a0[i])/(b0[i])
  }
  ans=list(dat=affy.trans,m=a0,se=b0)
  return(ans)
}

micro.trans.result=neighboring_sample_normalize(as.matrix(micro.uni))
micro.trans=micro.trans.result$dat



#### get the microarray test set after normalization
micro.norm=matrix(0,nrow(micro.raw),ncol(micro.raw))
a0=micro.trans.result$m
b0=micro.trans.result$se
for(i in 1:nrow(exon.match)){
  micro.norm[i,]=as.numeric((micro.raw[i,]-a0[i])/b0[i])
}




