#### This is the code for sample-quantile normalization

#### compute the quantile of exon array data
compute_quantile=function(train_exon){
  data_sorted=train_exon
  for(i in 1:dim(train_exon)[2]){
    data_sorted[,i]=sort(train_exon[,i],decreasing=TRUE)
  }
  data_quantile=apply(data_sorted,1,mean)
  return(data_quantile)
}

exon.quantile=compute_quantile(exon.match)
exon.quantile=as.numeric(sort(exon.quantile))

#### transform the microarray data by quantile normalization
quantile_map=function(mdata){
  qdata=matrix(0,nrow(mdata),ncol(mdata))
  mdata_sort=matrix(0,nrow(mdata),ncol(mdata))
  for (i in 1:ncol(mdata)){
    mdata_sort[,i]=sort(mdata[,i])
  }
  for (j in 1:ncol(mdata)){
    r=rank(mdata[,j])
    for(i in 1:nrow(mdata)){
      qdata[i,j]=exon.quantile[r[i]]
    }
  }
  return(qdata)
}
micro.sample.quantile=quantile_map(micro.raw)
write.table(micro.sample.quantile,"micro_sample_quantile.txt",row.names=F)
