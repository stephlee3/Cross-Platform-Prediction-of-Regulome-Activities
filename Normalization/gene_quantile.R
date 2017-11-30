#### Fig3: gene-quantile normalization: gene-quantile


#### read the data
exon.match=read.table("exon_match.txt",header=T)
micro.sel=read.table("micro_sel.txt",header=T)

#### loess_plot to generate the plot and get normalized microarray test data
micro.test.quantile=matrix(0,11820,7)
loess_plot=function(i){
  exon.sort=as.numeric(sort(exon.match[i,]))
  micro.sort=as.numeric(sort(micro.sel[i,]))
  plot(exon.sort~micro.sort,cex.lab=1.2,cex.axis=1.2)
  exon.model=loess(exon.sort~micro.sort,span=0.8)
  exon.fitted=exon.model$fitted
  micro.test.quantile[i,]=predict(exon.model,as.numeric(micro.raw[i,]))
  lines(exon.fitted~micro.sort,col="red",lwd=3,lty=1)
}

write.table(micro.test.quantile,"micro_gene_quantile.txt",row.names = F)

