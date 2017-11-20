# read the data
exon.data=read.table("exon_match.txt",header=T)
micro.data=read.table("micro_match.txt",header=T)

# test and train partitioning
id.test=c(8,9,24,43,45,47,52)
id.train=as.matrix(c(1:57))[-id.test]
write.table(id.test,"id_test.txt")

# select genes
cor.tmp=read.table("cross_cor.txt",header=T)
p=c(-1,0,0.1,0.2,0.3,0.4,0.5,0.6)
sel.gene=which(cor.tmp>p[5])

# train and test sets
exon.train=exon.data[sel.gene,id.train]
exon.test=exon.data[sel.gene,id.test]
micro.test=micro.data[sel.gene,]

write.table(exon.train,"exon_train.txt")
write.table(exon.test,"exon_test.txt")
write.table(micro.test,"micro_test.txt")

# train the model
# get model data
get_regression_info <- function(DNase_train,Exon_train_mean,top_n){
  y_cor <- cor(DNase_train,t(Exon_train_mean))
  y_cor[is.na(y_cor)] <- 0
  
  if(top_n > 1){
    max_idx <- sort(y_cor, decreasing=TRUE, index.return=TRUE)$ix[1:top_n]
    data_train <- data.frame(y=DNase_train, t(Exon_train_mean[max_idx,]))
    fit <- lm(y~.,data_train)
  }else{
    max_idx <- which(y_cor == max(y_cor))[1]
    data_train <- data.frame(y=DNase_train, x=Exon_train_mean[max_idx,])
    fit <- lm(y~x,data_train)
  }
  
  return(list(lm_fit=fit,predictor=max_idx))
}

##standardize data
standardize_row_train <- function (Data_train)
{
  train_mean_sd <- matrix(data = NA, nrow = dim(Data_train)[1],
                          ncol = 2)
  Data_train_sd <- matrix(data = NA, nrow = dim(Data_train)[1],
                          ncol = dim(Data_train)[2])
  for (i in 1:dim(Data_train)[1]) {
    train_mean_sd[i, 1] <- mean(Data_train[i, ])
    train_mean_sd[i, 2] <- sd(Data_train[i, ])
    if (train_mean_sd[i, 2] == 0) {
      Data_train_sd[i, ] <- rep(0, dim(Data_train)[2])
    }
    else {
      Data_train_sd[i, ] <- (Data_train[i, ] - train_mean_sd[i,
                                                             1])/train_mean_sd[i, 2]
    }
  }
  result <- list(mean_sd = train_mean_sd, train = Data_train_sd)
}

standardize_rows <- function(Data_train, Data_test){
  
  train_mean_sd <- matrix(data=NA, nrow=dim(Data_train)[1], ncol=2)
  Data_train_sd <- matrix(data=NA, nrow=dim(Data_train)[1], ncol=dim(Data_train)[2])
  Data_test_sd <- matrix(data=NA, nrow=dim(Data_test)[1], ncol=dim(Data_test)[2])
  for(i in 1:dim(Data_train)[1]){
    train_mean_sd[i,1] <- mean(Data_train[i,])
    train_mean_sd[i,2] <- sd(Data_train[i,])
    if(train_mean_sd[i,2]==0){
      Data_train_sd[i,] <- rep(0,dim(Data_train)[2])
      Data_test_sd[i,] <- rep(0,dim(Data_test)[2])
    }
    else{
      Data_train_sd[i,] <- (Data_train[i,] - train_mean_sd[i,1])/train_mean_sd[i,2]
      Data_test_sd[i,] <- (Data_test[i,] - train_mean_sd[i,1])/train_mean_sd[i,2]
    }
  }
  
  result <- list("mean_sd"=train_mean_sd, "train"=Data_train_sd, "test"=Data_test_sd)
  
}

##get cluster mean##
cluster_data_mean <- function(Exon_data,cluster){
  
  cluster_unique <- unique(cluster)
  Exon_mean <- matrix(data=NA, nrow=length(cluster_unique), ncol=dim(Exon_data)[2])
  
  
  for(j in 1:length(cluster_unique)){
    cluster_idx <- which(cluster==j)
    Exon_cluster <- Exon_data[cluster_idx,]
    Exon_cluster <- as.matrix(Exon_cluster)
    
    for (i in 1:dim(Exon_data)[2]){
      if(length(Exon_cluster)==dim(Exon_data)[2]){
        Exon_mean[j,i] <- Exon_cluster[i]
      }
      else{
        Exon_mean[j,i] <- mean(Exon_cluster[,i])
      }
    }
  }
  
  return(Exon_mean)
}

##get data quantile##
compute_quantile <- function(train_exon){
  data_sorted <- train_exon
  for(i in 1:dim(train_exon)[2]){
    data_sorted[,i] <- sort(train_exon[,i],decreasing=TRUE)
  }
  data_quantile <- apply(data_sorted,1,mean)
  return(data_quantile)
}

##main script##
Num_cluster <- 500   ##different number of cluster can be set here
Num_predictor <- 8  ##different number of cluster can be set here


##input data are Exon_data_sample.txt and DNase_data_sample.rda
Exon_train <- as.matrix(read.table(file="exon_train.txt",header=TRUE,row.names=1))  ##read exon train gene expression data
exon_test <- as.matrix(read.table(file="exon_test.txt",header=TRUE,row.names=1))  ##read exon test gene expression data
micro_test <- as.matrix(read.table(file="micro_test.txt",header=TRUE,row.names=1))  ##read micro test gene expression data

Exon_processed <- standardize_rows(Exon_train,exon_test) #standardize gene expression data
Exon_train_sd <- Exon_processed$train
exon_test_sd=Exon_processed$test

micro_processed=standardize_rows(Exon_train,micro_test)
micro_test_sd=micro_processed$test

set.seed(2017)
Exon_cluster <- kmeans(Exon_train_sd, centers=Num_cluster, nstart=10, iter.max=50)    ##cluster gene expression
cluster_idx=Exon_cluster$cluster
write.table(Exon_cluster$cluster,file="cluster_idx.txt",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

Exon_train_mean <- cluster_data_mean(Exon_train_sd,Exon_cluster$cluster) #calculated average gene expression within each cluster
micro_test_mean=cluster_data_mean(micro_test_sd,Exon_cluster$cluster)
exon_test_mean=cluster_data_mean(exon_test_sd,Exon_cluster$cluster)

write.table(Exon_train_mean,"exon_train_mean.txt",row.names = F)
write.table(micro_test_mean,"micro_test_mean.txt",row.names = F)
write.table(exon_test_mean,"exon_test_mean.txt",row.names = F)


##Locus-level prediction##

load("DNase_data.rda")  ##read DNase-seq data, first three columns contain the genomic locus information
DNase_processed <- standardize_row_train(as.matrix(DNase_data[,id.train+3])) #standardize DNase-seq data
DNase_mean_sd <- DNase_processed$mean_sd
DNase_train_sd <- DNase_processed$train
write.table(DNase_mean_sd,"DNase_mean_sd.txt",row.names = F)
write.table(DNase_data[,id.test+3],"DNase_c.txt",row.names = F)

#DNase_train_sd=as.matrix(read.table("DNase_train_sd.txt")

coef_all <- matrix(data=NA,nrow=dim(DNase_train_sd)[1],ncol=Num_predictor+1)
predictor_idx <- matrix(data=NA,nrow=dim(DNase_train_sd)[1],ncol=Num_predictor)

for (i in 1:dim(DNase_train_sd)[1]){
  regress_result <- get_regression_info(DNase_train_sd[i,],Exon_train_mean,Num_predictor) #build regression model for each DHS
  coef_all[i,] <- coefficients(regress_result$lm_fit) #get regression coefficients
  predictor_idx[i,] <- regress_result$predictor #get predictor index
}

write.table(coef_all[,-1],file="regress_coef.txt",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE) ##b0 is zero and excluded
write.table(predictor_idx,file="regress_predictor.txt",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)



