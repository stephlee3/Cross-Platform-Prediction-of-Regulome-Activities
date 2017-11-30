# Cross-Platform Prediction of Regulome Activities

## Overview
This is the project for Cross-Platform Predicition of Regulome Activities based on gene expression data. The goal is to use human exon array data to train the prediction model, then apply it to the gene expression data from different platforms.(e.g. affymetrix microarray data from GEO, or RNA-seq data from ENCODE project).

## Data Description
- DNase-seq data
The bowtie aligned(alignment based on hg19) DNase-seq data for 57 human cell types with normal karyotype were downloaded from the [ENCODE](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwDnase) in bam format.

- exon array gene expression data
The Affymetrix Human Exon 1.0 ST Array(i.e. exon array) data for the same 57 ENCODE cell types were downloaded from GEO([GEO accession number: GSE19090](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse19090)). All samples were processed using the GeneBASE software to compute gene-level expression. 

- microarray gene expression data
The microarray gene expression data is downloaded from the [Gene Expression BARCODE project](http://barcode.luhs.org/). We use the data from Affymetrix Human Genome U133A Array(GPL96), a compendium of 11778 samples. The processed data can be downloaded from [R Package](http://www.bioconductor.org/packages/release/data/experiment/html/Affyhgu133aExpr.html).

## Normalization
We need to normalize the data from different platforms. Here we introduce four different normalization methods to tackle platform effect, check the [Normalization](https://github.com/stephlee3/Cross-Platform-Prediction-of-Regulome-Activities/tree/master/Normalization) part for more details.

## Model Training
We build the prediction model based on [BIRD](https://github.com/WeiqiangZhou/BIRD), using gene expression data to predict DH level. Please check the [model_train](https://github.com/stephlee3/Cross-Platform-Prediction-of-Regulome-Activities/tree/master/model_train) part for more details.

## Prediction Result Analysis
To evaluate the prediction performance, please check the [Prediction](https://github.com/stephlee3/Cross-Platform-Prediction-of-Regulome-Activities/tree/master/Prediction) part for more details.



