## Get the data from Geoff
setwd("~/Desktop/CRUK")
## Installing paclages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicRanges")
install.packages("devtools")
library(devtools)
devtools::install_github("friend1ws/pmsignature")
install.packages("ggplot2")
install.packages("Rcpp")
# Load packages and functions
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(QDNAseq))
suppressMessages(library(flexmix))
suppressMessages(library(NMF))
suppressMessages(source("britroc-cnsignatures-bfb69cd72c50/main_functions.R"))
suppressMessages(source("britroc-cnsignatures-bfb69cd72c50/helper_functions.R"))
source("get10Mbfeatures.R")
source("getFeatureVectorAndCountData.R")
library(pmsignature)
library(Rcpp)

## Get data from Geoff
#Basic functions
samp_annotation<-read.table("britroc-cnsignatures-bfb69cd72c50/manuscript_Rmarkdown/data/britroc_sample_data.csv",stringsAsFactors = F,sep = ",",header=T)
chrlen<-read.table(paste(this_path,"data/hg19.chrom.sizes.txt",sep="/"),sep="\t",stringsAsFactors = F)[1:24,]

result<- samp_annotation %>% 
  filter(!samp_annotation$Failed=="Y") %>%
  dplyr::select(Britroc_No,IM.JBLAB_ID,star_rating) %>%
  dplyr::group_by(Britroc_No) %>%
  dplyr::slice(which.max(star_rating))

all_CN<-readRDS("britroc-cnsignatures-bfb69cd72c50/manuscript_Rmarkdown/data/britroc_absolute_copynumber.rds")
all_CN<-all_CN[,colnames(all_CN)%in%samp_annotation[!samp_annotation$Failed=="Y","IM.JBLAB_ID"]]
ids<-result %>% filter(star_rating==3)
ids<-ids$IM.JBLAB_ID
hq_CN<-all_CN[,colnames(all_CN)%in%ids]

## Use Shiraishi's modified functions
setwd("~/Desktop/CRUK")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicRanges")
install.packages("devtools")
install.packages("ggplot2")
install.packages("Rcpp")
library(devtools)
devtools::install_github("friend1ws/pmsignature")
library(pmsignature)
library(Rcpp)

## Get the necessary data to find the signatures
MF0<-get10Mbfeatures(all_CN,chrlen,0)
MFM<-MF0[,c(1,3:9)]
#Organize some of the columns into bins with fewer categories
Col2_new<-MFM$Col2
Col2_new[which(Col2_new>5)]<-6
MFM$Col2<-Col2_new
Col5_new<-MFM$Col5
Col5_new[which(Col5_new>6 & Col5_new<=10)]<-7
Col5_new[which(Col5_new>10)]<-8
MFM$Col5<-Col5_new
Col6_new<-MFM$Col6
Col6_new[which(Col6_new>6 & Col6_new<=10)]<-7
Col6_new[which(Col6_new>10)]<-8
MFM$Col6<-Col6_new
Col7_new<-MFM$Col7
Col7_new[which(Col6_new>6)]<-7
MFM$Col7<-Col7_new

#Force all entries to be bigger than 0 to create the MutationFeatureData
Final_MFM<-MFM
Final_MFM[,2:8]<-MFM[,2:8]+1

#Vector of Possible Features (POFE)
POFE<-c()
for(j in 2:8){
  POFE<-c(POFE,max(Final_MFM[,j]))
}

#Matrices of Feature Vector List (FEVE) and of Count Data (CODA)
FVCD<-getFeatureVectorAndCountData(Final_MFM)
FEVE<-FVCD[[1]]
CODA<-FVCD[[2]]

#Signature identification
MFD<-new(Class="MutationFeatureData",featureVectorList=FEVE,
    sampleList=as.character(unique(Final_MFM[,1])),countData=CODA, 
    possibleFeatures=as.integer(POFE))
NOS<-6 #Number of Signatures
MFD_Param<-getPMSignature(MFD,NOS)

#Signature by sample 
heatmap(MFD_Param@sampleSignatureDistribution)

#Signature by Component 
LMat<-list()
NZrows<-c()
NZcols<-c()
for(i in 1:NOS){
  Mat<-getSignatureValue(MFD_Param, i)
  Mat[which(Mat<0.05)]<-0
  LMat[[i]]<-Mat
  NZrows<-c(NZrows,which(Mat>0,arr.ind = TRUE)[,1])
  NZcols<-c(NZcols,which(Mat>0,arr.ind = TRUE)[,2])
}
#Selected components
selcomp<-unique(matrix(c(NZrows,NZcols),ncol=2))
#reorder in ascending
selcomp<-selcomp[order(selcomp[,1],selcomp[,2]),]
MATNZ<-matrix(0,nrow=NOS,ncol=max(POFE))
MATNZ[selcomp]<-1
MATNZ #This gives an idea of some feature components that could have been ignored
SignatureByComponent<-matrix(0,nrow=NOS,ncol=nrow(selcomp))
for(i in 1:NOS){
  Mat<-getSignatureValue(MFD_Param, i)
  for(j in 1:nrow(selcomp)){
    SignatureByComponent[i,j]<-Mat[selcomp[j,1],selcomp[j,2]]
  }
}
heatmap(SignatureByComponent,Rowv = NA,Colv = NA)