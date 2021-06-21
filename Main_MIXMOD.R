## Get the data from Geoff
setwd("~/Desktop/CRUK")
## Installing paclages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicRanges")
BiocManager::install("GWASTools")
install.packages("devtools")
library(devtools)
devtools::install_github("friend1ws/pmsignature")
install.packages("ggplot2")
install.packages("Rcpp")
install.packages('gplots')
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
library(GWASTools)
library(Rcpp)
library(gplots)

## Get data from Geoff
#Basic functions
samp_annotation<-read.table("britroc-cnsignatures-bfb69cd72c50/manuscript_Rmarkdown/data/britroc_sample_data.csv",stringsAsFactors = F,sep = ",",header=T)
#Chromosome lengths
chrlen<-read.table(paste(this_path,"data/hg19.chrom.sizes.txt",sep="/"),sep="\t",stringsAsFactors = F)[1:24,]
#Centromeres positions
C19<-centromeres.hg19
C19$centromere<-floor((C19$left.base+C19$right.base)/2)

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
#Get the necessary data to find the signatures
MF0<-get10Mbfeatures(all_CN,chrlen,0)
N<-10
MF5<-getNMbfeatures(all_CN,chrlen,C19,N,0)
MF10hq<-getNMbfeatures(hq_CN,chrlen,C19,N,0)
MFM<-MF10hq[,c(1,3:10)]
MFM<-MFM[-which(MFM$Col2==0&MFM$Col3==0&MFM$Col4==0),]


Feature_short_names<-c('ID','BPC','CNLP','CNGP','SCNCH','MAXCN','MINCN','WSS','CEND')
colnames(MFM)<-Feature_short_names
CN_SegmentFeatures<-list()
for(j in 2:ncol(MFM)){
  CN_SegmentFeatures[[j-1]]<-data.frame(MFM[,c(1,j)])
}
names(CN_SegmentFeatures)<-Feature_short_names[2:length(Feature_short_names)]
## Apply mixture models
seed=77777
min_prior=0.001
model_selection="BIC"
nrep=1
niter=10000

dat<-as.numeric(CN_SegmentFeatures[["BPC"]][,2])
BPC_mm<-fitComponent(dat,seed=seed,dist="pois",model_selection=model_selection,
                       min_prior=min_prior,niter=niter,nrep=nrep,min_comp=2,max_comp=8)

dat<-as.numeric(CN_SegmentFeatures[["SCNCH"]][,2])
SCNCH_mm<-fitComponent(dat,seed=seed,dist="pois",model_selection=model_selection,
                      min_prior=min_prior,niter=niter,nrep=nrep,min_comp=2,max_comp=8)

dat<-as.numeric(CN_SegmentFeatures[["MAXCN"]][,2])
MAXCN_mm<-fitComponent(dat,seed=seed,dist="pois",model_selection=model_selection,
                     min_prior=min_prior,niter=niter,nrep=nrep,min_comp=2,max_comp=8)

dat<-as.numeric(CN_SegmentFeatures[["MINCN"]][,2])
MINCN_mm<-fitComponent(dat,seed=seed,dist="pois",model_selection=model_selection,
                      min_prior=min_prior,niter=niter,nrep=nrep,min_comp=2,max_comp=8)

dat<-as.numeric(CN_SegmentFeatures[["WSS"]][,2])
WSS_mm<-fitComponent(dat,seed=seed,dist="pois",model_selection=model_selection,
                     min_prior=min_prior,niter=niter,nrep=nrep,min_comp=2,max_comp=8)

dat<-as.numeric(CN_SegmentFeatures[["CEND"]][,2])
CEND_mm<-fitComponent(dat,seed=seed,dist="pois",model_selection=model_selection,
                      min_prior=min_prior,niter=niter,nrep=nrep,min_comp=2,max_comp=8)

MFM$CNLP<-floor(MFM$CNLP*5)+1
MFM$CNGP<-floor(MFM$CNGP*5)+1
Col6_new<-MFM$MAXCN
Col6_new[which(Col6_new>6 & Col6_new<=10)]<-7
Col6_new[which(Col6_new>10)]<-8
MFM$MAXCN<-Col6_new+1
Col7_new<-MFM$MINCN
Col7_new[which(Col6_new>6)]<-7
MFM$MINCN<-Col7_new+1
Col9_new<-MFM$CEND
Col9_new[which(Col9_new<=3e7)]<-1
Col9_new[which(Col9_new>3e7 & Col9_new<=8e7)]<-2
Col9_new[which(Col9_new>8e7 & Col9_new<=1e8)]<-3
Col9_new[which(Col9_new>1e8)]<-4
MFM$CEND<-Col9_new

Final_MFM<-data.frame(MFM$ID,BPC_mm@cluster,MFM$CNLP,MFM$CNGP,SCNCH_mm@cluster,
                      MFM$MAXCN,MFM$MINCN,WSS_mm@cluster,MFM$CEND)








#MFM<-MF5[,c(1,3:10)]
#MF10<-MFM
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
Col8_new<-MFM$Col8
Col8_new[which(Col8_new<=10)]<-0
Col8_new[which(Col8_new>10&Col8_new<=20)]<-1
Col8_new[which(Col8_new>20&Col8_new<=30)]<-2
Col8_new[which(Col8_new>30&Col8_new<=40)]<-3
Col8_new[which(Col8_new>40&Col8_new<=50)]<-4
Col8_new[which(Col8_new>50)]<-5
MFM$Col8<-Col8_new
MFM$Col9<-round(MFM$Col9/2,1)*10

#Remove segments without aberrations
Final_MFM<-MFM[-which(MFM$Col2==0&MFM$Col3==0&MFM$Col4==0),]
#Force all entries to be bigger than 0 to create the MutationFeatureData
Final_MFM[,2:9]<-Final_MFM[,2:9]+1

#Vector of Possible Features (POFE)
POFE<-c()
for(j in 2:9){
  POFE<-c(POFE,max(Final_MFM[,j]))
}

### NOTE: Move to Lena's alternative instead 
# 
# #Matrices of Feature Vector List (FEVE) and of Count Data (CODA)
# FVCD<-getFeatureVectorAndCountData(Final_MFM)
# FEVE<-FVCD[[1]]
# CODA<-FVCD[[2]]
# 
# #Signature identification
# MFD<-new(Class="MutationFeatureData",featureVectorList=FEVE,
#     sampleList=as.character(unique(Final_MFM[,1])),countData=CODA, 
#     possibleFeatures=as.integer(POFE))
# NOS<-6 #Number of Signatures
# MFD_Param<-getPMSignature(MFD,NOS)
# 
# #Signature by sample 
# heatmap(MFD_Param@sampleSignatureDistribution)
# 
# #Signature by Component 
# LMat<-list()
# NZrows<-c()
# NZcols<-c()
# for(i in 1:NOS){
#   Mat<-getSignatureValue(MFD_Param, i)
#   Mat[which(Mat<0.05)]<-0
#   LMat[[i]]<-Mat
#   NZrows<-c(NZrows,which(Mat>0,arr.ind = TRUE)[,1])
#   NZcols<-c(NZcols,which(Mat>0,arr.ind = TRUE)[,2])
# }
# #Selected components
# selcomp<-unique(matrix(c(NZrows,NZcols),ncol=2))
# #reorder in ascending
# selcomp<-selcomp[order(selcomp[,1],selcomp[,2]),]
# MATNZ<-matrix(0,nrow=NOS,ncol=max(POFE))
# MATNZ[selcomp]<-1
# MATNZ #This gives an idea of some feature components that could have been ignored
# SignatureByComponent<-matrix(0,nrow=NOS,ncol=nrow(selcomp))
# for(i in 1:NOS){
#   Mat<-getSignatureValue(MFD_Param, i)
#   for(j in 1:nrow(selcomp)){
#     SignatureByComponent[i,j]<-Mat[selcomp[j,1],selcomp[j,2]]
#   }
# }
# heatmap(SignatureByComponent,Rowv = NA,Colv = NA)
# visMembership(MFD, MFD_Param)
# 

##Lena's alternative
MFM_combinations <- apply(Final_MFM[,2:9], 1, paste0, collapse='-')
num_combinations <- table(MFM_combinations)
featureVectorList_mod <- sapply(names(num_combinations), function(j){
  as.numeric(strsplit(j, split = "-")[[1]])})

countData_mod0 <- apply(Final_MFM[,2:9], 1, function(j) which(names(num_combinations) == paste0(j, collapse = '-')))
## adding the numerical id of each patient
countData_mod0 <- cbind(patientid=countData_mod0,
                        cat_combination=as.numeric(as.factor(Final_MFM[,1])),
                        count=1)
##' counting the number of instances in each category of combination of features
##' for each patient
countData_mod <- aggregate(count~patientid+cat_combination, countData_mod0, sum)
countData_Lena <- t(countData_mod)


MFD_Lena<-new(Class="MutationFeatureData",featureVectorList=featureVectorList_mod,
              sampleList=as.character(unique(Final_MFM[,1])),countData=countData_Lena, 
              possibleFeatures=as.integer(POFE))
NOS<-5 #Number of Signature
MFD_Param_Lena<-getPMSignature(MFD_Lena,K=NOS,numInit = 30)

#Signature by sample 
heatmap(MFD_Param_Lena@sampleSignatureDistribution,scale='none')

#Signature by Component 
LMat<-list()
NZrows<-c()
NZcols<-c()
for(i in 1:NOS){
  Mat<-getSignatureValue(MFD_Param_Lena, i)
  Mat[which(Mat<0.05)]<-0
  LMat[[i]]<-Mat
  NZrows<-c(NZrows,which(Mat>0,arr.ind = TRUE)[,1])
  NZcols<-c(NZcols,which(Mat>0,arr.ind = TRUE)[,2])
}
#Selected components
selcomp<-unique(matrix(c(NZrows,NZcols),ncol=2))
#reorder in ascending
selcomp<-selcomp[order(selcomp[,1],selcomp[,2]),]
selcompNames<-c()
Feature_names<-c('Bp count','CN loss %','CN gain %','Sum CN changes','Max CN',
                 'Min CN','Weighted Sement sum','Centromere distance')
for(i in 1:nrow(selcomp)){
  CurrentFeature<-Feature_names[selcomp[i,1]]
  selcompNames<-c(selcompNames,paste(CurrentFeature,selcomp[i,2]))
}
# MATNZ<-matrix(0,nrow=NOS,ncol=max(POFE))
# MATNZ[selcomp]<-1
# MATNZ #This gives an idea of some feature components that could have been ignored
SignatureByComponent<-matrix(0,nrow=NOS,ncol=nrow(selcomp))
for(i in 1:NOS){
  Mat<-getSignatureValue(MFD_Param_Lena, i)
  for(j in 1:nrow(selcomp)){
    SignatureByComponent[i,j]<-Mat[selcomp[j,1],selcomp[j,2]]
  }
}
colnames(SignatureByComponent)<-selcompNames
heatmap(SignatureByComponent,scale='none')
visMembership(MFD_Lena, MFD_Param_Lena)