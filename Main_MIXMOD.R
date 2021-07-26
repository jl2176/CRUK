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
#suppressMessages(library(dbplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(QDNAseq))
suppressMessages(library(flexmix))
suppressMessages(library(NMF))
suppressMessages(source("britroc-cnsignatures-bfb69cd72c50/main_functions.R"))
suppressMessages(source("britroc-cnsignatures-bfb69cd72c50/helper_functions.R"))
source("getNMbfeatures.R")
source("getFeatureVectorAndCountData.R")
library(pmsignature)
library(GWASTools)
library(Rcpp)
library(gplots)
library(corrplot)
library(RColorBrewer)

## Get data from Geoff
#Basic functions
samp_annotation<-read.table("data/britroc_sample_data.csv",stringsAsFactors = F,sep = ",",header=T)
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

all_CN<-readRDS("data/britroc_absolute_copynumber.rds")
all_CN<-all_CN[,colnames(all_CN)%in%samp_annotation[!samp_annotation$Failed=="Y","IM.JBLAB_ID"]]
ids<-result %>% filter(star_rating==3)
ids<-ids$IM.JBLAB_ID
hq_CN<-all_CN[,colnames(all_CN)%in%ids]

## Use Shiraishi's modified functions
#Get the necessary data to find the signatures
#MF0<-get10Mbfeatures(all_CN,chrlen,0)
N<-10
#MF5<-getNMbfeatures(all_CN,chrlen,C19,N,0)
MF10hq<-getNMbfeatures(hq_CN,chrlen,C19,N,0)
MFM<-MF10hq[,c(1,3:11)]
MFM<-MFM[-which(MFM$Col2==0&MFM$Col8==2e7),]
#MFM[is.nan(MFM[,5]),5]<-0


Feature_short_names<-c('ID','BPC','CNLP','CNGP','SCNCH','MAXCN','MINCN','WSS','CEND','osCN')
colnames(MFM)<-Feature_short_names
CN_SegmentFeatures<-list()
for(j in 2:ncol(MFM)){
  CN_SegmentFeatures[[j-1]]<-data.frame(MFM[,c(1,j)])
}
names(CN_SegmentFeatures)<-Feature_short_names[2:length(Feature_short_names)]

## Correlations between features
heatmap(abs(cor(MFM[,-c(1,9)])),scale='none')

## Apply mixture models
seed=77777
min_prior=0.001
model_selection="BIC"
nrep=1
niter=10000

dat<-as.numeric(CN_SegmentFeatures[["BPC"]][,2])
BPC_mm<-fitComponent(dat,seed=seed,dist="pois",model_selection=model_selection,
                       min_prior=min_prior,niter=niter,nrep=nrep,min_comp=2,max_comp=8)

BPC_com<-parameters(BPC_mm)
BPC_sortcom<-sort(BPC_com)
BPC_finalcluster<-BPC_mm@cluster
for(k in 1:length(BPC_com)){
  mk<-match(BPC_sortcom[k],BPC_com)
  BPC_finalcluster[which(BPC_finalcluster==mk)]<-k+length(BPC_com)
}
BPC_finalcluster<-BPC_finalcluster-length(BPC_com)
MFM$BPC[which(MFM$BPC>3)]<-4
MFM$BPC<-MFM$BPC+1

# dat<-as.numeric(CN_SegmentFeatures[["CNLP"]][,2])
# CNLP_mm<-fitComponent(dat,seed=seed,model_selection=model_selection,
#                        min_prior=min_prior,niter=niter,nrep=nrep,min_comp=2,max_comp=4)
# 
# dat<-as.numeric(CN_SegmentFeatures[["CNGP"]][,2])
# CNGP_mm<-fitComponent(dat,seed=seed,model_selection=model_selection,
#                        min_prior=min_prior,niter=niter,nrep=nrep,min_comp=2,max_comp=4)

dat<-as.numeric(CN_SegmentFeatures[["SCNCH"]][,2])
SCNCH_mm<-fitComponent(dat,dist="pois",seed=seed,model_selection=model_selection,
                       min_prior=min_prior,niter=niter,nrep=nrep,min_comp=2,max_comp=8)

SCNCH_com<-parameters(SCNCH_mm)
SCNCH_sortcom<-sort(SCNCH_com)
SCNCH_finalcluster<-SCNCH_mm@cluster
for(k in 1:length(SCNCH_com)){
  mk<-match(SCNCH_sortcom[k],SCNCH_com)
  SCNCH_finalcluster[which(SCNCH_finalcluster==mk)]<-k+length(SCNCH_com)
}
SCNCH_finalcluster<-SCNCH_finalcluster-length(SCNCH_com)

dat<-as.numeric(CN_SegmentFeatures[["MAXCN"]][,2])
MAXCN_mm<-fitComponent(dat,seed=seed,model_selection=model_selection,
                     min_prior=min_prior,niter=niter,nrep=nrep,min_comp=2,max_comp=6)
MFM$MAXCN[which(MFM$MAXCN>5)]<-6
MFM$MAXCN[which(MFM$MAXCN==0)]<-1

dat<-as.numeric(CN_SegmentFeatures[["MEANCN"]][,2])
MEANCN_mm<-fitComponent(dat,seed=seed,model_selection=model_selection,
                      min_prior=min_prior,niter=niter,nrep=nrep,min_comp=2,max_comp=6)

MFM$WSS[which(MFM$WSS<=1e7)]<-1
MFM$WSS[which(MFM$WSS>1e7&MFM$WSS<=2e7)]<-2
MFM$WSS[which(MFM$WSS>2e7&MFM$WSS<=3e7)]<-3
MFM$WSS[which(MFM$WSS>3e7&MFM$WSS<=4e7)]<-4
MFM$WSS[which(MFM$WSS>4e7&MFM$WSS<=5e7)]<-5
MFM$WSS[which(MFM$WSS>5e7)]<-6


dat<-as.numeric(CN_SegmentFeatures[["WSS"]][,2])
dat<-dat[which(dat>0)]
WSS_mm<-fitComponent(dat,seed=seed,model_selection=model_selection,
                     min_prior=min_prior,niter=niter,nrep=nrep,min_comp=2,max_comp=8)

# WSS_com<-parameters(WSS_mm)
# WSS_sortcom<-sort(WSS_com)
# WSS_finalcluster<-WSS_mm@cluster
# for(k in 1:length(WSS_com)){
#   mk<-match(WSS_sortcom[k],WSS_com)
#   WSS_finalcluster[which(WSS_finalcluster==mk)]<-k+length(WSS_com)
# }
# WSS_finalcluster<-WSS_finalcluster-length(WSS_com)

dat<-as.numeric(CN_SegmentFeatures[["CEND"]][,2])
CEND_mm<-fitComponent(dat,seed=seed,model_selection=model_selection,
                      min_prior=min_prior,niter=niter,nrep=nrep,min_comp=2,max_comp=8)


dat<-as.numeric(CN_SegmentFeatures[["osCN"]][,2])
MFM$osCN[which(MFM$osCN>=1)]<-2
MFM$osCN[which(MFM$osCN==0)]<-1

MFM$CNLP<-MFM$CNLP/1e7
MFM$CNLP[MFM$CNLP==1]<-3
MFM$CNLP[MFM$CNLP<1&MFM$CNLP>0]<-2
MFM$CNLP[MFM$CNLP==0]<-1

MFM$CNGP[MFM$CNGP>=0.3]<-2
MFM$CNGP[MFM$CNGP<0.3]<-1
# Col6_new<-MFM$MAXCN
# Col6_new[which(Col6_new<5)]<-1
# Col6_new[which(Col6_new>=5)]<-2
# MFM$MAXCN<-Col6_new
# Col7_new<-MFM$MINCN
# Col7_new[which(Col7_new>6)]<-7
# MFM$MINCN<-Col7_new+1
Col9_new<-MFM$CEND
Col9_new[which(Col9_new<=1e7)]<-1
Col9_new[which(Col9_new>1e7 & Col9_new<=5e7)]<-2
Col9_new[which(Col9_new>5e7 & Col9_new<=9e7)]<-3
Col9_new[which(Col9_new>9e7)]<-4
MFM$CEND<-Col9_new

# Final_MFM<-data.frame(MFM$ID,BPC_mm@cluster,MFM$CNLP,MFM$CNGP,SCNCH_mm@cluster,
#                       MFM$MAXCN,MFM$MINCN,WSS_mm@cluster,MFM$CEND)

set.seed(7777)

#Set the Matrix with the base features
NF<-3 #number of features
 # Final_MFM<-data.frame(MFM$ID,MFM$BPC,MFM$MAXCN,
 #                       MFM$WSS,MFM$CEND)
Final_MFM<-data.frame(MFM$ID,MFM$BPC,MFM$CNLP,
                      MFM$MAXCN)

heatmap(abs(cor(Final_MFM[,-1])),scale='none')

# #And the randomize version
# Random_Final_MFM<-Final_MFM
# Random_Final_MFM[,2:(NF+1)]<-randomize(Random_Final_MFM[,2:(NF+1)])

#Vector of Possible Features (POFE)
POFE<-c()
for(j in 2:(NF+1)){
  POFE<-c(POFE,max(Final_MFM[,j]))
}

#Sample by component matrix
SampleByComponent<-matrix(0,nrow=length(unique(Final_MFM[,1])),ncol=sum(POFE))
cur_row<-0
for(i in unique(Final_MFM[,1])){
  cur_row<-cur_row+1
  cur_pos<-1
  for(j in 1:NF){
    SampleByComponent[cur_row,cur_pos:(cur_pos+POFE[j]-1)]<-
      table(c(1:POFE[j],Final_MFM[which(Final_MFM[,1]==i),j+1]))-1
    cur_pos<-cur_pos+POFE[j]
  }
}
CompNames<-c()
Feature_names<-c('BPC','CNLP',
                 'MAX_CN')
for(i in 1:length(Feature_names)){
  CurrentFeature<-Feature_names[i]
  for(j in 1:POFE[i]){
    CompNames<-c(CompNames,paste(CurrentFeature,j))
  }
}
colnames(SampleByComponent)<-CompNames
rownames(SampleByComponent)<-unique(Final_MFM[,1])
heatmap(SampleByComponent,Rowv = NA,Colv = NA)

## Set the parameters for PMSignature (Base matrix)
MFM_combinations <- apply(Final_MFM[,2:(NF+1)], 1, paste0, collapse='-')
num_combinations <- table(MFM_combinations)
featureVectorList_mod <- sapply(names(num_combinations), function(j){
  as.numeric(strsplit(j, split = "-")[[1]])})

countData_mod0 <- apply(Final_MFM[,2:(NF+1)], 1, function(j) which(names(num_combinations) == paste0(j, collapse = '-')))
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

# ## Set the parameters for the PMSignature (randomized version)
# MFM_combinations_r <- apply(Random_Final_MFM[,2:(NF+1)], 1, paste0, collapse='-')
# num_combinations_r <- table(MFM_combinations_r)
# featureVectorList_mod_r <- sapply(names(num_combinations_r), function(j){
#   as.numeric(strsplit(j, split = "-")[[1]])})
# 
# countData_mod0_r <- apply(Random_Final_MFM[,2:(NF+1)], 1, function(j) which(names(num_combinations_r) == paste0(j, collapse = '-')))
# ## adding the numerical id of each patient
# countData_mod0_r <- cbind(patientid=countData_mod0_r,
#                         cat_combination=as.numeric(as.factor(Random_Final_MFM[,1])),
#                         count=1)
# ##' counting the number of instances in each category of combination of features
# ##' for each patient
# countData_mod_r <- aggregate(count~patientid+cat_combination, countData_mod0_r, sum)
# countData_Lena_r <- t(countData_mod_r)
# 
# 
# MFD_Lena_r<-new(Class="MutationFeatureData",featureVectorList=featureVectorList_mod_r,
#               sampleList=as.character(unique(Random_Final_MFM[,1])),countData=countData_Lena_r, 
#               possibleFeatures=as.integer(POFE))



## Select Number of Signatures
min_NS<-3
max_NS<-10
numcorsigs<-c()
meancorsigs<-c()
sparsigcomp<-c()
BIC_MFD<-c()
for(k in min_NS:max_NS){
  NOS<-k #Number of Signature
  MFD_og<-getPMSignature(MFD_Lena,K=NOS,numInit = 10)
  #MFD_r<-getPMSignature(MFD_Lena_r,K=NOS,numInit = 10)
  
  #Signature by Component 
  LMat<-list()
  NZrows<-c()
  NZcols<-c()
  #LMat_r<-list()
  #NZrows_r<-c()
  #NZcols_r<-c()
  for(i in 1:NOS){
    Mat<-getSignatureValue(MFD_og, i)
    #Mat_r<-getSignatureValue(MFD_r, i)
    LMat[[i]]<-Mat
    NZrows<-c(NZrows,which(Mat>0,arr.ind = TRUE)[,1])
    NZcols<-c(NZcols,which(Mat>0,arr.ind = TRUE)[,2])
    #LMat_r[[i]]<-Mat_r
    #NZrows_r<-c(NZrows_r,which(Mat_r>0,arr.ind = TRUE)[,1])
    #NZcols_r<-c(NZcols_r,which(Mat_r>0,arr.ind = TRUE)[,2])
  }
  #Selected components
  selcomp<-unique(matrix(c(NZrows,NZcols),ncol=2))
  #selcomp_r<-unique(matrix(c(NZrows_r,NZcols_r),ncol=2))
  #reorder in ascending
  selcomp<-selcomp[order(selcomp[,1],selcomp[,2]),]
  #selcomp_r<-selcomp_r[order(selcomp_r[,1],selcomp_r[,2]),]
  selcompNames<-c()
  Feature_names<-c('Bp count','MAX CN',
                   'Weighted Sement sum')
  for(i in 1:nrow(selcomp)){
    CurrentFeature<-Feature_names[selcomp[i,1]]
    selcompNames<-c(selcompNames,paste(CurrentFeature,selcomp[i,2]))
  }
  SignatureByComponent<-matrix(0,nrow=NOS,ncol=nrow(selcomp))
  #SignatureByComponent_r<-matrix(0,nrow=NOS,ncol=nrow(selcomp_r))
  for(i in 1:NOS){
    Mat<-getSignatureValue(MFD_og, i)
    #Mat_r<-getSignatureValue(MFD_r, i)
    for(j in 1:nrow(selcomp)){
      SignatureByComponent[i,j]<-Mat[selcomp[j,1],selcomp[j,2]]
      #SignatureByComponent_r[i,j]<-Mat_r[selcomp_r[j,1],selcomp_r[j,2]]
    }
  }
  colnames(SignatureByComponent)<-selcompNames
  #colnames(SignatureByComponent_r)<-selcompNames
  mrq<-rquery.cormat(t(SignatureByComponent),type='full',graph = FALSE)$r
  numcorsigs<-c(numcorsigs,length(which(abs(mrq[mrq<1])>=0.5))/2)
  meancorsigs<-c(meancorsigs,mean(abs(mrq[mrq<1])))
  sparsigcomp<-c(sparsigcomp,sparseness(SignatureByComponent))
  #sparsigcomp_r<-c(sparsigcomp_r,sparseness(SignatureByComponent_r))
  BIC_MFD<-c(BIC_MFD,-2*MFD_og@loglikelihood+k*(sum(POFE)-1)*log(as.numeric(nrow(hq_CN))))
}
par(mfrow=c(2,2))
plot(min_NS:max_NS,numcorsigs,col='blue',type='l',lwd=3,xlab='Number of Signatures',ylab='Number of high corrs',
     main='Number of high correlations (>0.5)\n between signature pairs')
plot(min_NS:max_NS,meancorsigs,col='blue',main='Mean correlation between signatures',type='l',lwd=3,xlab='Number of Signatures',ylab='Mean correlation')
plot(min_NS:max_NS,sparsigcomp,col='blue',xlab='Number of Signatures',ylab='Sparseness',
     main='Sparseness of Signature by Component matrix',type='l',lwd=3)
#lines(min_NS:max_NS,sparsigcomp_r,col='red',lwd=3)
plot(min_NS:max_NS,BIC_MFD,main='BIC of the algorithm',col='blue',type='l',lwd=3,xlab='Number of Signatures',ylab='BIC')

NOS<-4 #Number of Signature
MFD_Param_Lena<-getPMSignature(MFD_Lena,K=NOS,numInit = 50)

#Signature by sample 
SampleBySignature<-MFD_Param_Lena@sampleSignatureDistribution
rownames(SampleBySignature)<-unique(Final_MFM[,1])
heatmap(SampleBySignature)

#Signature by Component 
LMat<-list()
NZrows<-c()
NZcols<-c()
for(i in 1:NOS){
  Mat<-getSignatureValue(MFD_Param_Lena, i)
  LMat[[i]]<-Mat
  NZrows<-c(NZrows,which(Mat>0,arr.ind = TRUE)[,1])
  NZcols<-c(NZcols,which(Mat>0,arr.ind = TRUE)[,2])
}
#Selected components
selcomp<-unique(matrix(c(NZrows,NZcols),ncol=2))
#reorder in ascending
selcomp<-selcomp[order(selcomp[,1],selcomp[,2]),]
selcompNames<-c()
Feature_names<-c('BPC','CNLP','MaxCN')
for(i in 1:nrow(selcomp)){
  CurrentFeature<-Feature_names[selcomp[i,1]]
  selcompNames<-c(selcompNames,paste(CurrentFeature,selcomp[i,2]))
}
SignatureByComponent<-matrix(0,nrow=NOS,ncol=nrow(selcomp))
for(i in 1:NOS){
  Mat<-getSignatureValue(MFD_Param_Lena, i)
  for(j in 1:nrow(selcomp)){
    SignatureByComponent[i,j]<-Mat[selcomp[j,1],selcomp[j,2]]
  }
}
colnames(SignatureByComponent)<-selcompNames
heatmap(SignatureByComponent,Rowv = NA,Colv = NA,scale='none')
visMembership(MFD_Lena, MFD_Param_Lena,colourBrewer = 'Paired')

a1<-visMembership(MFD_Lena, MFD_Param_Lena,colourBrewer = 'Paired')
a2<-visMembership(MFD_Lena_SE, MFD_Param_Lena_SE,colourBrewer = 'Paired')



TSigs<-7+4+4

CombMat<-matrix(0,nrow=length(unique(Final_MFM[,1])),ncol=TSigs+2)
samp_names<-unique(Final_MFM[,1])
GeoffSig<-t(coef(sigs))
for(i in 1:dim(GeoffSig)[1]){
  GeoffSig[i,]<-GeoffSig[i,]/sum(GeoffSig[i,])
}
CombMat[,1:7]<-GeoffSig
CombMat[,8:11]<-MFD_Param_Lena@sampleSignatureDistribution
CombMat[,12:15]<-MFD_Param_Lena_SE@sampleSignatureDistribution
nm_a1_max<-max(table(Final_MFM[,1]))
nm_a2_max<-max(table(Final_MFM_SE[,1]))
for(j in 1:length(unique(Final_MFM[,1]))){
  pos_a1<-match(names(sort(table(Final_MFM[,1])))[j],unique(Final_MFM[,1]))
  nm_a1<-sort(table(Final_MFM[,1]))[j]
  pos_a2<-match(names(sort(table(Final_MFM_SE[,1])))[j],unique(Final_MFM_SE[,1]))
  nm_a2<-sort(table(Final_MFM_SE[,1]))[j]
  CombMat[pos_a1,TSigs+1]<-nm_a1/nm_a1_max
  CombMat[pos_a2,TSigs+2]<-nm_a2/nm_a2_max
}
rownames(CombMat)<-as.character(unique(Final_MFM[,1]))
colnames(CombMat)<-c(paste0('G',1:7),paste0('A',1:4),paste0('B',1:4),'NM_A1','NM_A2')
hmCM<-heatmap(CombMat,scale="none")
heatmap(cor(CombMat),scale='none' )

maxscors<-c()
for(i in 1:7){
  maxscors<-c(maxscors,sort(abs((cor(CombMat))[,i]))[16])
}
namesbpmc<-c()
for(j in 1:7){
  namesbpmc<-c(namesbpmc,names(maxscors)[j],paste0('G',j))
}
barplot(matrix(c(maxscors,rep(1,7)),nrow=2,byrow=TRUE),beside = TRUE,
        names.arg = namesbpmc,col=c('orange','red2'),las=2,
        main='Maximum correlation with Geoff signatures',
        ylab='Correlation (abs. value)')



hmCM$rowInd
mat3<-matrix(0,nrow=nrow(mat2),ncol=ncol(mat2))
mat3names<-c()
for(i in 1:length(hmCM$rowInd)){
  mat3[length(hmCM$rowInd)-i+1,]<-mat2[hmCM$rowInd[i],]
  mat3names<-c(mat3names,samps[hmCM$rowInd[i]])
}
rownames(mat3)<-rev(mat3names)
pheatmap(mat3,show_rownames=TRUE,show_colnames=FALSE,
         scale = "none",cluster_rows = FALSE,clustering_method="ward.D2",
         clustering_distance_cols="euclidean")

##Mutations information
mut_info<-readRDS("data/mutation_matrix.RDS")
mut_info<-mut_info[which(rownames(mut_info)%in%unique(Final_MFM[,1])),]
mut_info<-mut_info[,which(colSums(1*as.matrix(mut_info))>3)]
SBS<-MFD_Param_Lena@sampleSignatureDistribution
rownames(SBS)<-unique(Final_MFM[,1])
mut_mat<-matrix(0,nrow=nrow(mut_info),ncol=ncol(mut_info)+ncol(SBS))
mut_mat[,1:ncol(mut_info)]<-1*as.matrix(mut_info)
mut_mat[,(ncol(mut_info)+1):ncol(mut_mat)]<-SBS[match(rownames(mut_info),rownames(SBS)),]
rownames(mut_mat)<-rownames(mut_info)
colnames(mut_mat)<-c(colnames(mut_info),'S1','S2','S3','S4')
heatmap(mut_mat,scale='none')


