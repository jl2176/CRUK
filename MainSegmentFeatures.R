MF20<-SegmentFeatures(all_CN,chrlen,C19)
MF20_hq<-SegmentFeatures(hq_CN,chrlen,C19)
#MF2<-MF20[,-2]
MF2<-MF20_hq[,-2]
MF2<-MF2[which(MF2$Col2!=2),]
# MF2$Col4[which(MF2$Col4<(-1))]<--1
# MF2$Col4[which(MF2$Col4>1)]<-1
# MF2$Col5[which(MF2$Col5<(-1))]<--1
# MF2$Col5[which(MF2$Col5>1)]<-1
MF2$Col7<-MF2$Col4+MF2$Col5
# MF2$Col7<-MF2$Col7+2
MF2<-MF2[,c(1,2,3,6,7)]
Feature_short_names<-c('ID','CNval','SegSize','DistToCent',
                 'RelAdjCN')
colnames(MF2)<-Feature_short_names
CN_SegmentFeatures<-list()
for(j in 2:ncol(MF2)){
  CN_SegmentFeatures[[j-1]]<-data.frame(MF2[,c(1,j)])
}
names(CN_SegmentFeatures)<-Feature_short_names[2:length(Feature_short_names)]
## Apply mixture models
seed=77777
min_prior=0.001
model_selection="BIC"
nrep=1
niter=1000

# dat<-as.numeric(CN_SegmentFeatures[["CNval"]][,2])
# CNval_mm<-fitComponent(dat,seed=seed,dist="pois",model_selection=model_selection,
#                          min_prior=min_prior,niter=niter,nrep=nrep,min_comp=2,max_comp=5)

dat<-as.numeric(CN_SegmentFeatures[["SegSize"]][,2])
SegSize_mm<-fitComponent(dat,seed=seed,dist="pois",model_selection=model_selection,
                         min_prior=min_prior,niter=niter,nrep=nrep,min_comp=2,max_comp=3)

SegSize_com<-parameters(SegSize_mm)
SegSize_sortcom<-sort(SegSize_com)
SegSize_finalcluster<-SegSize_mm@cluster
for(k in 1:length(SegSize_com)){
  mk<-match(SegSize_sortcom[k],SegSize_com)
  SegSize_finalcluster[which(SegSize_finalcluster==mk)]<-k+length(SegSize_com)
}
SegSize_finalcluster<-SegSize_finalcluster-length(SegSize_com)

dat<-as.numeric(CN_SegmentFeatures[["DistToCent"]][,2])
DistToCent_mm<-fitComponent(dat,dist="pois",seed=seed,model_selection=model_selection,
                      min_prior=min_prior,niter=niter,nrep=nrep,min_comp=2,max_comp=3)

DistToCent_com<-parameters(DistToCent_mm)
DistToCent_sortcom<-sort(DistToCent_com)
DistToCent_finalcluster<-DistToCent_mm@cluster
for(k in 1:length(DistToCent_com)){
  mk<-match(DistToCent_sortcom[k],DistToCent_com)
  DistToCent_finalcluster[which(DistToCent_finalcluster==mk)]<-k+length(DistToCent_com)
}
DistToCent_finalcluster<-DistToCent_finalcluster-length(DistToCent_com)


# dat<-as.numeric(CN_SegmentFeatures[["RelAdjCN"]][,2])
# RelAdjCN_mm<-fitComponent(dat,seed=seed,dist="pois",model_selection=model_selection,
#                           min_prior=min_prior,niter=niter,nrep=nrep,min_comp=2,max_comp=3)
# 
# RelAdjCN_com<-parameters(RelAdjCN_mm)
# RelAdjCN_sortcom<-sort(RelAdjCN_com)
# RelAdjCN_finalcluster<-RelAdjCN_mm@cluster
# for(k in 1:length(RelAdjCN_com)){
#   mk<-match(RelAdjCN_sortcom[k],RelAdjCN_com)
#   RelAdjCN_finalcluster[which(RelAdjCN_finalcluster==mk)]<-k+length(RelAdjCN_com)
# }
# RelAdjCN_finalcluster<-RelAdjCN_finalcluster-length(RelAdjCN_com)

MF2$CNval[which(MF2$CNval<2)]<-1
MF2$CNval[which(MF2$CNval>=2&MF2$CNval<=6)]<-2
MF2$CNval[which(MF2$CNval>6)]<-3

MF2$RelAdjCN[which(MF2$RelAdjCN<2)]<-1
MF2$RelAdjCN[which(MF2$RelAdjCN>=2&MF2$RelAdjCN<6)]<-2
MF2$RelAdjCN[which(MF2$RelAdjCN>=6)]<-3

Final_MFM<-data.frame(MF2$ID,MF2$CNval,SegSize_finalcluster,DistToCent_finalcluster,MF2$RelAdjCN)

POFE<-c()
for(j in 2:5){
  POFE<-c(POFE,max(Final_MFM[,j]))
}

##Lena's alternative
MFM_combinations <- apply(Final_MFM[,2:5], 1, paste0, collapse='-')
num_combinations <- table(MFM_combinations)
featureVectorList_mod <- sapply(names(num_combinations), function(j){
  as.numeric(strsplit(j, split = "-")[[1]])})

countData_mod0 <- apply(Final_MFM[,2:5], 1, function(j) which(names(num_combinations) == paste0(j, collapse = '-')))
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

## Select Number of Signatures
min_NS<-3
max_NS<-8
numcorsigs<-c()
meancorsigs<-c()
maxcorsigs<-c()
for(k in min_NS:max_NS){
  NOS<-k #Number of Signature
  MFD_Param_Lena_SE<-getPMSignature(MFD_Lena,NOS,numInit = 10)
  #Signature by Component 
  LMat<-list()
  NZrows<-c()
  NZcols<-c()
  for(i in 1:NOS){
    Mat<-getSignatureValue(MFD_Param_Lena_SE, i)
    #Mat[which(Mat<0.01)]<-0
    LMat[[i]]<-Mat
    NZrows<-c(NZrows,which(Mat>0,arr.ind = TRUE)[,1])
    NZcols<-c(NZcols,which(Mat>0,arr.ind = TRUE)[,2])
  }
  #Selected components
  selcomp <- unique(matrix(c(NZrows,NZcols),ncol=2))
  #reorder in ascending
  selcomp <- selcomp[order(selcomp[,1],selcomp[,2]),]
  
  selcompNames<-c()
  Feature_names<-c('CN value','Segment length','Distance to centromere',
                   'Relation adjacent CN')
  for(i in 1:nrow(selcomp)){
    CurrentFeature<-Feature_names[selcomp[i,1]]
    selcompNames<-c(selcompNames,paste(CurrentFeature,selcomp[i,2]))
  }
  SignatureByComponent<-matrix(0,nrow=NOS,ncol=nrow(selcomp))
  for(i in 1:NOS){
    Mat<-getSignatureValue(MFD_Param_Lena_SE, i)
    for(j in 1:nrow(selcomp)){
      SignatureByComponent[i,j]<-Mat[selcomp[j,1],selcomp[j,2]]
    }
  }
  colnames(SignatureByComponent) <- selcompNames
  mrq<-rquery.cormat(t(SignatureByComponent),type='full')$r
  numcorsigs<-c(numcorsigs,length(which(mrq[mrq<1]>=0.5)))
  meancorsigs<-c(meancorsigs,mean(mrq[mrq<1]))
  maxcorsigs<-c(maxcorsigs,max(mrq[mrq<1]))
}
par(mfrow=c(2,2))
plot(min_NS:max_NS,numcorsigs,col='blue',pch=19,xlab='Number of Signatures',ylab='Number of high corrs')
plot(min_NS:max_NS,meancorsigs,col='blue',pch=19,xlab='Number of Signatures',ylab='Mean correlation')
plot(min_NS:max_NS,maxcorsigs,col='blue',pch=19,xlab='Number of Signatures',ylab='Max correlation')

## Plots for chosen number of signatures
NOS<-5 #Number of Signature
MFD_Param_Lena_SE<-getPMSignature(MFD_Lena,NOS,numInit = 10)

#Signature by sample 
heatmap(MFD_Param_Lena_SE@sampleSignatureDistribution,scale = "none")

#Signature by Component 
LMat<-list()
NZrows<-c()
NZcols<-c()
for(i in 1:NOS){
  Mat<-getSignatureValue(MFD_Param_Lena_SE, i)
  #Mat[which(Mat<0.01)]<-0
  LMat[[i]]<-Mat
  NZrows<-c(NZrows,which(Mat>0,arr.ind = TRUE)[,1])
  NZcols<-c(NZcols,which(Mat>0,arr.ind = TRUE)[,2])
}
#Selected components
selcomp <- unique(matrix(c(NZrows,NZcols),ncol=2))
#reorder in ascending
selcomp <- selcomp[order(selcomp[,1],selcomp[,2]),]

selcompNames<-c()
Feature_names<-c('CN value','Segment length','Distance to centromere',
                 'Relation adjacent CN')
for(i in 1:nrow(selcomp)){
  CurrentFeature<-Feature_names[selcomp[i,1]]
  selcompNames<-c(selcompNames,paste(CurrentFeature,selcomp[i,2]))
}
selcompNames

SignatureByComponent<-matrix(0,nrow=NOS,ncol=nrow(selcomp))
for(i in 1:NOS){
  Mat<-getSignatureValue(MFD_Param_Lena_SE, i)
  for(j in 1:nrow(selcomp)){
    SignatureByComponent[i,j]<-Mat[selcomp[j,1],selcomp[j,2]]
  }
}
colnames(SignatureByComponent) <- selcompNames
heatmap(SignatureByComponent, Rowv = NA, Colv = NA,scale='none')
visMembership(MFD_Lena, MFD_Param_Lena_SE)

