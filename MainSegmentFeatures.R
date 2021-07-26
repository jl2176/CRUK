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

dat<-as.numeric(CN_SegmentFeatures[["CNval"]][,2])
CNval_mm<-fitComponent(dat,seed=seed,model_selection=model_selection,
                         min_prior=min_prior,niter=niter,nrep=nrep,min_comp=2,max_comp=6)
MF2$CNval[which(MF2$CNval==0)]<-1
MF2$CNval[which(MF2$CNval>5)]<-6
MF2$CNval[which(MF2$CNval>1)]<-MF2$CNval[which(MF2$CNval>1)]-1

dat<-as.numeric(CN_SegmentFeatures[["SegSize"]][,2])
SegSize_mm<-fitComponent(dat,seed=seed,model_selection=model_selection,
                         min_prior=min_prior,niter=niter,nrep=nrep,min_comp=2,max_comp=10)

Col3_new<-MF2$SegSize
Col3_new[which(Col3_new<=5e5)]<-1
Col3_new[which(Col3_new>5e5 & Col3_new<=1e6)]<-2
Col3_new[which(Col3_new>1e6 & Col3_new<=5e6)]<-3
Col3_new[which(Col3_new>5e6 & Col3_new<=1e7)]<-4
Col3_new[which(Col3_new>1e7 & Col3_new<=5e7)]<-5
Col3_new[which(Col3_new>5e7)]<-6
MF2$SegSize<-Col3_new

dat<-as.numeric(CN_SegmentFeatures[["DistToCent"]][,2])
DistToCent_mm<-fitComponent(dat,seed=seed,model_selection=model_selection,
                      min_prior=min_prior,niter=niter,nrep=nrep,min_comp=2,max_comp=6)

Col9_new<-MF2$DistToCent
Col9_new[which(Col9_new<=1e7)]<-1
Col9_new[which(Col9_new>1e7 & Col9_new<=5e7)]<-2
Col9_new[which(Col9_new>5e7 & Col9_new<=9e7)]<-3
Col9_new[which(Col9_new>9e7)]<-4
MF2$DistToCent<-Col9_new

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

# MF2$CNval[which(MF2$CNval<2)]<-1
# MF2$CNval[which(MF2$CNval>=2&MF2$CNval<=6)]<-2
# MF2$CNval[which(MF2$CNval>6)]<-3

# MF2$RelAdjCN[which(MF2$RelAdjCN<2)]<-1
# MF2$RelAdjCN[which(MF2$RelAdjCN>=2&MF2$RelAdjCN<6)]<-2
# MF2$RelAdjCN[which(MF2$RelAdjCN>=6)]<-3

MF2$RelAdjCN[which(MF2$RelAdjCN<(-3))]<-(-3)
MF2$RelAdjCN[which(MF2$RelAdjCN>3)]<-3
MF2$RelAdjCN<-MF2$RelAdjCN+4

Final_MFM_SE<-data.frame(MF2$ID,MF2$CNval,MF2$SegSize,MF2$DistToCent,MF2$RelAdjCN)

## Correlations between features
colnames(Final_MFM_SE)<-c('ID','CN','SegSize','DistCent','RelAdjCN')
heatmap(abs(cor(Final_MFM_SE[,-1])),scale='none')

NF<-4 #Number of features
POFE<-c()
for(j in 2:(NF+1)){
  POFE<-c(POFE,max(Final_MFM_SE[,j]))
}

#Sample by component matrix
SampleByComponent<-matrix(0,nrow=length(unique(Final_MFM_SE[,1])),ncol=sum(POFE))
cur_row<-0
for(i in unique(Final_MFM_SE[,1])){
  cur_row<-cur_row+1
  cur_pos<-1
  for(j in 1:NF){
    SampleByComponent[cur_row,cur_pos:(cur_pos+POFE[j]-1)]<-
      table(c(1:POFE[j],Final_MFM_SE[which(Final_MFM_SE[,1]==i),j+1]))-1
    cur_pos<-cur_pos+POFE[j]
  }
}
CompNames<-c()
Feature_names<-c('CNval','SegSize','DistToCent',
                 'RelAdjCN')
for(i in 1:length(Feature_names)){
  CurrentFeature<-Feature_names[i]
  for(j in 1:POFE[i]){
    CompNames<-c(CompNames,paste(CurrentFeature,j))
  }
}
colnames(SampleByComponent)<-CompNames
rownames(SampleByComponent)<-unique(Final_MFM_SE[,1])
heatmap(SampleByComponent,Rowv = NA,Colv = NA)


##Lena's alternative
MFM_combinations <- apply(Final_MFM_SE[,2:5], 1, paste0, collapse='-')
num_combinations <- table(MFM_combinations)
featureVectorList_mod <- sapply(names(num_combinations), function(j){
  as.numeric(strsplit(j, split = "-")[[1]])})

countData_mod0 <- apply(Final_MFM_SE[,2:5], 1, function(j) which(names(num_combinations) == paste0(j, collapse = '-')))
## adding the numerical id of each patient
countData_mod0 <- cbind(patientid=countData_mod0,
                        cat_combination=as.numeric(as.factor(Final_MFM_SE[,1])),
                        count=1)
##' counting the number of instances in each category of combination of features
##' for each patient
countData_mod <- aggregate(count~patientid+cat_combination, countData_mod0, sum)
countData_Lena <- t(countData_mod)


MFD_Lena_SE<-new(Class="MutationFeatureData",featureVectorList=featureVectorList_mod,
              sampleList=as.character(unique(Final_MFM_SE[,1])),countData=countData_Lena, 
              possibleFeatures=as.integer(POFE))

## Select Number of Signatures
min_NS<-3
max_NS<-12
numcorsigs<-c()
numcorsigs_6<-c()
numcorsigs_7<-c()
meancorsigs<-c()
sparsigcomp<-c()
BIC_MFD<-c()
for(k in min_NS:max_NS){
  NOS<-k #Number of Signature
  MFD_og<-getPMSignature(MFD_Lena_SE,K=NOS,numInit = 10)
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
  Feature_names<-c('CN value','Segment length','Distance to centromere',
                   'Relation adjacent CN')
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
  mrq<-cor(t(SignatureByComponent))
  numcorsigs<-c(numcorsigs,length(which(abs(mrq[mrq<1])>=0.5))/2)
  meancorsigs<-c(meancorsigs,mean(abs(mrq[mrq<1])))
  sparsigcomp<-c(sparsigcomp,sparseness(SignatureByComponent))
  numcorsigs_6<-c(numcorsigs_6,length(which(abs(mrq[mrq<1])>=0.6))/2)
  numcorsigs_7<-c(numcorsigs_7,length(which(abs(mrq[mrq<1])>=0.7))/2)
  #sparsigcomp_r<-c(sparsigcomp_r,sparseness(SignatureByComponent_r))
  BIC_MFD<-c(BIC_MFD,-2*MFD_og@loglikelihood+k*(sum(POFE)-1)*log(as.numeric(nrow(hq_CN))))
}
par(mfrow=c(2,2))
plot(min_NS:max_NS,meancorsigs,col='blue',main='Mean correlation between signatures',
     type='l',lwd=3,xlab='Number of Signatures',ylab='Mean correlation')
plot(min_NS:max_NS,numcorsigs,col='blue',type='l',lwd=3,xlab='Number of Signatures',ylab='Pairs of signatures',
     main='Number of correlations >= 0.5\n between signature pairs')
#plot(min_NS:max_NS,sparsigcomp,col='blue',xlab='Number of Signatures',ylab='Sparseness',
#      main='Sparseness of Signature by Component matrix',type='l',lwd=3)

plot(min_NS:max_NS,numcorsigs_6,col='blue',type='l',lwd=3,xlab='Number of Signatures',ylab='Pairs of signatures',
     main='Number of correlations >= 0.6\n between signature pairs')

plot(min_NS:max_NS,numcorsigs_7,col='blue',type='l',lwd=3,xlab='Number of Signatures',ylab='Pairs of signatures',
     main='Number of correlations >= 0.7\n between signature pairs')
#plot(min_NS:max_NS,BIC_MFD,main='BIC of the algorithm',col='blue',type='l',lwd=3,xlab='Number of Signatures',ylab='BIC')


## Plots for chosen number of signatures
NOS<-4 #Number of Signature
MFD_Param_Lena_SE<-getPMSignature(MFD_Lena_SE,NOS,numInit = 30)

#Signature by sample 
SampleBySignature_SE<-MFD_Param_Lena_SE@sampleSignatureDistribution
rownames(SampleBySignature_SE)<-unique(Final_MFM_SE[,1])
heatmap(SampleBySignature_SE)

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
Feature_names<-c('CNval','SegSize','DistToCent',
                 'RelAdjCN')
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
# cont<-1
# for(j in 1:NF){
#   heatmap(SignatureByComponent[,cont:(cont+POFE[j]-1)], Rowv = NA, Colv = NA,scale='none')
#   cont<-cont+POFE[j]
# }
visMembership(MFD_Lena_SE, MFD_Param_Lena_SE,colourBrewer = 'Paired')

TSigs<-7+NOS

CombMat_SE<-matrix(0,nrow=length(unique(Final_MFM_SE[,1])),ncol=TSigs+1)
samp_names<-unique(Final_MFM_SE[,1])
GeoffSig<-t(coef(sigs))
for(i in 1:dim(GeoffSig)[1]){
  GeoffSig[i,]<-GeoffSig[i,]/sum(GeoffSig[i,])
}
CombMat_SE[,1:7]<-GeoffSig
CombMat_SE[,8:TSigs]<-MFD_Param_Lena_SE@sampleSignatureDistribution
nm_a2_max<-max(table(Final_MFM_SE[,1]))
for(j in 1:length(unique(Final_MFM[,1]))){
  pos_a2<-match(names(sort(table(Final_MFM_SE[,1])))[j],unique(Final_MFM_SE[,1]))
  nm_a2<-sort(table(Final_MFM_SE[,1]))[j]
  CombMat_SE[pos_a2,TSigs+1]<-nm_a2/nm_a2_max
}
rownames(CombMat_SE)<-as.character(unique(Final_MFM_SE[,1]))
colnames(CombMat_SE)<-c(paste0('G',1:7),paste0('B',1:NOS),'NM_A2')
hmCM<-heatmap(CombMat_SE,scale="none")
heatmap(cor(CombMat_SE),scale='none' )

maxscors<-c()
for(i in 1:7){
  maxscors<-c(maxscors,sort(abs((cor(CombMat_SE))[8:(TSigs),i]))[NOS])
}
namesbpmc<-c()
for(j in 1:7){
  namesbpmc<-c(namesbpmc,names(maxscors)[j],paste0('G',j))
}
par(mfrow=c(1,1))
barplot(matrix(c(maxscors,rep(1,7)),nrow=2,byrow=TRUE),beside = TRUE,
        names.arg = namesbpmc,col=c('navyblue','lightblue'),las=2,
        main='Maximum correlation with Geoff signatures',
        ylab='Correlation (abs. value)')

maxscors<-c()
for(i in 1:7){
  maxscors<-c(maxscors,sort(abs((cor(CombMat_SE))[8:(TSigs+1),i]))[NOS+1])
}
namesbpmc<-c()
for(j in 1:7){
  namesbpmc<-c(namesbpmc,names(maxscors)[j],paste0('G',j))
}
par(mfrow=c(1,1))
barplot(matrix(c(maxscors,rep(1,7)),nrow=2,byrow=TRUE),beside = TRUE,
        names.arg = namesbpmc,col=c('navyblue','lightblue'),las=2,
        main='Maximum correlation with Geoff signatures',
        ylab='Correlation (abs. value)')



