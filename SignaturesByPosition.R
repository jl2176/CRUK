chr_pos<-c(1,1,3,4,4,5,8,8,9,11,12,13,14,16,16,17,18,19,19,20,22)
base_pos<-c(14912,15095,23476,16906,11272,18230,6851,3930,12948,
            13990,6888,2530,11446,2117,3040,8330,3790,7560,
            330,3030,3080,4870)
samps<-getSampNames(hq_CN)
SampleByPosition_matrix<-matrix(0,nrow=length(samps),ncol=length(chr_pos))

row_idx<-0
for(i in samps){
  row_idx<-row_idx+1
  for(j in 1:length(chr_pos)){
    c<-chr_pos[j]
    CurrentChrLen<-chrlen[chrlen[,1]==paste0("chr",c),2]
    if(class(abs_profiles)=="QDNAseqCopyNumbers")
    {
      segTab<-getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
    }else
    {
      segTab<-abs_profiles[[i]]
      colnames(segTab)[4]<-"segVal"
    }
    # Round the segment values to integers
    segTab$segVal[as.numeric(segTab$segVal) < 0] <- 0
    segTab$segVal<-round(as.numeric(segTab$segVal))
    currseg<-segTab[segTab$chromosome==c,]
    currseg$start<-round(as.numeric(currseg$start)/10000)
    currseg$end<-round(as.numeric(currseg$end)/10000)
    pos_int<-base_pos[j]
    CN_int<-currseg$segVal[max(which(currseg$start<=pos_int))]
    SampleByPosition_matrix[row_idx,j]<-CN_int
  }
}
SampleByPosition_matrix[which(SampleByPosition_matrix>2)]<-3
SampleByPosition_matrix[which(SampleByPosition_matrix<2)]<-1
##SampleByPosition acts as the Sample By Component matrix
POFE_alt<-rep(3,length(chr_pos))
MFM_combinations_alt <- apply(SampleByPosition_matrix, 1, paste0, collapse='-')
num_combinations_alt <- table(MFM_combinations_alt)
featureVectorList_mod_alt <- sapply(names(num_combinations_alt), function(j){
  as.numeric(strsplit(j, split = "-")[[1]])})

countData_mod0_alt <- apply(SampleByPosition_matrix, 1, function(j) which(names(num_combinations_alt) == paste0(j, collapse = '-')))
## adding the numerical id of each patient
countData_mod0_alt <- cbind(patientid=countData_mod0_alt,
                        cat_combination=as.numeric(as.factor(samps)),
                        count=1)
##' counting the number of instances in each category of combination of features
##' for each patient
countData_mod_alt <- aggregate(count~patientid+cat_combination, countData_mod0_alt, sum)
countData_Lena_alt <- t(countData_mod_alt)

MFD_alt<-new(Class="MutationFeatureData",featureVectorList=featureVectorList_mod_alt,
              sampleList=as.character(samps),countData=countData_Lena_alt, 
              possibleFeatures=as.integer(POFE_alt))
NOS<-3 #Number of Signature
MFD_Param_Lena<-getPMSignature(MFD_alt,NOS)

heatmap(MFD_Param_Lena@sampleSignatureDistribution)

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
SignatureByComponent<-matrix(0,nrow=NOS,ncol=nrow(selcomp))
for(i in 1:NOS){
  Mat<-getSignatureValue(MFD_Param_Lena, i)
  for(j in 1:nrow(selcomp)){
    SignatureByComponent[i,j]<-Mat[selcomp[j,1],selcomp[j,2]]
  }
}
heatmap(SignatureByComponent,Rowv = NA,Colv = NA)
visMembership(MFD_Lena, MFD_Param_Lena)


##Feature map
install.packages('pheatmap')
library(pheatmap)

mat2<-SampleByPosition_matrix

pheatmap(mat2,show_rownames=FALSE,show_colnames=FALSE,
         scale = "none",clustering_method="ward.D2",
         clustering_distance_cols="euclidean")
rownames(mat2)<-samps
