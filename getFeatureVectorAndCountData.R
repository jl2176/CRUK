##DESCRIPTION
#Given a matrix in the format of the output of 'get10Mbfeatures', this function
#returns the matrices 'Feature_Vector' of dimensions: 
#(number of features) x (number of observed combinations of features).
#And also 'Count_Data', a matrix of dimensions 3 rows, each containing: 
#combination of feature id, sample id and number of repetitions.

#The input must be a dataframe and only include a first column with the sample 
#names followed by columns with the information of the features.

getFeatureVectorAndCountData<-function(SegmentFeatures){
  SampNames<-unique(SegmentFeatures[,1])
  Feature_Vector<-unique(SegmentFeatures[,2:dim(SegmentFeatures)[2]])
  FeatureDictionary<-do.call('paste', c(Feature_Vector, sep = '\r'))
  CD_1<-c()
  CD_2<-c()
  CD_3<-c()
  i<-0
  for(mutfeat in FeatureDictionary){
    i<-i+1
    j<-0
    for(SampID in SampNames){
      j<-j+1
      SampRows<-which(SegmentFeatures[,1]==SampID)
      #Combinations of Features
      SampCombs<-do.call('paste', c(SegmentFeatures[SampRows,2:dim(SegmentFeatures)[2]], sep = '\r'))
      NumReps<-sum(SampCombs==mutfeat)
      if(NumReps>0){
        CD_1<-c(CD_1,i)
        CD_2<-c(CD_2,j)
        CD_3<-c(CD_3,NumReps)
      }
    }
  }
  Count_Data<-matrix(c(CD_1,CD_2,CD_3),nrow=3,byrow=TRUE)
  return(list(t(Feature_Vector),Count_Data))
}
