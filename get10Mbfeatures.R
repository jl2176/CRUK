##DESCRIPTION

# This function takes as input a QDNAseq file, describing the CN information
# for all the positions in a genome and transforms it into a data frame with
# information about CN features grouped by 10Mb segments.

# The selected features for each 10Mb segment are:
# Break point counts, Copy number loss percentage, Copy number gain percentage,
# Sum of CN change points (in absolute value), Maximum CN value, 
# Minimum CN value and Weighted segment sum. 

get10Mbfeatures<-function(abs_profiles,chrlen,keeptrack)
{
  samps<-getSampNames(abs_profiles)
  ## Features of segments
  Col0<-c() #sample name
  Col1<-c() #chromosome name
  Col2<-c() #break point count
  Col3<-c() #copy loss % 
  Col4<-c() #copy gain % 
  Col5<-c() #sum CN change points
  Col6<-c() #Maximum CN
  Col7<-c() #Minimum CN
  Col8<-c() #Weighted segment sum 
  for(i in samps)
  {
    if(class(abs_profiles)=="QDNAseqCopyNumbers")
    {
      segTab<-getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
    }else
    {
      segTab<-abs_profiles[[i]]
      colnames(segTab)[4]<-"segVal"
    }
    if(keeptrack==1){
      print(i)
    }
    # Round the segment values to integers
    segTab$segVal[as.numeric(segTab$segVal)<0]<-0
    segTab$segVal<-round(as.numeric(segTab$segVal))
    chrs<-unique(segTab$chromosome)
    for(c in chrs)
    {
      currseg<-segTab[segTab$chromosome==c,]
      # Add information for missing gaps
      k=1
      while(k<dim(currseg)[1]){
        if((as.numeric(currseg$start[k+1])-as.numeric(currseg$end[k]))!=1){
          newcurrseg<-currseg[1:k,]
          rowtoadd<-data.frame(c,as.numeric(currseg$end[k])+1,
                               as.numeric(currseg$start[k+1])-1,currseg$segVal[k])
          names(rowtoadd)<-colnames(currseg)
          newcurrseg[k+1,]<-rowtoadd
          newcurrseg[(k+2):(dim(currseg)[1]+1),]<-currseg[(k+1):dim(currseg)[1],]
          currseg<-newcurrseg
        }
        k<-k+1
      }
      # Define intervals and create a new data frame with the additional 
      # breakpoints of the intervals
      CurrentChrLen<-chrlen[chrlen[,1]==paste0("chr",c),2]
      intervals<-seq(1,(floor(CurrentChrLen/1e7)+1)*1e7,1e7)
      NewStart<-unique(sort(c(as.numeric(currseg$start),intervals)))
      NewStart<-NewStart[NewStart<intervals[length(intervals)]]
      NewEnd<-sort(c(as.numeric(currseg$end),intervals[-1]-1,
                     as.numeric(currseg$start)[1]-1))
      NewEnd<-NewEnd[NewEnd<intervals[length(intervals)]]
      NewEnd<-unique(NewEnd[NewEnd>0])
      upperlimit<-NewEnd[length(NewEnd)]
      currseg<-currseg[1:min(which(as.numeric(currseg$end)>upperlimit)),]
      SVa<-as.numeric(currseg$segVal)
      NewsegVal<-c()
      if(as.numeric(currseg$start[1])>1){
        NewsegVal<-c(NewsegVal,SVa[1])
      }
      previousEndPos<-0
      for(j in 1:(length(intervals)-1)){
        EndPos<-min(which(((intervals[j+1]-as.numeric(currseg$end))<=1)==TRUE))
        if(EndPos>previousEndPos){
          NewsegVal<-c(NewsegVal,SVa[(previousEndPos+1):EndPos])
          if(EndPos<length(SVa)&&(min(abs((intervals[j+1]-as.numeric(currseg$end)))-1)!=0)){
            NewsegVal<-c(NewsegVal,SVa[EndPos])
          }
        }
        if(EndPos==previousEndPos&&(min(abs((intervals[j+1]-as.numeric(currseg$end)))-1)!=0)){
          NewsegVal<-c(NewsegVal,SVa[EndPos])
        }
        previousEndPos<-EndPos
      }
      # Newcurrseg defines the new data frame with updated intervals and an 
      # additional column for the size of the intervals. 
      Newcurrseg<-data.frame(rep(c,length(NewStart)),NewStart,NewEnd,
                             NewEnd-NewStart+1,NewsegVal)
      colnames(Newcurrseg)<-c('chromosome','start','end','size','segVal')
      segment_start<-0
      # Now we can calculate the information for the features
      for(j in 1:(length(intervals)-1)){
        segment_end<-which(Newcurrseg$end==(intervals[j+1]-1))
        segment_rows<-Newcurrseg[(segment_start+1):segment_end,]
        segment_start<-segment_end
        Col0<-c(Col0,i) #sample name
        Col1<-c(Col1,c) #current chromosome
        Col2<-c(Col2,length(which(diff(segment_rows$segVal)!=0))) #break points
        segment_losses<-which(segment_rows$segVal<2)
        segment_gains<-which(segment_rows$segVal>2)
        #Col3 = copy loss % 
        Col3<-c(Col3,sum(segment_rows$size[segment_losses])/sum(segment_rows$size))
        #Col4 = copy gain %
        Col4<-c(Col4,sum(segment_rows$size[segment_gains])/sum(segment_rows$size))
        Col5<-c(Col5,sum(abs(diff(segment_rows$segVal)))) #sum CN changes
        Col6<-c(Col6,max(segment_rows$segVal)) #max CN 
        Col7<-c(Col7,min(segment_rows$segVal)) #min CN
        #Col8 = (weighted segment sum) 
        Col8<-c(Col8,sum(segment_rows$size*segment_rows$segVal))
      }
    }
  }
  # Group the percentages and the weighted sum in bins
  Col3<-floor(Col3*20)
  Col4<-floor(Col4*20)
  Col8<-floor(Col8/5e7)
  MFM<-data.frame(Col0,Col1,Col2,Col3,Col4,Col5,Col6,Col7,Col8)
  return(MFM)
}
