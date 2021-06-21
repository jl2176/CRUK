##DESCRIPTION

# This function takes as input a QDNAseq file, describing the CN information
# for all the positions in a genome and transforms it into a data frame with
# information about CN features grouped by NMb segments.

# The selected features for each NMb segment are:
# Break point counts, Copy number loss percentage, Copy number gain percentage,
# Sum of CN change points (in absolute value), Maximum CN value, 
# Minimum CN value, Weighted segment sum and Distance to closest chromosome end.

getNMbfeatures<-function(abs_profiles,chrlen,centromere_pos,N,keeptrack=0)
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
  Col9<-c() #Distance to centromere
  Col10<-c() #Length Oscilating CN
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
    segTab$segVal[as.numeric(segTab$segVal) < 0] <- 0
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
      CurrentCentPos<-centromere_pos[centromere_pos[,1]==c,4]
      intervals<-seq(1,(floor(CurrentChrLen/(N*10^6))+1)*N*10^6,N*10^6)
      NewStart<-unique(sort(c(as.numeric(currseg$start),intervals)))
      NewStart<-NewStart[NewStart<intervals[length(intervals)]]
      NewEnd<-sort(c(as.numeric(currseg$end),intervals[-1]-1,
                     as.numeric(currseg$start)[1]-1))
      NewEnd<-NewEnd[NewEnd<intervals[length(intervals)]]
      NewEnd<-unique(NewEnd[NewEnd>0])
      upperlimit<-NewEnd[length(NewEnd)]
      if(sum(as.numeric(currseg$end)>=upperlimit)==0){
        NewEnd<-NewEnd[-(length(NewEnd)-1)]
        currseg$end[dim(currseg)[1]]<-upperlimit
      }
      currseg<-currseg[1:min(which(as.numeric(currseg$end)>=upperlimit)),]
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
        if((EndPos==previousEndPos)&&
           ((min(abs((intervals[j+1]-as.numeric(currseg$end)))-1)!=0)||(j==(length(intervals)-1)))){
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
        segment_centre<-Newcurrseg[segment_end,2]-floor(N*10^6/2)
        distCent<-min(abs(segment_start-CurrentCentPos),abs(segment_end-CurrentCentPos))
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
        #Col9 = (distance to centromere)/(chromosome length)
        Col9<-c(Col9,distCent)
        #Number of oscillating CN segments
        oscCounts<-0
        if(nrow(segment_rows)>3)
        {
          prevval<-segment_rows$segVal[1]
          count=0
          for(j in 3:nrow(segment_rows))
          {
            if(segment_rows$segVal[j]==prevval&segment_rows$segVal[j]!=segment_rows$segVal[j-1])
            {
              count<-count+1
            }else{
              oscCounts<-oscCounts+count
              count=0
            }
            prevval<-segment_rows$segVal[j-1]
          }
        }
        Col10<-c(Col10,oscCounts)
      }
    }
  }
  MFM<-data.frame(Col0,Col1,Col2,Col3,Col4,Col5,Col6,Col7,Col8,Col9,Col10)
  return(MFM)
}