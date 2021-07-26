SegmentFeatures<-function(abs_profiles,chrlen,centromere_pos)
{
  samps<-getSampNames(abs_profiles)
  ## Features of segments
  Col0<-c() #sample name
  Col1<-c() #chromosome name
  Col2<-c() #CN
  Col3<-c() #segment length
  Col4<-c() #Relation to CN of previous segment
  Col5<-c() #Relation to CN of following segment
  Col6<-c() #normalized distance to centromere
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
    # Round the segment values to integers
    segTab$segVal[as.numeric(segTab$segVal) < 0] <- 0
    segTab$segVal<-round(as.numeric(segTab$segVal))
    chrs<-unique(segTab$chromosome)
    for(c in chrs)
    {
      CurrentChrLen<-chrlen[chrlen[,1]==paste0("chr",c),2]
      CurrentCentPos<-centromere_pos[centromere_pos[,1]==c,4]
      currseg<-segTab[segTab$chromosome==c,]
      curr_c<-currseg$chromosome
      curr_starts<-currseg$start
      curr_ends<-currseg$end
      curr_sv<-currseg$segVal
      dels_sv<-c()
      if(nrow(currseg)>1){
        for(qt in 2:length(curr_sv)){
          if(curr_sv[qt]==curr_sv[qt-1]){
            dels_sv<-c(dels_sv,qt)
          }
        }
      }
      if(length(dels_sv)>0){
        curr_c<-curr_c[-dels_sv]
        curr_starts<-curr_starts[-dels_sv]
        curr_ends<-curr_ends[-(dels_sv-1)]
        curr_sv<-curr_sv[-dels_sv]
        currseg_n<-colnames(currseg)
        currseg<-data.frame(curr_c,curr_starts,curr_ends,curr_sv)
        colnames(currseg)<-currseg_n
      }
      for(n in 1:nrow(currseg)){
        Col0<-c(Col0,i) #Sample name
        Col1<-c(Col1,c) #Chromosome name
        Col2<-c(Col2,currseg$segVal[n]) #CN of the segment
        Col3<-c(Col3,as.numeric(currseg$end[n])-as.numeric(currseg$start[n])) #Segment length
        #CN with respect to adjacent positions
        if(n==1&&nrow(currseg)>1){
          Col4<-c(Col4,0)
          Col5<-c(Col5,abs(currseg$segVal[n]-currseg$segVal[n+1]))
        }
        if(n>1&&n==nrow(currseg)){
          Col4<-c(Col4,abs(currseg$segVal[n]-currseg$segVal[n-1]))
          Col5<-c(Col5,0)
        }
        if(n>1&&n<nrow(currseg)){
          Col4<-c(Col4,currseg$segVal[n]-currseg$segVal[n-1])
          Col5<-c(Col5,currseg$segVal[n]-currseg$segVal[n+1])
        }
        if(n==1&&n==nrow(currseg)){
          Col4<-c(Col4,0)
          Col5<-c(Col5,0)
        }
        #Distance to centromere
        segment_ref<-min(abs(as.numeric(currseg$end[n])-CurrentCentPos),
                         abs(as.numeric(currseg$start[n])-CurrentCentPos))
        Col6<-c(Col6,segment_ref)
      }
    }
  }
  MFM<-data.frame(Col0,Col1,Col2,Col3,Col4,Col5,Col6)
  return(MFM)
}
