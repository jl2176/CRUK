## Mutational Hotspots
## Load all_CN from Geoff
samps<-getSampNames(all_CN)

par(mfrow=c(2,2))
for(c in 5:8){
  CurrentChrLen<-chrlen[chrlen[,1]==paste0("chr",c),2]
  hotspots_gains<-rep(0,round(CurrentChrLen/10000))
  hotspots_gains_weighted<-rep(0,round(CurrentChrLen/10000))
  hotspots_losses<-rep(0,round(CurrentChrLen/10000))
  hotspots_losses_weighted<-rep(0,round(CurrentChrLen/10000))
  for(i in samps){
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
    for(k in 1:nrow(currseg)){
      if(currseg$segVal[k]>2){
        hotspots_gains[currseg$start[k]:(currseg$end[k]-1)]<-
          hotspots_gains[currseg$start[k]:(currseg$end[k]-1)]+1
        hotspots_gains_weighted[currseg$start[k]:(currseg$end[k]-1)]<-
          hotspots_gains_weighted[currseg$start[k]:(currseg$end[k]-1)]+currseg$segVal[k]-2
      }
      if(currseg$segVal[k]<2){
        hotspots_losses[currseg$start[k]:(currseg$end[k]-1)]<-
          hotspots_losses[currseg$start[k]:(currseg$end[k]-1)]+1
        hotspots_losses_weighted[currseg$start[k]:(currseg$end[k]-1)]<-
          hotspots_losses_weighted[currseg$start[k]:(currseg$end[k]-1)]+2-currseg$segVal[k]
      }
    }
  }
  plot(hotspots_gains_weighted,col='blue',main=paste('Chr.',c),xlab='Position',
       ylab='Counts')
  points(hotspots_gains,col='lightblue')
  points(hotspots_losses_weighted,col='red')
  points(hotspots_losses,col='orange')
}

