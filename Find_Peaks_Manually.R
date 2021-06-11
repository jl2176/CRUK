c=22
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
which(hotspots_gains_weighted==max(hotspots_gains_weighted))
which(hotspots_gains_weighted[12000:14000]==max(hotspots_gains_weighted[12000:14000]))+12000
which(hotspots_losses_weighted==max(hotspots_losses_weighted))
which(hotspots_losses_weighted[9000:12000]==max(hotspots_losses_weighted[9000:12000]))+9000

