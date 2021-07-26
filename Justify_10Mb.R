MF5_j<-getNMbfeatures(hq_CN,chrlen,C19,5,0)
MF8_j<-getNMbfeatures(hq_CN,chrlen,C19,8,0)
MF10_j<-getNMbfeatures(hq_CN,chrlen,C19,10,0)
MF15_j<-getNMbfeatures(hq_CN,chrlen,C19,15,0)
MF20_j<-getNMbfeatures(hq_CN,chrlen,C19,20,0)
MF25_j<-getNMbfeatures(hq_CN,chrlen,C19,25,0)

#Breakpoint count
par(mfrow=c(2,2))
mean_bp_count<-c(mean(MF5_j[,3]),mean(MF8_j[,3]),mean(MF10_j[,3]),
                 mean(MF15_j[,3]),mean(MF20_j[,3]),mean(MF25_j[,3]))
plot(1:6,mean_bp_count,col='red',xlab='Segment size (Mb)',
     main='Mean number of breakpoints\n per segment',ylab='Mean number of BP',
     xaxt='n',lwd=3)
axis(1,at=1:6,labels=c(5,8,10,15,20,25))

#Segments per sample
segs_per_exp<-c(dim(MF5_j)[1]/91,dim(MF8_j)[1]/91,dim(MF10_j)[1]/91,
                dim(MF15_j)[1]/91,dim(MF20_j)[1]/91,dim(MF25_j)[1]/91)

plot(1:6,segs_per_exp,col='red',xlab='Segment size (Mb)',
     main='Number of segments per sample',ylab='Number of segments',
     xaxt='n',lwd=3)
axis(1,at=1:6,labels=c(5,8,10,15,20,25))

#Information lost
inf_lost<-c(sum(chrlen[-c(8,21),2]%%5e6),sum(chrlen[-c(8,21),2]%%8e6),sum(chrlen[-c(8,21),2]%%10e6),
            sum(chrlen[-c(8,21),2]%%15e6),sum(chrlen[-c(8,21),2]%%20e6),sum(chrlen[-c(8,21),2]%%25e6))
plot(1:6,inf_lost/1e6,col='red',xlab='Segment size (Mb)',
     main='Total information lost \n per sample (Mb)',ylab='Information lost (Mb)',
     xaxt='n',lwd=3)
axis(1,at=1:6,labels=c(5,8,10,15,20,25))

#Unique/Total segments
uniq_tot<-c(dim(distinct(MF5_j[,-c(1,2,10)]))[1]/dim(MF5_j)[1],dim(distinct(MF8_j[,-c(1,2,10)]))[1]/dim(MF8_j)[1],
            dim(distinct(MF10_j[,-c(1,2,10)]))[1]/dim(MF10_j)[1],dim(distinct(MF15_j[,-c(1,2,10)]))[1]/dim(MF15_j)[1],
            dim(distinct(MF20_j[,-c(1,2,10)]))[1]/dim(MF20_j)[1],dim(distinct(MF25_j[,-c(1,2,10)]))[1]/dim(MF25_j)[1])
plot(1:6,uniq_tot,col='red',xlab='Segment size (Mb)',
     main='Ratio of unique/total segments',ylab='Ratio',
     xaxt='n',lwd=3)
axis(1,at=1:6,labels=c(5,8,10,15,20,25))


#EXTRA: FEAT DIST EXAMPLE
par(mfrow=c(1,1))
hist(c(rnorm(1000,10,3),rnorm(1000,0,3),rnorm(1000,20,2)),col='lightblue',
     xlab='Value',main='Histogram for the feature distribution',xlim=c(-10,25),
     breaks=30)
x <- seq(-10, 25, length=100)
hx1 <- dnorm(x,10,3)
hx2 <- dnorm(x,0,3)
hx3 <- dnorm(x,20,2)
plot(x,hx1,col='blue2',lwd=3,type='l',ylim=c(0,0.2),ylab='Probability',xlab='Value',main='Result after Mixture Models')
lines(x,hx2,col='red3',lwd=3)
lines(x,hx3,col='orange',lwd=3)







