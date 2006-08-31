# Copyright Per Broberg #
Zfreq<-function(z,Z,B=Bs,outz=TRUE){
library(stats)
limes<-max(abs(cbind(z,Z)))
n1<- hist(z, breaks=seq(-1*limes,limes,by=2*limes/20), plot=FALSE)
n2<- hist(Z, breaks=seq(-1*limes,limes,by=2*limes/20), plot=FALSE)
pi<-n2$counts/(n1$counts+n2$counts)
pi[1:10]<-rev(pava(rev(pi[1:10])))
pi[11:20]<-pava(pi[11:20])
ip<-cbind(n1$mids,pi);ip<-na.omit(ip[(ip[,"pi"] >0)&(ip[,"pi"] <1),])
logit<-log(ip[,"pi"]/(1-ip[,"pi"]));logit
t<-spline(ip[,1],logit,n=201)
pihat<-exp(t$y)/(exp(t$y)+1)
p0<-min((1/pihat-1)^(-1)*B);p1<-1-p0
library(boot) # Here is a bug! #
lgt<-1-p0*(1/pihat-1)/B
lgt[(lgt<0)]<-0
plt<-cbind(t$x,logit(lgt)) 
plt<-plt[!(is.na(plt[,2])|is.infinite(plt[,2])),]
plt.lo<-loess(plt[,2]~plt[,1])
plot(plt.lo,ylab="logit(p1(Z))",xlab="Z",type="n",ylim=c(0,2))
lines(plt.lo,col="red")
plt.lo<-loess(plt[,2]~plt[,1])
if(outz==TRUE){temp<-predict.loess(plt.lo,Z);return(data.frame(p1<-exp(temp)/(1+exp(temp))))}
}
