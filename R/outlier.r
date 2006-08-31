# This program identifies outliers in a multivariate distribution #
outlier<-function(M){
   library(stats)
   tr.pca<-prcomp(M, retx=TRUE,center=TRUE,scale.=TRUE) # SVD DIRECTLY ON DATA MATRIX #
   x<-as.vector(tr.pca$sdev^2);y<-cumsum(x)/sum(x)
   sign<-max(length(z<-y[y<0.9]),1) 
# NUMBER OF SIGNIFICANT COMPONENTS (EIGEN VALUES GREATER THAN 1) #
   print(paste("The number of significant dimensions is",sign))
   vars<-as.vector(tr.pca$sdev^-2)[1:sign]  # INVERSE OF COMPONENTS COVARIANCE MATRIX #                    
   sc<-as.matrix(tr.pca$x)[,1:sign] # $scores from significant factors #
   load<-as.matrix(tr.pca$rotation)[,1:sign] # loadings #
   dmodx<-(M -as.matrix(sc) %*% t(as.matrix(load)))^2 %*% (numeric(length(x))+1)
# CALCULATES THE Q STATISTIC #
   nt<-sign+1
   hot2 <-as.matrix(sc)^2 %*% as.matrix(vars) 
   ut<-data.frame(cbind(hot2,dmodx))
   colnames(ut)<-c("Hotelling","DMODX")
return(ut)}
