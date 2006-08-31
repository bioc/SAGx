# A goodness of clustering measure based on the PRESS statistic #
# Copyright Per Broberg 20APR01 #
cluster.q<-function(data,cluster){
internal<-function(test){
adds<-0
for (i in 1:ncol(test)){
res<-0
   for (j in 1:nrow(test)){
        fmla <- as.formula(paste(colnames(test)[i],"~", paste(colnames(test)[-i], collapse= "+")))
   res<-(test[j,i]-predict(lm.test<-lm(fmla,
         data=test[-j,]),test[j,-i]))^2+res
}
q2<-1-res/((nrow(test)-1)*var(test[,i]));adds<-adds+q2
}
return(q2)
}
mean(by(data,cluster,internal))
}
