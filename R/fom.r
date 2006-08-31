
# FOM =Figure of Merit, a measure that decreases with improved clustering #
# Copyright Per Broberg 20APR01 #
fom<-function(data,cluster){
temp<-numeric(ncol(data))
for (i in 1:ncol(data)){
   temp[i]<-sqrt(sum(na.omit(aggregate(data[,-i],list(cluster),ssqvar<-function(x) var(x)*(length(x)-1))[,-1])))
}
return(sum(temp)/sqrt((nrow(data)-length(table(cluster)))/nrow(data)))
}



