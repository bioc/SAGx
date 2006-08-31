# From Venables and Ripley 
myclus<-function(data = swiss,k = 3){
hc <- hclust(dist(data), "ave")
 memb <- cutree(hc, k)
 cent <- NULL
for(i in 1:k){
   cent <- rbind(cent, apply(data[memb==i,,drop=FALSE], 2, mean))
 }
return(kmeans(data,cent))
}
