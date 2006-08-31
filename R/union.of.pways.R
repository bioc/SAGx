# Author : PBG 09NOV05, based on function merge.list in package RCurl
# X = list(a=c(1,2),c=c(1,2));Y = list(a=c(3,4),d=c(12,2))

union.of.pways <- function(x,y){
# x = X;y = Y
   t1 = !is.list(x);t2 = !is.list(y);if(any(c(t1,t2))) stop("Arguments must be lists.")
   if (length(x) == 0) return(y)
    if (length(y) == 0) return(x)
    i = match(names(y), names(x))
    i = is.na(i)
    if (any(i)) x[names(y)[which(i)]] = y[which(i)]
    i = match(names(y), names(x))
    if (any(i)) for( j in i) x[[j]] = union(unlist(y[names(x)[j]]),unlist(x[[j]]))
    x
}
