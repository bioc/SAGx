GSEA.mean.t <- function(data = X, samroc = samroc.res, probeset = probeset, pway = kegg, type = c("original","absolute"),  two.side = FALSE, cutoff = c(5,Inf), B = 1000, smooth = FALSE){
# Author : PBG 26OCT05, based on Tian et al. PNAS 2005.Modified 19MAY06. Use of design slot added 16OCT06, restandardisation 27OCT06
# data = X; samroc = samroc.res; probeset = probeset; pway = kegg; type = "original"; two.side = FALSE; cutoff = c(5,Inf); B = 1000
if(smooth) stop("Smoothing not implemented yet.")
type <- match.arg(type) 
if(type=="maxmean") stop("maxmean not implemented yet.") 

# if maxmean restand should be done separately for s- and s+

if(type=="maxmean") {
   maxplus <- function(x) pmax(x, 0)
   maxminus <- function(x) -pmin(x, 0)
   s.plus <- maxplus(statistic)
   s.minus <- max.minus(statistic)
   mean.plus.statistic <- mean(s.plus)
   mean.minus.statistic <- mean(s.minus)
   sd.plus.statistic <- sd(s.plus)
   sd.minus.statistic <- sd(s.minus) 
   s.plus.null <- apply(samroc@d0, 2, maxplus)
   s.minus <- max.minus(samroc@d)
   mean.plus.statistic <- mean(s.plus)
   mean.minus.statistic <- mean(s.minus)
   sd.plus.statistic <- sd(s.plus)
   sd.minus.statistic <- sd(s.minus) 
   maxmeanf <- function(x) { plusstat <- maxplus(x);minusstat <- maxminus(x);
                          plusmean <- colMeans(plustat);
                          minusmean <- colMeans(minusstat);
                          rowMax(cbind(plusmean, minusmean))}
} else{
   statistic <- switch(EXPR = type, absolute = abs(samroc@d), original = samroc@d)
   null.statistic <- switch(EXPR = type, absolute = abs(samroc@d0), original = samroc@d0)
   mean.statistic <- mean(statistic);sd.statistic <- sd(statistic) 
   mean.null.statistic <- mean(null.statistic);
   sd.null.statistic <- mean(apply(null.statistic, 2, sd))
   restand.stat <- function(x) (x - mean.statistic)/sd.statistic
   restand.null.stat <- function(x) (x - mean.null.statistic)/sd.null.statistic   
}

contrast <- samroc@contrast
if(is(data, "ExpressionSet")) data <- exprs(data)
path.n = sapply(pway, length)
pway = pway[(path.n > min(cutoff))&(path.n < max(cutoff))]
used.probesets <- unique(unlist(pway))
used.probesets <- used.probesets[used.probesets %in% probeset]
n.used <- sapply(pway, function(x) sum(used.probesets %in% x))
pway <- pway[n.used > min(cutoff)]
data <- as.matrix(data);rownames(data) <- probeset
used.data <- data[probeset %in% used.probesets,]
statistic = statistic[probeset %in% used.probesets]
statistic <- restand.stat(statistic)
null.statistic = null.statistic[probeset %in% used.probesets,]
null.statistic <- restand.null.stat(null.statistic)
probeset <- probeset[probeset %in% used.probesets]
b <- B - ncol(null.statistic);
xmat <- samroc@design
xhat <- solve(t(xmat)%*%xmat)%*% t(xmat)
ss <- samroc@se
smoothinteger <- smooth+0
npar <- qr(xmat)$rank      
scalek <- (contrast%*%solve(t(xmat)%*%xmat)%*%contrast)^-1

analysis <- function(x){
  obs.stat <- mean(statistic[probeset %in% x]);nulls <- null.statistic[probeset %in% x,]
   x.data <- used.data[probeset %in% x,]
#c-code  = bootstraploop
     result<-.C("newboot",
                as.double(x.data),
                        as.integer(nrow(x.data)),
                as.integer(ncol(x.data)),
                        as.double(xhat),
                        as.integer(nrow(xhat)),
                as.integer(ncol(xhat)),
                        as.double(xmat),
                        as.double(contrast),
                as.integer(smoothinteger),
                  as.integer(b),
                as.double(npar),
                as.double(scalek),
                as.double(ss),
                resultdiffs=as.double(rep(0,(nrow(x.data))*b)),
                        resultsses=as.double(rep(0,(nrow(x.data))*b)),
                        length(contrast = 0),
                PACKAGE = "SAGx" )                
    diffs <- matrix(result$resultdiffs,ncol=b)
    sses <- matrix(result$resultsses,ncol=b)                
    dstari <- diffs/(samroc@s0 + sses);if(type=="absolute") dstari <- abs(dstari)
    dstari <- restand.null.stat(dstari)
    all.nulls <- cbind(dstari, nulls)
    null.stats <- colMeans(all.nulls)
    Fn <- ecdf(-abs(null.stats))
    pj <- (B*Fn(-abs(obs.stat))+1)/(B+1)
    pj
}

out <- sapply(pway, analysis)
out
}

 
 

 

