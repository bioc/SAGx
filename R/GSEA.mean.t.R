GSEA.mean.t <- function(samroc = samroc.res, probeset = probeset, pway = kegg, type = c("original","absolute", "maxmean"),  two.side = FALSE, cutoff = c(10,Inf), restand = TRUE){
# Author : PBG 26OCT05, based on Tian et al. PNAS 2005.Modified 19MAY06. Use of design slot added 16OCT06, restandardisation 20OCT06
# Version 07FEB07
# data = golub; samroc = samroc.res; probeset = probeset; pway = kegg; type = "maxmean"; two.side = FALSE; cutoff = c(10,Inf); B = 1000;smooth = FALSE;restand = TRUE
if(type == "maxmean") stop("maxmean not implemented yet.") 
type <- match.arg(type) 
if(is(data, "ExpressionSet")) data <- exprs(data)
statistic <- samroc@d
null.statistic <- samroc@d0
# if maxmean restand should be done separately for s- and s+
if(type == "maxmean") {
   maxplus <- function(x) pmax(x, 0)
   maxminus <- function(x) -pmin(x, 0)
   s.plus <- maxplus(statistic)
   s.minus <- maxminus(statistic)
   mean.plus.statistic <- mean(s.plus)
   mean.minus.statistic <- mean(s.minus)
   sd.plus.statistic <- sd(s.plus)
   sd.minus.statistic <- sd(s.minus) 
   s.plus.null <- maxplus(null.statistic)
   s.minus.null <- maxminus(null.statistic)
   mean.plus.statistic <- mean(s.plus)
   mean.minus.statistic <- mean(s.minus)
   sd.plus.statistic <- sd(s.plus)
   sd.minus.statistic <- sd(s.minus) 
   mean.plus.null.statistic <- mean(s.plus.null)
   mean.minus.null.statistic <- mean(s.minus.null)
   sd.plus.null.statistic <- sd(s.plus.null)
   sd.minus.null.statistic <- sd(s.minus.null) 
   maxmeanf <- function(x) { plusstat <- maxplus(x);minusstat <- maxminus(x);
                          plusmean <- (colMeans(plustat)-mean.plus.statistic)/sd.plus.statistic;
                          minusmean <- (colMeans(minusstat)-mean.minus.statistic)/sd.minus.statistic;
                          rowMax(cbind(plusmean, minusmean))}
} else{
   statistic <- switch(EXPR = type, absolute = abs(statistic), original = statistic)
   null.statistic <- switch(EXPR = type, absolute = abs(null.statistic), original = null.statistic)
   mean.statistic <- mean(statistic);sd.statistic <- sd(statistic) 
   mean.null.statistic <- mean(null.statistic);
   sd.null.statistic <- mean(apply(null.statistic, 2, sd))
   restand.stat <- function(x) (x - mean.statistic)/sd.statistic
   restand.null.stat <- function(x) (x - mean.null.statistic)/sd.null.statistic   
}

path.n = sapply(pway, length)
pway = pway[(path.n > min(cutoff))&(path.n < max(cutoff))]
used.probesets <- unique(unlist(pway))
used.probesets <- used.probesets[used.probesets %in% probeset]
n.used <- sapply(pway, function(x) sum(used.probesets %in% x))
pway <- pway[n.used > min(cutoff)]
statistic = statistic[probeset %in% used.probesets]
if(restand) statistic <- restand.stat(statistic)
null.statistic = null.statistic[probeset %in% used.probesets,]
if(restand) null.statistic <- restand.null.stat(null.statistic)
probeset <- probeset[probeset %in% used.probesets]

analysis <- function(x){
# Note that statistic has been restandardised
  gs.statistics <- statistic[probeset %in% x] 
  N <- length(gs.statistics)
  mean.stat <- mean(gs.statistics)
  median.stat <- median(gs.statistics)
  nulls <- null.statistic[probeset %in% x,]
    null.stats <- colMeans(nulls);null.mean <- mean(null.stats);null.sd <- sd(null.stats)
  if(type == "original"){ranked.nulls  <-  apply(nulls, 2, rank);
  ranked.nulls <- ranked.nulls*(nulls > 0);rank.sums <- colSums(ranked.nulls)
  rank.null.mean  <- mean(ranked.nulls)*N;rank.null.sd <- sd(rank.sums)
  rank.stat <- sum(rank(gs.statistics)[gs.statistics > 0]);pw <- pnorm(rank.stat, mean = rank.null.mean, sd = rank.null.sd)
  pw <- 2*pmin(pw, 1 - pw)}  
#    pj <- (B*Fn(-abs(obs.stat))+1)/(B+1)
    pj <- pnorm(mean.stat, mean = null.mean, sd = null.sd)
    pj <- 2*pmin(pj, 1 - pj)
    if(type == "original") list(c(pj, mean.stat, pw, median.stat))
    else list(c(pj, mean.stat, median.stat))
}
out <- sapply(pway, analysis)
names(out) <- names(pway)
if(type == "original"){x1 <- sapply(out, function(x) x[1]);x2 <- sapply(out, function(x) x[2]);
x3 <- sapply(out, function(x) x[3]);x4 <- sapply(out, function(x) x[4])
out <- cbind(x1, x2, x3, x4);colnames(out) <- c("normal p-value", "mean statistic", "Wilcoxon p-value", "median statistic")}
else {
x1 <- sapply(out, function(x) x[1]);x2 <- sapply(out, function(x) x[2]);
x3 <- sapply(out, function(x) x[3]);out <- cbind(x1, x2, x3)
colnames(out) <- c("normal p-value", "mean statistic", "median statistic")
}

return(out)
}

 

 
 

 

