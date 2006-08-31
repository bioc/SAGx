# The Lehmann Nonparametrics: Statistical Methods Based on Ranks p. 233 #
# Assumes that groups are coded in increasing numerical order #

JT.test <- function(data = utmat,class = g, labels = NULL){

n.obs <- ncol(data)

class.tab <- unique(class);class.tab <- class.tab[order(class.tab)]

wilcox.stat <- function(x,y){
        r <- rank(c(x, y))
        n.x <- length(x)
        n.y <- length(y)
        STATISTIC <- sum(r[seq(along = x)]) - n.x * (n.x + 1) / 2  
        return(STATISTIC)
}


sum.stat <- function(x,class){
   sums <- 0
   upper <- length(class.tab)-1
   for(i in 1:upper){
      for(j in i+1:length(class.tab)){
         sums <- sums+wilcox.stat(x[class==class.tab[i]],x[class==class.tab[j]])
      }
   }
sums                 
}

EH <- (n.obs^2 - sum(table(class)^2))/4

STDH <-sqrt((n.obs^2*(2*n.obs+3)-sum(table(class)^2*(2*table(class)+3)))/72 )

utres <- t(apply(data, 1, function(x) c(2*min(pnorm((sum.stat(x,class)-EH)/STDH),1-pnorm((sum.stat(x,class)-EH)/STDH)),tapply(x,class,median),cor(rank(x),rank(class) ) )))

colnames(utres)[c(1,ncol(utres))] <- c("trend p-value", "rank correlation");

ifelse(is.null(labels),colnames(utres)[-c(1,ncol(utres))] <- class.tab, colnames(utres)[-c(1,ncol(utres))] <- labels)
is.null(labels)
#list(pvalues=utres[,1],medians=utres[,-c(1,ncol(utres))],trend = utres[,ncol(utres)])
return(utres)
}
