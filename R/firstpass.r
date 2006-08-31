### Per Broberg 25JAN01, modified 16JAN02, ######
### 10NOV03: data.frame(cbind( deleted and replaceds by data.frame( #############
### 18JAN05: A speed-up for Wilcoxon based on Kronecker products was put in place #
firstpass<-function( data = D, probes = probes, g, log = FALSE, present = NULL, labels = NULL, output.data = FALSE ){
   utdat<-data
   library(stats)
   kru<-function(x)kruskal.test(x,g)$p.val
   grps <- unique(g)
   grp.tab <- table(g)
   ifelse((length(grps) > 2), pval <- apply(data,1,kru),{
   result <- apply(data, 1, function(x) kronecker(x[g==grps[1]],x[g==grps[2]], FUN = "-") );
   result <- t(result);
   result <- rowSums(!(result +0 < 0));
   test <- pwilcox(result, grp.tab[2], grp.tab[1] );
   pval <- apply(as.data.frame(test), 1, function(x) 2*min(x, 1-x) )}
   )
   if(log==TRUE)data<-log2(data*(data>0)+1)
   rowm <- function(x)tapply(x,g,mean)
   rowsd <- function(x)tapply(x,g,sd)
   gmeans <- t(apply(data,1,rowm))
   gsd <- t(apply(data,1,rowsd))
   ns <- table(g)
   tet<-function(z){z*qt(.975,df=ns)/sqrt(ns)}
   clim<-apply(gsd,1,tet)
   lcl<-gmeans-t(clim);ucl<-gmeans+t(clim)
   ng<-dim(gsd)[2]
   colnames(gmeans)<-paste("Mean",1:ng);if(log==TRUE)gmeans<-2**gmeans
   colnames(lcl)<-paste("LCL",1:ng);if(log==TRUE)lcl<-2**lcl
   colnames(ucl)<-paste("UCL",1:ng);if(log==TRUE)ucl<-2**ucl
   if(log==TRUE)cat("Geometric means have been calculated\n")
   ifelse(length(table(g))==2, {diff<-gmeans %*% c(1,-1);
                                utmat <- data.frame(probes, gmeans, lcl, ucl, pval, diff);
                                colnames(utmat)[ncol(utmat)] <- c("Difference")},
                                utmat <- data.frame(probes, gmeans, lcl, ucl, pval) )
   if (!is.null(present)){pres<-apply(present,1,function(x) mean(x=3))
                       utmat <- data.frame(utmat,pres)}
   if(output.data) utmat <- data.frame(utmat, utdat)
   colnames(utmat)[1]<-c("Probesets")
   tvect <- c(1:ncol(utmat))
   ifelse(length(table(g))==2,colnames(utmat)[tvect[colnames(utmat) == "pval"]] <- c("Wilcoxon p-value"),
                              colnames(utmat)[tvect[colnames(utmat) == "pval"]] <- c("Kruskal-Wallis p-value"))
   if(!is.null(labels)){nnn <- ncol(gmeans)+1; colnames(utmat)[2:nnn] <- labels}  
   return(utmat)
}
