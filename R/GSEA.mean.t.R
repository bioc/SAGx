GSEA.mean.t <- function(data = X, samroc = samroc.res, probeset = probeset, pway = kegg, absolute = TRUE, two.side = FALSE, cutoff = c(5,Inf), B = 1000, smooth = FALSE){
# Author : PBG 26OCT05, based on Tian et al. PNAS 2005.Modified 19MAY06
# data = X.Y; samroc = samroc.res;statistic = samroc.res$d;null.statistic = samroc.res$d0;pway = kegg;two.side = TRUE;probeset = rownames(X.Y);B <- 10000
# data = matrix(rnorm(nrow(X.Y)*ncol(X.Y)), ncol = ncol(X.Y))
# samroc.res <- samrocNboot(data = data, formula = ~as.factor(cl), B = 20)
# samroc = samroc.res;samroc$formula <- ~as.factor(cl); samroc$contrast <- c(0,1);smooth = FALSE
# data = dats.COPD.HS.path; samroc = samroc.COPD.HS.path; probeset = probeset.sel; pway = all.gsets; formula = ~ late.C.HS+g.C.HS; contrast = c(0,0,0,1); absolute = TRUE; two.side = FALSE; cutoff = c(5,200); B = 1000
#formula <- samroc$formula 
#contrast <- samroc$contrast
if(smooth) stop("Smoothing not implemented yet.")
if(!absolute) stop("Only absolute value of test statistic of test statistic currently used.")
if(two.side) stop("Only one-sided test currently used base on absolute value of test statistic.")
formula <- samroc@formula
contrast <- samroc@contrast
path.n = sapply(pway, length)
pway = pway[(path.n > min(cutoff))&(path.n < max(cutoff))]
used.probesets <- unique(unlist(pway));used.probesets <- used.probesets[used.probesets %in% probeset]
n.used <- sapply(pway, function(x) sum(used.probesets %in% x))
pway <- pway[n.used > min(cutoff)]
data <- as.matrix(data);rownames(data) <- probeset
used.data <- data[probeset %in% used.probesets,]
statistic = samroc@d[probeset %in% used.probesets]
null.statistic = samroc@d0[probeset %in% used.probesets,]
# probeset <- probeset[probeset %in% used.probesets]
b <- B - ncol(null.statistic);
sam <- model.matrix(formula)
xmat <- model.matrix(formula)
utest.design <- model.matrix(formula) 
xhat <- solve(t(xmat)%*%xmat)%*% t(xmat)
utest <- Xprep(indata = data, formula = formula, contrast = contrast)
ss <- sqrt(utest$Vest/utest$k)
smoothinteger <- smooth+0
npar <- qr(xmat)$rank      
scalek <- (contrast%*%solve(t(xmat)%*%xmat)%*%contrast)^-1


analysis <- function(x){
   obs.stat <- mean(abs(statistic[used.probesets %in% x]));nulls <- null.statistic[used.probesets %in% x,]
   x.data <- used.data[used.probesets %in% x,]
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
                
    diffs<-matrix(result$resultdiffs,ncol=b)
    sses<-matrix(result$resultsses,ncol=b)                
    dstari <- diffs/(samroc@s0 + sses)
    all.nulls <- cbind(dstari, nulls)
    if(absolute) null.stats <- colMeans(abs(all.nulls)) else null.stats <- colMeans(all.nulls)
    Fn <- ecdf(-abs(null.stats))
    pj <- (B*Fn(-abs(obs.stat))+1)/(B+1)
    pj
}

out <- sapply(pway, analysis)
out
}
