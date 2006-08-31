# samrocNboot, author Per Broberg, version 31AUG02, modified 13NOV02, 27MAY04, 31MAY06  #
samrocNboot <- function (data = M, formula = ~as.factor(g), contrast = c(0,1), N = c(50, 100, 200, 300), B = 100, perc = 0.6, smooth = FALSE, w = 1, measure = "euclid", probeset = NULL) 
{
    data <- as.matrix(data)
    nrows <- nrow(data);ncols = ncol(data)
#  data <- utmat;formula <- ~as.factor(g);contrast <- c(0, 1);N <- c(50,100,200,300)
#  B <- 100;perc <- 0.6;smooth <- F;measure <- c("euclid");w <- 1
    utest <- Xprep(indata = data, formula = formula, contrast = contrast)
    ss <- sqrt(utest$Vest/utest$k)
    minss <- min(ss)
    library(stats)
    if (smooth) {
        totmeans <- c(rep(1/ncol(data), ncol(data)))%*%t(data)
        tv <- smooth.spline(totmeans, ss, spar = 1.5)
        ss <- predict(tv, totmeans)$y
        ss[ss < minss] <- minss
    }
    ssq <- quantile(ss, probs = seq(0, perc, by = 0.05))
    ssq <- c(0, ssq)
    ss <- sqrt(utest$Vest/utest$k)
    diffs <- matrix(nrow = nrows, ncol = B)
    diffs1<-t(diffs)
    sses <- matrix(nrow = nrows, ncol = B)
    varses <- matrix(nrow = nrows, ncol = B)
    pj <- numeric(length = nrows)
    pts <- matrix(nrow = length(ssq), ncol = length(N))
    p.is <- matrix(nrow = length(ssq), ncol = length(N))
    target <- numeric(length(ssq))
    xmat <- model.matrix(formula)
    npar <- qr(xmat)$rank       
    xhat <- solve(t(xmat)%*%xmat)%*% t(xmat)
    scalek <- (contrast%*%solve(t(xmat)%*%xmat)%*%contrast)^-1
    ifelse(smooth, smoothinteger<-1,smoothinteger<-0)
        
#c-code  = bootstraploop
     result<-.C("samrocNboot",
                as.double(data),
                        as.integer(nrow(data)),
                as.integer(ncol(data)),
                        as.double(xhat),
                        as.integer(nrow(xhat)),
                as.integer(ncol(xhat)),
                        as.double(xmat),
                        as.double(contrast),
                as.integer(smoothinteger),
                  as.integer(B),
                as.double(npar),
                as.double(scalek),
                as.double(ss),
                resultdiffs=as.double(rep(0,(nrow(data))*B)),
                        resultsses=as.double(rep(0,(nrow(data))*B)),
                PACKAGE = "SAGx" )
    
       diffs<-matrix(result$resultdiffs,ncol=B)
       sses<-matrix(result$resultsses,ncol=B)
  
    for (i in 1:length(ssq)) {
        dstari <- diffs/(ssq[i] + sses)
        di <- t(utest$Mbar/(sqrt(utest$Vest/utest$k) + ssq[i]))
        di.min <- quantile(abs(di), probs = (nrow(data)-N)/nrow(data), na.rm = TRUE)
        library(stats)
        Fn <- ecdf(abs(as.vector(dstari)))
        alpha <- 1 - Fn(di.min)
        pj <- 1 - Fn(abs(di))
#        test <- sweep(abs(dstari),1, abs(di), FUN = "-")
#        pj <- t(crossprod(c(rep(1/ncol(test), ncol(test))), t(test > 0)))
#        pj <- apply(as.matrix(di), 1, function(x) mean(abs(dstari) >  abs(x)))
        pt <- numeric(length(N))
#        p.i <- apply(as.matrix(alpha), 1, function(x) max(pj[pj <  x]))
#        p.is[i, ] <- p.i
        p.is[i, ] <- alpha
        for (j in 1:length(alpha)) pt[j] <- mean((pj < alpha[j]))
           pts[i, ] <- pt
    }
    pfn <- matrix(nrow = length(ssq), ncol = length(alpha))
    tone <- t(utest$Mbar/ss)
    tnull <- diffs/sses
    Fn <- ecdf(abs(as.vector(tnull)))
    ps <- 1 - Fn(abs(tone))
#    test <- sweep(abs(tnull),1, abs(tone), FUN = "-")
#    ps <- t(crossprod(c(rep(1/ncol(test), ncol(test))), t(test > 0)))
#    ps <- apply(as.matrix(tone), 1, function(x) mean(abs(tnull) > abs(x)))
    p0 <- p0.mom(ps = ps)$PRE
    if(p0 >= 0.999) p0 <- 0.999
    for (j in 1:length(ssq)) pfn[j, ] <- 1 - p0 * (1 - alpha) -  pts[j, ]
    pfn[pfn < 0] <- 0
    pfp <- p0 * p.is
    pfp[pfp == -Inf] <- 0
    if(measure == "euclid") criterion <- function(x, y , w) sqrt(x^2 + w*y^2)
    if(measure=="sum")   criterion <- function(x, y, w) x + w*y  
    ifelse(pfn == 0, distances <- p0 * alpha/pts, distances <- criterion(pfn, pfp, w) )
    target <- min(distances)
    indices <- which(distances == target, arr.ind = TRUE)
    indices <- indices[1, ]
    snew <- ssq[indices[1]]
    N.opt <- N[indices[2]]
    di <- t(utest$Mbar/(sqrt(utest$Vest/utest$k) + snew))
    dstari <- diffs/(snew + sses)
    Vest0 <- as.numeric(utest$k) * sses^2
#    pj <- apply(as.matrix(t(di)), 1, function(x) mean(abs(dstari) > abs(x)))
    Fn <- ecdf(abs(as.vector(dstari)))
#    test <- sweep(abs(dstari),1, abs(di), FUN = "-")
    pj <- 1 - Fn(abs(di))
#    pj <- t(crossprod(c(rep(1/ncol(test), ncol(test))), t(test > 0)))
        fp <- p0 * pj
        fn <- 1 - p0 * (1 - pj) - rank(pj)/length(pj)
        errors <- fp + fn
    if(is.null(probeset)) probeset <- paste("probeset",1:nrows)
res <-   new("samroc.result", 
             d = as.vector(di), 
             diff = as.vector(utest$Mbar),
             se = as.vector(sqrt(utest$Vest/utest$k)),
             d0 = dstari,
             p0 = p0,
             s0 = snew,
             pvalues = pj,
             N.list = as.integer(N.opt),
             errors = errors,
             formula = formula,
             contrast = contrast,
             annotation = as.character(date()),
             N.sample = ncols,
             B = as.integer(B),
             call = as.character(match.call()),
             id = as.character(probeset),
             error.df = as.integer(utest$f)
)

res 
}
