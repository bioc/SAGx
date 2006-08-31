# Calculation of F statistic #
# rewritten 14OCT02, modified 13NOV02, xmat0 corrected 29OCT03, re-written with linear algebra 08MAR05  #
Fstat <- function (indata = M, formula1 = ~as.factor(g), formula0 = "mean",  design1 = NULL, design0 = NULL, B = NULL) 
{
# indata <- dats;formula1 <- ~lavage + smoke;design0 <- model.matrix(~lavage);formula0 <- "mean";B <- 50;design1 <- NULL 
    nobs <- ncol(indata)
    ifelse(is.null(design1), xmat1 <- model.matrix(formula1),   xmat1 <- design1)
    if(is.null(design0)&is.null(formula0)) stop("Please provide null model")
    ifelse((formula0 == "mean"),  xmat0 <- model.matrix(1:nobs~1),
                                  ifelse(is.null(design0),xmat0 <- model.matrix(formula0), xmat0 <- design0))                                                
    npar1 <- qr(xmat1)$rank
    ifelse(formula0 == "mean", npar0 <- 1, npar0 <- qr(xmat0)$rank)

 # New 08MAR05  
   xhat1 <- solve(t(xmat1)%*%xmat1)%*% t(xmat1)
   beta <- crossprod(t(xhat1), t(indata))
   yhat1 <- crossprod(beta, t(xmat1))
   resid1 <- indata - yhat1    
   totSS.out <- crossprod(c(rep(1, nobs)),t(resid1 * resid1)) 
   xhat0 <- solve(t(xmat0)%*%xmat0)%*% t(xmat0)
   beta0<- crossprod(t(xhat0), t(indata))
   yhat0 <- crossprod(beta0, t(xmat0))
   resid2 <- yhat1 - yhat0
   SSmodel.out <- crossprod(c(rep(1, nobs)),t(resid2 * resid2))
   ##########################

    Fstatv <- (nobs - npar1) * SSmodel.out/((npar1 - npar0) * totSS.out)
    p.value <- numeric(nrow(indata))
    if(is.null(B)) p.value <- NULL
    null.resids <- indata - yhat0
    if (!is.null(B)) {
    Fstat0.res <- matrix(nrow = nrow(indata), ncol = B)
      for(i in 1:B){
        simdata <- yhat0 + null.resids[, sample(1:ncol(indata))]
        beta.star <- crossprod(t(xhat1), t(simdata))
        yhat1.star <- crossprod(beta.star, t(xmat1))
        resid1.star <- simdata - yhat1.star 
        totSS.star <- crossprod(c(rep(1, nobs)),t(resid1.star * resid1.star))
        beta0.star <- crossprod(t(xhat0), t(simdata))
        yhat0.star <- crossprod(beta0.star, t(xmat0))
        resid2.star <- yhat1.star - yhat0.star
        SSmodel <- crossprod(c(rep(1, nobs)),t(resid2.star * resid2.star))   
        Fstat0 <- (nobs - npar1) * SSmodel/((npar1 - npar0) * totSS.star)
        Fstat0.res[,i] <- Fstat0       
    }}
# as.vector introduced 12OCT05
    if (!is.null(B)) {Fn <- ecdf(as.vector(Fstat0.res)); p.value <- 1 - Fn(Fstatv)}
    list(Fstat = Fstatv, fnum = nobs - npar1, fdenom = npar1 -npar0, design1 = xmat1, design0 = xmat0, SS1 = totSS.out, SS0 = SSmodel.out, pvalue = p.value)}
