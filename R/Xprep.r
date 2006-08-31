# Preparation of data for RBayEst starting from linear model, e.g. effect~group #
# Author: Per Broberg 14AUG01 #
# History : modified 23SEP01, 20JAN02, 22FEB04 #

Xprep<-function(indata=M,formula=~as.factor(g),contrast=c(0,1),design=NULL){
# indata <- utmat;formula <-  ~as.factor(g);contrast <- c(0,1);design <- NULL
   if(is.null(design)) xmat <- model.matrix(formula) else xmat <- design
   npar <- qr(xmat)$rank;nobs <- ncol(indata)
   xhat <- solve(t(xmat)%*%xmat)%*% t(xmat)
   beta <- crossprod(t(xhat), t(indata))
   est.contrast <- crossprod(contrast, beta)
   yhat <- crossprod(beta, t(xmat))
   resid <- indata - yhat
   var <- crossprod(c(rep(1, ncol(resid))),t(resid * resid))/(nobs - npar) 
   scalek <- (contrast%*%solve(t(xmat)%*%xmat)%*%contrast)^-1
  list(Mbar = est.contrast,Vest = as.vector(var),k = scalek,f = nobs-npar,design = xmat, xhat = xhat)
}
