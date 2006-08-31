# Author PBG/28APR06
Xprep.resid <- function (data = M, formula = ~as.factor(g), design = NULL) 
{
    if (is.null(design)) 
        xmat <- model.matrix(formula)
    else xmat <- design
    npar <- qr(xmat)$rank
    nobs <- ncol(data)
    xhat <- solve(t(xmat) %*% xmat) %*% t(xmat)
    beta <- crossprod(t(xhat), t(data))
    yhat <- crossprod(beta, t(xmat))
    resid <- data - yhat
    resid
}           
