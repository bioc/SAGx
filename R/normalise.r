# Normalisation according to Tibshirani et al #
# Author : PBG 2002; corrected 20JAN05
normalise <- function (x, linear = TRUE) {
#   x <- dats
    x <- as.matrix(x)
#    y <- apply(x, 1, mean)
    y <- t(c(rep(1/ncol(x), ncol(x)))%*%t(x))
    x <- sign(x) * (abs(x))^(1/3)
    y <- sign(y) * (abs(y))^(1/3)
    for (i in 1:ncol(x)) {
        x[, i] <- predict(lm(y ~ x[, i]), as.data.frame(x[, i]))
    }
    if (linear) 
        return(x^3)
    else return(x)
}
