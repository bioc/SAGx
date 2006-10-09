# Author : Per Broberg 20APR04
# Based on Broberg P, Estimation of Relative Standard Deviation,(1999) in Drug Development and Industrial Pharmacy, Vol 25 no 1 37-43

rsd.test <- function(data1 = x, data2= y, B = NULL){
# data1 <- dats.ns;data2 <- dats.copd;B <- 500
nrows <- nrow(data1) 
rsd.func <- function(x, y) {
#  x <- data1;y <- data2
  N1 <- ncol(x);N2 <- ncol(y)
 # rowMeans introduced
  m1 <- rowMeans(x)
  m2 <- rowMeans(y)
  resid1 <- sweep(x, 1, m1)
sd1 <- t(sqrt(crossprod(c(rep(1, N1)), t(resid1 * resid1))/(N1 - 1)))
  resid2 <- sweep(y, 1, m2)
sd2 <- t(sqrt(crossprod(c(rep(1, N2)), t(resid2 * resid2))/(N2 - 1)))
RSD1 <- 100*sd1/m1;RSD2 <- 100*sd2/m2
t.stat <- (RSD1 - RSD2)/sqrt(RSD1**4/N1+RSD1**2/(2*N1) + RSD2**4/N2+RSD2**2/(2*N2))
list(cv1 = RSD1, cv2 = RSD2, t.stat = t.stat)
}

   rsd.obs <- rsd.func(data1, data2)
   if(!is.null(B)){
      t.stat.null <- matrix(nrow = nrow(data1), ncol = B) 
      dats <- data.frame(data1, data2)
      for(i in 1:B){
         samp <- sample(1:ncol(dats), ncol(data1))
         t.stat.null[,i] <- rsd.func(dats[, samp], dats[,-samp])$t.stat
         }
         abs.stat.null <- abs(t.stat.null)
         Fn <- ecdf(-abs(as.vector(abs.stat.null)))
         p.vals <- (nrows*B*Fn(-abs(rsd.obs$t.stat))+1)/(nrows*B+1)
         } else p.vals <- NULL
         list(cv1  = rsd.obs$cv1, cv2 = rsd.obs$cv2, t.stat = rsd.obs$t.stat, p.vals = p.vals)
}
