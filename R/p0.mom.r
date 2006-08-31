# Author Per Broberg 07JAN04, Copyright Per Broberg
# Changed 03MAY04, 11AUG04, 15NOV04 #################################

p0.mom <- function(ps = pvalues){
# ps <- pvals
   R.s <- function(ps, s){
             sum(exp(s*ps))/length(ps)
                         }
   ss <- seq(0.01, 0.1, by = 0.001)
   R.s.s <- vector(length = length(ss) ) 
   R.s.s <- apply(as.data.frame(ss), 1, function(x) R.s(ps, x))
   g <- function(x) (exp(x) - 1)/x
r.n <- vector(length = length(ss))
g.n <- g(ss)
# Since r.n < R.s.s < g.n it is reasonable to choose r.n[1] between 1 and R.s.s[1]
# r.n[1] <-  mean(c(1,R.s.s[1]))

r.n.est <- function(x){
   temp.vect <- vector(length =length(ss))
   temp.vect[1] <- x
   for(i in 2:length(ss)) temp.vect[i] <- (R.s.s[i]*g.n[i-1] - R.s.s[i-1]*g.n[i] + temp.vect[i-1]*(g.n[i] - R.s.s[i]))/(g.n[i-1] - R.s.s[i-1])
   temp.vect
}   

# p0.0 <- min(mean((R.s.s-r.n)/(g.n-r.n)),1)
#s1 <- ss;s2 <- ss**2;s3 <- ss**3;s4 <- ss**4;s5 <- ss**5
#lm.ps <- lm(log(r.n) ~ s1 + s2 + s3 + s4 + s5 -1) 

#a.s <- coefficients(lm.ps)

#pred <- function(x) exp(x[1]*s1+x[2]*s2+x[3]*s3+x[4]*s4+x[5]*s5)

#target <- function(x) var((R.s.s - pred(x))/(g.n - pred(x)))

target <- function(x) sd((R.s.s - r.n.est(x))/(g.n - r.n.est(x)))/mean((R.s.s - r.n.est(x))/(g.n - r.n.est(x)))

optimal.value <- optimize(target, c(1,min(R.s.s)))$minimum

r.n <- r.n.est(optimal.value)

# p0.0 <- mean((R.s.s - pred(optimal.value))/(g.n - pred(optimal.value)))

p0.0 <- min(mean((R.s.s-r.n)/(g.n-r.n)),1)

# First smooth the density for ps 
   counts <- hist(ps, plot = FALSE, breaks = 'FD')

   K <- length(counts$mids)
   M <- matrix(ncol = K, nrow = K)
   lambda <- 0.1
   for(i in 1:K) for(j in 1:K) M[i,j] <- dnorm((counts$mids[i]-counts$mids[j])/lambda)/lambda
   row.sum <- matrix(rep(1, ncol(M))%*%M)
   M <- sweep(M, 1, row.sum, FUN = "/")
   ##########################################
      mu <- M%*%counts$counts
   x1 <- counts$mids;x2  <- x1^2;x3 <- x1^3;x4 <- x1^4
   glm.test <- glm(counts$counts ~ x1 + x2 + x3 + x4, family  = poisson(), offset = log(mu))
   spl <- smooth.spline(counts$mids, glm.test$fitted, df = 4)
   x1 <- ps;x2  <- x1^2;x3 <- x1^3;x4 <- x1^4
   fs <- length(counts$counts)*predict(spl, ps)$y/sum(counts$counts)
   ###########################################

p0.1 <- min(fs)

# ifelse(p0.1 < 0.8, r.n[1] <-  sum(exp(ss[1]*ps)*(1 - p0.1/fs))/sum((1 - p0.1/fs))), r.n[1] <-  mean(c(1,R.s.s[1])))

p0.4 <- min(min(pava(rev(counts$density))),1)


if(p0.0 > 0.95) p0.1 <- min(1, quantile(fs, probs = c(0.25)))

if(p0.1 < 1){
r.n[1] <- (R.s.s[1] - p0.1 * g.n[1])/(1 - p0.1)
for(i in 2:length(ss)) r.n[i] <- (R.s.s[i]*g.n[i-1] - R.s.s[i-1]*g.n[i] + r.n[i-1]*(g.n[i] - R.s.s[i]))/(g.n[i-1] - R.s.s[i-1])
}

p.test <- seq(0, 1, by = 0.005)
error.obs <- apply(as.data.frame(p.test), 1, function(x) mean((R.s.s - x*g.n - (1-x)*r.n)^2)) 
p0.2 <- p.test[which.min(error.obs)]

# w <- 0.5
# if(p0.1 < 1) w <- p0.1*(R.s.s[1]-g.n[1])/(1-R.s.s[1])/(1-p0.1)

# r.n[1] <- w + (1 - w) * R.s.s[1]
# for(i in 2:length(ss)) r.n[i] <- (R.s.s[i]*g.n[i-1] - R.s.s[i-1]*g.n[i] + r.n[i-1]*(g.n[i] - R.s.s[i]))/(g.n[i-1] - R.s.s[i-1])
# p0.2 <- min(mean((R.s.s-r.n)/(g.n-r.n)),1)
ifelse(p0.0 < 0.8, x <- p0.1, x <- p0.0)
r.n[1] <- (R.s.s[1] - x * g.n[1])/(1 - x)
for(i in 2:length(ss)) r.n[i] <- (R.s.s[i]*g.n[i-1] - R.s.s[i-1]*g.n[i] + r.n[i-1]*(g.n[i] - R.s.s[i]))/(g.n[i-1] - R.s.s[i-1])
p0.3 <- min(mean((R.s.s-r.n)/(g.n-r.n)),1);if(is.na(p0.3)) p0.3 <- 1
list(mgf = p0.0, PRE = p0.1, p0.2 = p0.2, p0.3 = p0.3, p0.4 = p0.4)
}
