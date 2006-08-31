list.intersection.p <- function(N = 14000, N1 = 100, N2 = 200, common = 30){
# N <- 14000; N1 <- 100; N2 <- 200; common <- 30
lambda <- N1*N2/N
p.int <- function(x) exp(lchoose(N,x) + ((lchoose(N-x,N1-x) - lchoose(N,N1)) + (lchoose(N-N1,N2-x)  - lchoose(N,N2))))
poisson.approx <- ppois(common-1,lambda)
if((poisson.approx > 0.99999)&(lambda < 15)) out <- 1-poisson.approx else 
   {sum.p <- 0;
      for(i in 0:common-1) sum.p <- sum(sum.p,p.int(i));out <- 1 - sum.p}
if(lambda > 15) out <- 1 - pnorm((common - lambda -0.5)/sqrt(lambda))
out

} 
