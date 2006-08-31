# After J.Storey A Direct Approach to the False Discovery Rate (2001)
estimatep0 <- function(ps = pp, B = 500, range = seq(0,0.95, by = 0.05)){
# B <- 500;range <- seq(0,0.95, by = 0.05)
Bp <- matrix(nrow = B, ncol = length(ps))
esti1 <- numeric(length(range))
ests <- matrix(nrow = B, ncol = length(range))

for (i in 1:B) Bp[i,] <- ps[sample(1:length(ps), replace = TRUE)]
for (j in 1:length(range)) esti1[j] <- mean(ps > range[j])/(1-range[j])
est1 <- min(esti1) # The plug-in estimate 

for (i in 1:B) {for (j in 1:length(range)) 
                ests[i,j] <- mean(Bp[i,] > range[j])/(1-range[j])} 

SS <- apply(ests, 2, function(x) sum((x - est1)^2))

lambda <- range[which.min(SS)] # find which cut-off minimises the SS #

p0 <- mean(ps > lambda)/(1-lambda)
list(p0 = p0)}
