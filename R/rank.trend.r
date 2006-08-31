# Trend analysis
rank.trend <- function(data = x, pattern = c(1:ncol(data)), har = FALSE){
# data <- medians;pattern <- c(1:ncol(data));har <- FALSE
    rs <- t(apply(data,1,function(x) rank(x)))
    diff <- scale(rs,center = pattern, scale = FALSE)
    if(har) hardy <- apply(rs, 1, function(x) pattern %*% x) else hardy <- NULL
    sds <- apply(rs,1, function(x) {ds <- table(x);sds <- sum(ds^3-ds);return(sds)})
    score <- apply(diff,1, function(x) sum(x^2)) 
    N <- ncol(data)
    ED <- (N^3-N)/6-sds/12
    SD <- sqrt(((N-1)*N^2*(N+1)^2/36)*(1-sds/(N^3-N)))
    pvals <- pnorm((score-ED)/SD)
    ut <- cbind(score,hardy, pvals)
    colnames(ut) <- c("Trend score","Hardy score","p-value for no trend")
    return(ut)
}
