# Calculation of fp and fn from a vector of p-values #
# Bug fix 21NOV03#
fp.fn <- function(ps = pvals, B = 100){
        p0 <- estimatep0(ps = ps, B = B)$p0
        fp <- p0 * ps
        fn <- 1 - p0 * (1 - ps) - rank(ps)/length(ps)
        errors <- fp + fn
        N <- which.min(errors[sort.list(ps)])
        list(p0 = p0, fp = fp, fn = fn, errors = errors, N = N)
}
