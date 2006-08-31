# Author PBG 11NOV04; updated 09MAY05; updated 14JUN05
pava.fdr <- function(ps = pvalues, p0 = NULL){
if(is.null(p0)) p0 <- p0.mom(ps = ps)$PRE
order.ps <- order(ps)
pvals.ord <- ps[order.ps]
diff.p <- diff(pvals.ord)
local.fdr <- length(pvals.ord)*p0*isoreg(c(pvals.ord[1], diff.p))$yf
local.fdr <- pmin(local.fdr, 1)
i <- 1:length(local.fdr)
q.value <- 1/i * cumsum(local.fdr)
back.order <- order(order.ps) 
list(pava.fdr = q.value[back.order], p0 = p0, pava.local.fdr = local.fdr[back.order])
}
