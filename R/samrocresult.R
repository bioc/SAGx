# definition of class samroc.result
setClass("samroc.result", 
  representation(
    d = "numeric", 
    diff = "numeric",
    se = "numeric",
    d0 = "matrix",
    p0 = "numeric",
    s0 = "numeric",
    pvalues = "numeric",
    N.list = "integer",
    errors = "numeric",
    formula = "formula",
    contrast = "vector",
    annotation = "character",
    N.sample = "integer",
    B = "integer",
    call = "character",
    id = "character",
    error.df = "integer",
    design = "matrix"
  )
)

setMethod("show", "samroc.result", function(object)
{
  nsamples <- object@N.sample
  ngenes <- length(object@d)
  cat("Samroc result:\n")
    cat("Data:", nsamples, "samples with", ngenes, "genes.\n")
  cat("Model:", as.character(object@formula), "\n")
  cat("Using", object@B, "permutations\n")
  cat("Fudge factor:",object@s0,". Estimated proportion unchanged genes:",format.pval(object@p0, digits = 2),".\n")
  cat("Annotation:",as.character(object@annotation),"\n")
  cat("Call:",object@call,"\n")
  cat("\n")
})

setMethod("plot", "samroc.result", function(x,y, ...)
{
plot.samroc.result <- function(x,...){
 dens0 <- density(x@d0)
 densobs <- density(x@d)
 x.coord <- quantile(dens0$x, probs = 0.55)
 y.coord <- quantile(c(dens0$y,densobs$y), probs = 1)
 plot(dens0, main = "", col = "green", xlab = "Statistic")
 lines(densobs, col = "red")
 legend(x.coord,y.coord, legend = c("Observed statistic", "Null distribution"), col = c("red","green"),
            lty = 1, cex = .8, y.int = 1)
 }
plot.samroc.result(x,...)
})

