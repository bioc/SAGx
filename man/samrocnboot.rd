\name{samrocNboot}
\alias{samrocNboot}

\title{Calculate ROC curve based SAM statistic}

\description{A c-code version of samrocN. Calculation of the regularised t-statistic which minimises 
the false positive and false negative rates.}

\usage{samrocNboot(data=M,formula=~as.factor(g), contrast=c(0,1), N = c(50, 100, 200, 300), B=100, perc = 0.6, 
smooth=FALSE, w = 1, measure = "euclid", probeset = NULL)}

\arguments{
\item{data}{The data matrix}
\item{formula}{a linear model formula}
\item{contrast}{the contrast to be estimnated }
\item{N}{the size of top lists under consideration}
\item{B}{the number of bootstrap iterations}
\item{perc}{the largest eligible percentile of SE to be used as fudge factor}
\item{smooth}{if TRUE, the std will be estimated as a smooth function of expression level}
\item{w}{the relative weight of false positives}
\item{measure}{the goodness criterion}
\item{probeset}{probeset ids;if NULL then "probeset 1", "probeset 2", ... are used.}
}

\author{Per Broberg and Freja Vamborg}

\value{An object of class samroc.result.}

\details{The test statistic is based on the one in Tusher et al (2001):

\deqn{\frac{d = diff}{s_0+s}}{d = diff / (s_0 + s)}

where \eqn{diff} is a the estimate of a constrast, \eqn{s_0} is the regularizing constant 
and \eqn{s} the standard error.  At the heart of the method lies an estimate of the false negative and false positive rates. The
test is calibrated so that these are minimised. For calculation of \eqn{p}-values a bootstrap procedure is invoked. Further details are given in Broberg (2003).

The p-values are calculated through permuting the rows of the design matrix. NB This is not adequate for all linear models.

samrocNboot uses C-code to speed up the bootstrap loop. } 

\examples{
library(multtest)
#Loading required package: genefilter 
#Loading required package: survival 
#Loading required package: splines 
#Loading required package: reposTools 
data(golub)
 # This makes the expression data from Golub et al available
 # in the matrix golub, and the sample labels in the vector golub.cl
set.seed(849867)
samroc.res <- samrocNboot(data = golub, formula = ~as.factor(golub.cl))
# The proportion of unchanged genes is estimated at
samroc.res@p0
# The fudge factor equals
 samroc.res@s0
# A histogram of p-values
 hist(samroc.res@pvalues)
 # many genes appear changed
}

\keyword{htest}

\references{
Tusher, V.G., Tibshirani, R., and Chu, G. (2001) Significance analysis of microarrays applied to
the ionizing radiation response. \emph{PNAS} Vol. 98, no.9, pp. 5116-5121 \cr

Broberg, P. (2002) Ranking genes with respect to differential expression , \url{http://genomebiology.com/2002/3/9/preprint/0007} \cr

Broberg. P: Statistical methods for ranking differentially expressed genes. Genome Biology 2003, 4:R41
\url{ http://genomebiology.com/2003/4/6/R41}


}
