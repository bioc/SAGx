\name{samrocN}
\alias{samrocN}

\title{Calculate ROC curve based SAM statistic}

\description{Calculation of the regularised t-statistic which minimises 
the false positive and false negative rates.}

\usage{samrocN(data=M,formula=~as.factor(g), contrast=c(0,1), N = c(50, 100, 200, 300),B=100, perc = 0.6, 
 smooth = FALSE, w = 1, measure = "euclid", probeset = NULL)}

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

\author{Per Broberg}

\value{An object of class samroc.result.}
}

\details{The test statistic is based on the one in Tusher et al (2001):

\deqn{\frac{d = diff}{s_0+s}}{d = diff / (s_0 + s)}

where \eqn{diff} is a the estimate of a constrast, \eqn{s_0} is the regularizing constant 
and \eqn{s} the standard error.  At the heart of the method lies an estimate of the false negative and false positive rates. The
test is calibrated so that these are minimised. For calculation of \eqn{p}-values a bootstrap procedure is invoked. Further details are given in Broberg (2003).
Note that the definition of p-values follows that in Davison and Hinkley (1997), in order to avoid p-values that equal zero.

The p-values are calculated through permuting the residuals obtained from the null model, assuming that this corresponds to the full model 
except for the parameter being tested, coresponding to the contrast 
coefficient not equal to zero. This means that factors not tested are kept fixed. NB This may be adequate for testing a factor with two levels or a
regression coefficient (correlation), but it is not adequate for all linear models.
 } 

\references{
Tusher, V.G., Tibshirani, R., and Chu, G. (2001) Significance analysis of microarrays applied to
the ionizing radiation response. \emph{PNAS} Vol. 98, no.9, pp. 5116-5121

Broberg, P. (2002) Ranking genes with respect to differential expression , \url{http://genomebiology.com/2002/3/9/preprint/0007}

Broberg. P: Statistical methods for ranking differentially expressed genes. Genome Biology 2003, 4:R41
\url{ http://genomebiology.com/2003/4/6/R41}

Davison A.C. and Hinkley D.V. (1997) Bootstrap Methods and Their Application. Cambridge University Press

}

\keyword{htest}
