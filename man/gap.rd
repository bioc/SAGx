\name{gap}
\alias{gap}

\title{GAP statistic clustering figure of merit}

\description{Calculates a goodness of clustering measure based on the average 
dispersion compared to a reference distribution.}

\usage{gap(data = swiss,class = g, B = 500)}

\arguments{
\item{data}{The data matrix}
\item{class}{a vector descibing the cluster memberships of the rows of data}
\item{B}{the number of bootstrap samples}
}

\author{Per Broberg}

\value{The GAP statistic}

\examples{
library(MASS)
data(swiss)
cl <- myclus(data = swiss, k = 3)
gap(swiss,cl$cluster)
}

\references{
Tishirani, R., Walther, G. and Hastie, T. (2000) Estimating the number of clusters in a dataset
via the Gap statistic. \emph{ Technical Report} Stanford
}

\keyword{multivariate}
