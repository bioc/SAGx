
\name{JT.test}
\alias{JT.test}
\title{Jonckheere-Terpstra trend test}

\description{The test is testing for a monotone trend in terms of the class parameter.
             The number of times that an individual of a higher class has a higher gene expression forms a basis for the inference.}

\usage{trendA <- JT.test(data = datsA, class = g, labels = c("NS", "HS", "COPD0", "COPD1", "COPD2"))}

\arguments{
\item{data}{A matrix with genes in rows and subjects in columns}
\item{class}{the column labels}
\item{labels}{the labels of the categories coded by class}
}

\value{a matrix with the columns
\item{pvalue}{the p-value for the two-sided test of no trend.}
\item{median1 ... mediann}{the medians for the n groups}
\item{trend}{the rank correlation with category}
}

\details{Assumes that groups are coded in incresing numerical order}


\examples{
# Enter the data as a vector
A <- as.matrix(c(99,114,116,127,146,111, 125,143,148,157,133,139, 149, 160, 184))
# create the class labels
g <- c(rep(1,5),rep(2,5),rep(3,5))
# The groups have the medians
tapply(A, g, median)
# JT.test indicates that this trend is significant at the 5% level
JT.test(data = t(A), class = g, labels = c("GRP 1", "GRP 2", "GRP 3",))
}

\author{Per Broberg}

\references{
Lehmann, EH (1975) \emph{Nonparametrics: Statistical Methods Based on Ranks} p. 233. Holden Day}

\keyword{nonparametric}
