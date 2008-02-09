

\name{JT.test}
\alias{JT.test}
\title{Jonckheere-Terpstra trend test}

\description{The test is testing for a monotone trend in terms of the class parameter.
             The number of times that an individual of a higher class has a higher gene expression forms a basis for the inference.}

\usage{trendA <- JT.test(data, class, labs = c("NS", "HS", "COPD0", "COPD1", "COPD2"), alternative = c("increasing", "decreasing", "two-sided"))}

\arguments{
\item{data}{A matrix with genes in rows and subjects in columns}
\item{class}{the column labels, if not an ordered fctor it will be redefined to be one.}
\item{labs}{the labels of the categories coded by class}
}

\value{an object of class JT-test, which extends the class htest, and includes the following slots
\item{statistic}{the observed JT statistic}
\item{parameter}{the null hypothesis parameter, if other value than 0.}
\item{p.value}{the p-value for the two-sided test of no trend.}
\item{method}{Jonckheere-Terpstra}
\item{alternative}{The relations between the levels: decreasing, increasing or two-sided}
\item{data.name}{the name of the input data}
\item{median1 ... mediann}{the medians for the n groups}
\item{trend}{the rank correlation with category}
\item{S1}{Predictive strength}
}

\details{Assumes that groups are given in increasing order, if the class variable is not an ordered factor, it will be redefined to be one.

The implementation owes to suggestions posted by to R list. The definition of predictive strength appears in Flandre and O'Quigley.}


\examples{
# Enter the data as a vector
A <- as.matrix(c(99,114,116,127,146,111, 125,143,148,157,133,139, 149, 160, 184))
# create the class labels
g <- c(rep(1,5),rep(2,5),rep(3,5))
# The groups have the medians
tapply(A, g, median)
# JT.test indicates that this trend is significant at the 5% level
JT.test(data = A, class = g, labs = c("GRP 1", "GRP 2", "GRP 3"), alternative = "two-sided")
}

\author{Per Broberg, acknowledging input from Christopher Andrews at SUNY Buffalo}

\references{
Lehmann, EH (1975) \emph{Nonparametrics: Statistical Methods Based on Ranks} p. 233. Holden Day}
Flandre, Philippe and O'Quigley, John, \emph{Predictive strength of Jonckheere's test for trend: an application to genotypic scores in HIV infection},
Statistics in Medicine, 2007, 26, 24, 4441-4454  

\keyword{nonparametric}
