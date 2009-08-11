\name{rank.trend}
\alias{rank.trend}
\title{Trend analysis based on ranks}

\description{Ranks are used to score genes with respect to degree of agreement to a given trend or pattern, Lehmann (1974) p.294.}

\usage{
rank.trend(data = x, pattern = c(1:ncol(data)), har = FALSE)
}
\arguments{
\item{data}{A data frame with one array in each column}
\item{pattern}{A permutation of the integers 1:ncol(data)} 
\item{har}{logical parameter indicating whether or not a score based on Hardy's theorem shall be calculated.}
}
\value{
A list with the components
\item{score}{the rank score for each gene}
\item{hardy}{if har = TRUE the hardy score, NULL otherwise}
\item{pvals}{the p-values for the null hypothesis of no trend}
}
\details{The rank scores gives a higher weight to a deviation from trend in more distant obseveations than a deviation between neighbouring observations.
The p-values are calculated through a normal approximation.}
\references{ Lehmann, E.L. (1975) Nonparametrics: Statistical Methods Based on Ranks, Holden-Day}

\author{Per Broberg}

\examples{
# not run
D <- c(123, 334, 578, 762, 755, 890)
rank.trend(data = t(as.matrix(D)), har = TRUE)
#     Trend score Hardy score p-value for no trend
# [1,]           2          90           0.01750284


}

\keyword{robust}
