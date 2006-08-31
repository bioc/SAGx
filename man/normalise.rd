\name{normalise}
\alias{normalise}

\title{Normalise arrays}

\description{Normalises arrays against a calculated average array, and calibrated 
linearly in a cube-root scatter plot.}

\usage{normalise(x,linear=TRUE)}

\arguments{
\item{x}{The data matrix}
\item{linear}{if linear=TRUE then the matrix elements are raised to the power of 3.}
}

\author{Per Broberg}

\value{normalised version of indata}

\references{
Tusher, V.G., Tibshirani, R., and Chu, G. (2001) Significance analysis of microarrays applied to
the ionizing radiation response. \emph{PNAS} Vol. 98, no.9, pp. 5116-5121
}

\keyword{multivariate}
