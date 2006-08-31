\name{estimatep0}
\alias{estimatep0}

\title{Estimate proportion unchanged genes}

\description{The function uses the vector of p-values to estimate p0.}

\usage{estimatep0(ps = pp, B = 500, range = seq(0,0.95, by = 0.05))}

\arguments{
\item{ps}{the vector of p-values, e.g. from firstpass}
\item{B}{the number of Bootstrap samples}
\item{range}{the values considered}
}

\author{Per Broberg}

\value{the value of p0, the proportion unchanged  genes}



\references{
Storey, J. A Direct Approach to the False Discovery Rate, Technical Report Stanford (2001)
}

\keyword{htest}
