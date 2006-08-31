\name{pava}
\alias{pava}

\title{Pooling of Adjacent Violators}

\description{The PAVA algorithm}


\usage{pava(x, wt = rep(1, length(x)))}

\arguments{
\item{x}{A numeric sequence}
\item{wt}{observation weights; 1 by default.}
}

\details{The algorithm will turn a non-increasing into a non-decreasing one. pava is an internal function used to force monotonicity, e.g. of p1 in function
Zfreq}

\author{R.F. Raubertas, code from S list}

\value{A non-decreasing sequence}

\keyword{regression}

\examples{
pava(c(1,2,4,3,5))
# [1] 1.0 2.0 3.5 3.5 5.0
}
