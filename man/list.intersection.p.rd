\name{list.intersection.p}
\alias{list.intersection.p}
\title{p-value for intersection of two gene lists.}

\description{Calculates a p-value for observing a number of probe sets common to two lists drawn
from the same chip.
}

\usage{list.intersection.p(N = 14000, N1 = 100, N2 = 200, common = 30)}

\arguments{
\item{N}{The selectable number of probe sets}
\item{N1}{the number of probe sets on the first list.}
\item{N2}{the number of probe sets on the second list}
\item{common}{the number of probe sets in common to the two lists.}
}

\author{Per Broberg}

\value{the p-value giving the probability of observing by chance at least as many in common as was actually observed.}

\keyword{distribution}
