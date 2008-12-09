\name{fp.fn}
\alias{fp.fn}
\title{Calculation of fp and fn based on a vector of p-values}

\description{Based on a vector of p-values the proportion false positive (fp) and the
proportion false negative are calculated for each entry, assuming that one to be the last to
be called significant. The sum of fp and fn is also calculated (errors). Furthermore, an
estimate of the proportion unchanged together with the number of the entry with minimum errors. }



\usage{
fp.fn(ps = pvals, B = 100)
}
\arguments{
\item{ps}{a vector of p-values}
\item{B}{the number of bootstrap loops done by the function estimatep0 called by fp.fn} 
}
\value{
A list with components
\item{p0}{the estimated proportion unchamged}
\item{fp}{the estimated proportion false positives}
\item{fn}{the estimated proportion false negatives}
\item{N}{the number of the p-value (significance level) that gives minimum fp + fn} 

}
\author{Per Broberg}

\keyword{htest}
