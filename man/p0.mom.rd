\name{p0.mom}
\alias{p0.mom}

\title{Estimate proportion unchanged genes}

\description{The function uses the vector of p-values to estimate p0.}

\usage{p0.mom(ps = pvalues)}

\arguments{
\item{ps}{the vector of p-values, e.g. from firstpass}
}

\value{the value of p0, the proportion unchanged  genes as a list with
components 
\item{mgf}{estimate from the mgf method}
\item{PRE}{estimate from the PRE method}
\item{experimental1}{ }
\item{experimental2}{ }
}


\author{Per Broberg}

\references{
Broberg, P. A new estimate of the proportion unchanged genes, 2005, \emph{Genome Biology} 5:p10 \cr
Broberg, P. A comparative review of estimates of the proportion unchanged genes and the false discovery rate, submitted (2004)
}

\keyword{nonparametric}
