\name{pava.fdr}
\alias{pava.fdr}

\title{Estimate of the FDR and the proportion unchanged genes}

\description{Estimates tail area and local false discovery rate using isotonic regression}

\usage{pava.fdr(ps = pvalues, p0 = NULL)}

\arguments{
\item{ps}{the vector of p-values, e.g. from firstpass}
\item{p0}{an estimate of the proportion unchanged genes}
}

\author{Per Broberg}

\value{a list with components
\item{pava.fdr}{estimate of the FDR}
\item{p0}{estimate of p0}
\item{pava.local.fdr}{estimate of the local fdr}
}

\details{If p0 = NULL the PRE estimate of p0 is calculated.}

\keyword{htest}

\references{
Broberg, P : A comparative review of estimates of the proportion unchanged genes and the false discovery rate, \emph{BMC Bioinformatics} 2005, 5(1):199 \cr
Aubert J, Bar-Hen A, Daudin J-J, Robin S: Determination of the differentially expressed genes in microarray experiments using local FDR. \emph{BMC Bioinformatics} 2004, 6(1):125
}
