\name{fom}
\alias{fom}

\title{Clustering Figure of Merit}

\description{Goodness of clustering measure based on prediction error.}

\usage{fom(data,cluster)}

\arguments{
\item{data}{The data matrix}
\item{cluster}{a vector descibing the cluster memberships}
}

\author{Per Broberg}

\value{The Figure of Merit measure of the current clustering}

\details{The criterion in the Reference is not correct in the article (i.e. does not follow from the premises), but has been corrected here.}

\references{
Yeung, K.Y., Haynor, D.R. and Ruzzo, W.L. (2001) Validating clustering for gene expression data. 
\emph{Bioinformatics} Vol. 17, pp. 309-318
}

\keyword{multivariate}
