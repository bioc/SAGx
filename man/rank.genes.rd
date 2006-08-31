
\name{rank.genes}
\alias{rank.genes}

\title{Rank genes with respect to multiple criteria}

\description{It is assumed that genes come in rows and the criteria in columns. 
             Furthermore, high values should be good.
 After ranking the genes with respect to each criterion, the function does a PCA on the ranks,
 uses the firsta PC to obtain the final ranks. In principle it could happen that genes are ranked
 in the opposite direction to the one intended, but that should be evident from a quick glance at the
 results.}

\usage{rank.genes(data = indats)}

\arguments{
\item{data}{A matrix with the criteria in columns.}
}

\author{Per Broberg}

\value{The total ranks of the genes.}

\keyword{multivariate}
