\name{GSEA.mean.t}
\alias{GSEA.mean.t}

\title{Gene Set Enrichment Analysis using output from samroc}

\description{Based on a list of gene sets, e.g. pathways, in terms Affymtrix identifiers, these
sets are ranked with respect to regulation as measured by an effect in a linear model using the SAM statistic.}

\usage{GSEA.mean.t(data = X, samroc = samroc.res, probeset = probeset, 
pway = kegg, absolute = TRUE, two.side = FALSE, cutoff = c(5,Inf), B = 1000, smooth = FALSE)}

\arguments{
\item{data}{data matrix}
\item{samroc}{an object of class samroc.result}
\item{probeset}{the Affymetrix identifiers}
\item{pway}{a list of pathways or gene sets}
\item{absolute}{if TRUE the absolute value of the absolute value of the samroc test statistic is used. Currently this is
the only option}
\item{two.side}{if TRUE a two-sided test is performed. Currently only one-sided test allowd}
\item{cutoff}{Gene sets with the number of members not falling within the interval given by \emph{cutoff} are excluded}
\item{B}{the number of permutations}
\item{smooth}{if TRUE smoothing of SEs is performed. Currently this is not permitted}
}

\author{Per Broberg}

\value{A vector of p-values}

\note{Experimental function}

\references{Tian, Lu and Greenberg, Steven A. and Kong, Sek Won and Altschuler, Josiah and Kohane, Isaac S. and Park, Peter J. 
(2005) Discovering statistically significant pathways in expression profiling studies, \emph{PNAS} Vol. 102, nr. 38, pp. 13544-13549

}


\keyword{multivariate}
