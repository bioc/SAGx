\name{outlier}
\alias{outlier}
\title{Identify outliers in the multivariate distribution}
\description{A PCA model is fitted to data and two statistics as measures of extremity are calculated. These are the Hotelling t-square and DMODX, the first is a measure of how far away from the centre of the projection subspace the projection of the observation is. The second one measures how remote from the projection the actual observation is. SVD is done directly on the data matrix. The number of significant dimensions is defined as the number of eigenvalues greater than 1. Typically arrays are in different columns.}

\usage{
outlier(M)
}
\arguments{
\item{M}{matrix}
}
\value{Dataframe with columns Hotelling and DMODX}
\author{Per Broberg}
\references{Jackson, J.E. (1991) \emph{A User's Guide to Principal Components}. Wiley}
\examples{
\dontrun{
# not run
ut<-outlier(M)
#[1] "The number of significant dimensions is 19"
colnames(ut)
#[1] "Hotelling" "DMODX"  
}
}

\keyword{multivariate}
