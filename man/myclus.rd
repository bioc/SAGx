\name{myclus}
\alias{myclus}

\title{A clustering function}

\description{Uses a hierarchical clustering to initiate a kmeans clustering.}

\usage{myclus(data = swiss, k = 3)}

\arguments{
\item{data}{The data matrix}
\item{k}{the number of clusters}
}

\author{From Ripley and Venables}

\value{a list from function kmeans}

\examples{
library(MASS)
data(swiss)
cl <- myclus(data = swiss, k = 3)
gap(swiss,cl$cluster)
}

\references{
Venables, W.N. and Ripley, B.D (2000) \emph{Modern Applied Statistics with S-PLUS}, Springer
}

\keyword{multivariate}
