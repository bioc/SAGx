
\name{rsd.test}
\alias{rsd.test}

\title{Compare two groups with respect to their RSD (CV)}

\description{A by row comparison of the Relative Standard Deviation (RSD), asa Coefficient of Variation (CV), is done
using a bootstrap}

\usage{rsd.test(data1 = x, data2= y, B = NULL)}

\arguments{
\item{data1}{A matrix with the samples for group 1  in columns.}
\item{data2}{A matrix with the samples for group 2  in columns.}
\item{B}{the number of bootstrap iterations. If NULL no bootstrap is performed.}
}

\value{
A list with the components
\item{cv1}{A vector of the RSD's for sample 1}
\item{cv2}{A vector of the RSD's for sample 2}
\item{t.stat}{the test statistic}
\item{p.vals}{A vector of p-values for the comparison between \eqn{cv1} and \eqn{cv2}}

}
\author{Per Broberg}

\references{
Broberg P, Estimation of Relative Standard Deviation,(1999) in \emph{Drug Development and Industrial Pharmacy}, Vol 25 no 1 37-43
}

\keyword{nonparametric}
