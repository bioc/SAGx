\name{Xprep.resid}
\alias{Xprep.resid}

\title{Calculation of input of residuals from linear model}

\description{The function fits a linear model to a microarray data matrix and calculates the residuals.}

\usage{Xprep.resid(data=M, formula=~as.factor(g), design=NULL)}

\arguments{
\item{data}{The data matrix}
\item{formula}{a linear model formula in the lm format}
\item{design}{the design matrix}
}

\author{Per Broberg}

\value{
A matrix with the residuals
}


\keyword{models}
