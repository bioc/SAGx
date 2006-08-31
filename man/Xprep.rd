\name{Xprep}
\alias{Xprep}

\title{Fitting of a linear model}

\description{The function fits a linear model to a microarray data matrix.}

\usage{Xprep(indata=M, formula=~as.factor(g), contrast=c(0,1), design=NULL)}

\arguments{
\item{indata}{The data matrix}
\item{formula}{a linear model formula in the lm format}
\item{contrast}{a vector defining the contrast of interest}
\item{design}{the design matrix}
}

\author{Per Broberg}

\value{a list with the entries 
\item{Mbar}{estimate of the contrast}
\item{Vest}{the error variance}
\item{k}{inverse of the scale factor turning Vest into a standard error}
\item{f}{the degrees of freedom of Vest}
\item{design}{the design matrix}
}


\keyword{models}
