\name{one.probeset.per.gene}
\alias{one.probeset.per.gene}

\title{Select the best probeset per gene}

\description{This function takes a vector of probeset identifiers, a vector of gene identifiers and a vector of present rates,
and outputs the probeset id per gene that corresponds to the highest present rate.}

\usage{one.probeset.per.gene(probeset = probeset, present = present, symbol = symbol)}

\arguments{
\item{probeset}{a vector of probeset id's}
\item{present}{a vector of present rates}
\item{symbol}{a vector of gene symbols}
}

\author{Per Broberg}

\value{A vector of probeset id's.}

\details{It is assumed that missing gene symbol is coded as ""}

\note{Experimental function. Feedback appreciated.}


\keyword{multivariate}
