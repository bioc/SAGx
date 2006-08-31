\name{union.of.pways}
\alias{union.of.pways}

\title{Create the union of two pathway lists}

\description{This function takes two lists where each component is a vector of probe sets ids and create
a new such list that contains all probe sets and pathways from the two lists.}

\usage{union.of.pways(x,y)}

\arguments{
\item{x}{the first list}
\item{y}{the second list}
}

\author{Per Broberg}

\value{A list which is the union of the two input lists.}

\examples{
X = list(a=c(1,2),c=c(1,2));Y = list(a=c(3,4),d=c(12,2))
union.of.pways(X,Y)
}

\note{Experimental function. Feedback appreciated.}

\details{The function \emph{merge.list} in package \emph{RCurl} forms a basis for this function which adds the 
ability to add new probe sets to existing pathways.}


\keyword{multivariate}
