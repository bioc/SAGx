\name{R2mim}
\alias{R2mim}
\title{Output a script file to WinMIM}


\description{Given a candidate probe, it produces a script file for WinMIM by calculating means and covariances and 
for the N most highly correlated probes (in absolute value). Here N is an input parameter, but a recommended value 10. WinMIM can find
a relevant graphical model for the dependencies between the probes.}

\usage{
R2mim(probe="12345_at",N=10,data=inm,out="u:\\study\\copd\\mimscr.txt")
}
\arguments{
\item{probe}{The name of the candidate probe}
\item{N}{The number of highly correlated probes to be studied} 
\item{data}{The input data set}
\item{out}{The MIM script file}
}
\value{
The correlation matrix
}
\author{Per Broberg}
\note{David Edwards' program WinMIM can be found on 
StatLib (\url{http://lib.stat.cmu.edu/graphmod/}). In MIM issue \code{input mimscr.txt}
and the calculations to find a model will start. When finished go to the Graphics menu and 
click on \code{Independence Graph}. The resulting graph can be exported both to WMF and LaTeX.}
\references{
Edwards, David (1995) \emph{Introduction to Graphical Modelling}. Springer-Verlag \cr
Lauritzen, Steffen (1996) \emph{Graphical Models}. Oxford University Press \cr
Whittaker, Joe (1990) \emph{Graphical Models in Multivariate Analysis}. Wiley
}

\keyword{multivariate}
