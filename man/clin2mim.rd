\name{clin2mim}
\alias{clin2mim}
\title{Output a script file to WinMIM, linking clinical data and gene expression}


\description{Given a clinical variable, it produces a script file for WinMIM by calculating means and covariances and 
for the N most highly correlated probes (in absolute value). Here N is an input parameter, but a recommended value 10. WinMIM can find
a relevant graphical model for the dependencies between the probes and the clinical variable.}
\usage{
clin2mim(variable="FEV1.ACTUAL",data=dbs,clindat=clinical,probes=probes,N=10,out="mimscr.txt")
}

\arguments{
\item{variable}{Clinical variable to be examined}
\item{data}{The input data set, with subject id in first column.}
\item{clindat}{The input clinical data, with subject id in first column}
\item{probes}{The name of the probes in the order of {\it data}}
\item{N} {The number of highly correlated probes to be studied} 
\item{out}{The MIM script file}
}

\value{
The correlation matrix
}

\author{Per Broberg}

\note{David Edwards' program WinMIM can be found on 
StatLib (\url{http://lib.stat.cmu.edu/graphmod/}). In MIM issue \code{input mimscript.txt}
and the calculations to find a model will start. When finished go to the Graphics menu and 
click on \code{Independence Graph}. The resulting graph can be exported both to WMF and LaTeX.}

\references{
Edwards, David (1995) \emph{Introduction to Graphical Modelling}. Springer-Verlag \cr
Lautitzen, Steffen (1996) \emph{Graphical Models}. Oxford University Press \cr
Whittaker, Joe (1990) \emph{Graphical Models in Multivariate Analysis}. Wiley \cr
}

\keyword{multivariate}
