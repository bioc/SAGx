\name{Zfreq}
\alias{Zfreq}

\title{Bayesian analysis of a Comparative experiment with two groups}

\description{The posterior probabilty of being differentially expressed is calculated as a 
function of the statistic Z from etgc1. A graph is output to the default, and a dataframe with
the probailities}


\usage{ut<-cbind(probes,Zfreq(z,Z,B=50))}

\arguments{
\item{z}{the simulated null distribution from etgc2}
\item{Z}{the the statistic output from etgc1}
}

\author{Per Broberg}

\value{A dataframe with the columns PROBE, and P, a Bayesian a-posteriori probability of being expressed}

\keyword{htest}
