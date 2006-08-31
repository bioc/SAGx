\name{samroc.result-class}
\docType{class}
\alias{samroc.result-class}
\alias{show,samroc.result-method}
\alias{plot,samroc.result-method}


\title{Class "samroc.result" for results of the function samrocN}

\description{The class samroc.result is the output of a call to
\code{\link{samrocN}} and the input of various other
functions. }

\section{Slots}{
  \describe{
    \item{\code{di}:}{Object of class "numeric". Observed test statistic.}
    \item{\code{diff}:}{Object of class "numeric". Estimate of effect, e.g. difference between group means.}
    \item{\code{se}:}{Object of class "numeric". Standard error of \code{diff}.}
    \item{\code{d0}:}{Object of class "matrix". Permutation test statistics.}
    \item{\code{p0}:}{Object of class "numeric". The estimated proportion unaffceted genes.}
    \item{\code{s0}:}{Object of class "numeric". The fudge factor.}
    \item{\code{pvalues}:}{Object of class "numeric". The p-values.}
    \item{\code{N.list}:}{Object of class "integer". The optimal top list size among the sizes suggested.}
    \item{\code{errors}:}{Object of class "numeric". The sum of false postives and false negatives given a list that includes the
    current gene.}
    \item{\code{formula}:}{Obeject of class "formula". The linear model formula used.}
    \item{\code{contrast}:}{Object of class "numeric". The contrast estimated.}
    \item{\code{annotation}:}{Object of class "character". Annotation or comments regarding the analysis. By default the date.} 
    \item{\code{N.sample}:}{Object of class "integer". The number of samples.}
    \item{\code{B}:}{Object of class "integer". The number of premutations.}
    \item{\code{call}:}{Object of class "character". The call to the function.} 
    \item{\code{id}:}{Object of class "character". The probeset ids.} 
    \item{\code{error.df}:}{Object of class "integer". The error degrees of freedom.} 
    }
}


\section{Methods}{
  \describe{
    \item{show}{(samroc.result): Summarizes the test result.}
    \item{plot}{(samroc.result): Plots the density of the observed test statistic and that of the corresponding null distribution}
    }
}



\author{Per Broberg}

\seealso{\code{\link{samrocN}}}

\keyword{methods}
