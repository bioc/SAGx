\name{firstpass}
\alias{firstpass}
\title{First pass description of GeneChip data}

\description{Does a first-pass analysis for a comparative experiment. This
includes the calculation of means and confidence intervals for the groups, and finally
a Kruskal-Wallis p-value for the null hypothesis of no difference }

\usage{
firstpass(data = D, probes = probes , g, log = FALSE, present = NULL, labels = NULL, output.data = FALSE)
}
\arguments{
\item{data}{A data frame with one array in each column}
\item{probes}{a vector containing the names of the probes in the same order as rows in D} 
\item{g}{A vector with the groups for the arrays, eg. TREATMENT and CONTROL}
\item{present}{A dataframe with the Present calls, 3 = P, 2 = M, 1 = A.}
\item{log}{if TRUE then data are log transformed through t(x) = log(1+x) and geometric means are calculated}
\item{labels}{a vector of labels given the group means}
\item{output.data}{if T the raw data are included in the output}
}
\value{
A dataframe with the coumns PROBES, followed by group means and sd's, lower confidence intervals and then,
 upper  confidence interval (confidence level 95\%), and followed a Kruskal-Wallis p-value, and finally  the input data,. If present names a dataframe holding
the present calls the proportion present is calculated. Furthermore, if there are two groups
the difference in group means is added.
}
\details{ A speed-up for Wilcoxon based on Kronecker products was put in place with SAGx v.1.4.5. Ties are currently not taken into account in Wilcoxon.}
\examples{
\dontrun{
# not run
 g <- c(rep(1,4),rep(2,4)); labs <- c("Mean Diet","Mean Control"); probes <- paste("Probe",1:1000)
 firstpass(data = utmat[1:2,], probes = probes[1:2], g, log = FALSE, labels = labs)
#  Probesets         Mean Diet      Mean Control             LCL.1             LCL.2             UCL.1             UCL.2               pval
#1   Probe 1 -12.3444460036497 -11.7495704973055 -12.9047961446666 -12.2832657957485 -11.7840958626327 -11.2158751988625 0.0433081428107922
#2   Probe 2 -7.99773926405627 -8.02799133391929 -8.47704512876227 -8.19487551919835 -7.51843339935028 -7.86110714864023  0.772829992684449
#          Difference  Subject 1  Subject 2  Subject 3  Subject 4  Subject 5  Subject 6  Subject 7  Subject 8
#1 -0.594875506344176 -12.345150 -11.805071 -12.776232 -12.451332 -11.595748 -12.320430 -11.482349 -11.599755
#2 0.0302520698630131  -7.660097  -8.157944  -8.404433  -7.768484  -7.979951  -8.017327  -8.197361  -7.917326
}
}

\keyword{nonparametric}
