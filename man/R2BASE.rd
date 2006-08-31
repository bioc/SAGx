\name{R2BASE}
\alias{R2BASE}

\title{Produces a BASE file}

\description{The function produces a BASE file for import to Gene Data Viewer.}

\usage{R2BASE(context.data = clingen, sample.ids = AZID, expression.data = dats, 
annotation = annots, out = "u:\\temp\\temp.base")}

\arguments{
\item{context.data}{e.g. a clinical database}
\item{sample.ids}{Sample Ids, that names the columns of the expression data.}
\item{expression.data}{a matrix with the gene expression data, samples correspond to columns and probesets to rows. It is assumed that probeset identifiers are found in the first column.}
\item{annotation}{annotations of the probesets, i.e. the rows in the expression.data. It is assumed that probeset identifiers are found in the first column.}
\item{out}{the output file including path}
}

\author{Per Broberg}

\value{The file produced complies with an old BASE format. However, none of these formats are documented
, as far as I know. So, essentially this function defines a data format that can be read by e.g. Gene Data Viewer.}

\keyword{models}
