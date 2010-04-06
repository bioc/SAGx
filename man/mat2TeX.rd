\name{mat2TeX}

\alias{mat2TeX}

\title{Ouput matrix to LaTeX}

\description{The function outputs a matrix to a LaTeX table}

\usage{mat2TeX(mat, digits = 4, rowNameTitle = "", file = "",
 roundNum = NULL, rowNameAlign = "l", matAlign = "r",
 prtHead = TRUE, prtEnd = TRUE, extraTitle = NULL,
 rowNameCols = 1, append = FALSE)}

\arguments{
\item{mat}{a matrix}
\item{digits}{number of digits}
\item{rowNameTitle}{title above row names}
\item{file}{output file}
\item{roundNum}{integer indicating the precision}
\item{rowNameAlign}{alignment of row names, default is "l"}
\item{matAlign}{alignment of columns, default is "r"}
\item{prtHead}{if TRUE the begin\{tabular\} line is produced}
\item{prtEnd}{if TRUE the end\{tabular\} line is produced}
\item{extraTitle}{extra title}
\item{rowNameCols}{the row name column, default is 1}
\item{append}{if TRUE the output is appended to file, deafult is FALSE}


}

\author{Juerg Kindermann; code found on R list}

\keyword{IO}
