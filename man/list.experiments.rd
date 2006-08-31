\name{list.experiments}
\alias{list.experiments}
\title{Display all experiment names and id's}
\description{Display all experiment names and id's in the GATC database}
}

\usage{
list.experiments(channel, chip = "HG_U95Av2")
}
\arguments{
\item{channel}{the ODBC channel set up through RODBC}
\item{chip}{the chip type}
}
\value{dataframe with column EXPERIMENT 
}

\details{The GATC database has caused some problems by switching between upper and lower case in an erratic manner. To solve
this all names are changed to upper case in the identification of experiments. Thus the function will not distinguish between the
experiments 'A' and 'a', but with any sensible naming strategy, the restriction is without consequence}

\examples{
# Not run
\dontrun{library(Rodbc)
channel <- odbcConnect(DBN, USRID, PWD) 
ut <- list.experiments(channel, chip = "hu6800")
colnames(ut)
#[1] "EXPERIMENT"  
}
}

\keyword{database}
