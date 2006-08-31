\name{fetchSignal}
\alias{fetchSignal}
\title{Fetch data from the GATC database}
\description{Fetch FILENAME, PROBESET, SIGNAL and ABS\_CALL from the GATC database}
}

\usage{
fetchSignal(experiment="AZ33 ALI", channel, chip="HG_U95Av2")
}
\arguments{
\item{experiment}{The name of the experiment corresponding to an individual chip}
\item{channel}{The channel to the database}
\item{chip}{the chip type}
}
\value{dataframe with columns 
}
\examples{
\dontrun{
# Do not run example 1. Fetch Probeset, Signal, ABS_CALL and CHIP for one sample.
library(RODBC)
(channel<-odbcConnect("DSN",uid="USERID",pwd="PASSWORD"))
ali.data <-fetchSignal(experiment="AZ33 ALI", channel, chip="hg_u95a")
 colnames(ali.data)
#[1] "FILENAME" "PROBESET" "SIGNAL" "ABS_CALL" "CHIP"

# Do not run example 2
t1 <- paste("select q1.name as name from experiment q1, physical_chip q2, chip_design q3")
t2 <- paste("where q1.physical_chip_id=q2.id and q3.id=q2.design_id and ")
t3 <- paste("upper(q1.name) like '%ALI%' and upper(q3.name)='HG_U95A'") 
Ids <- sqlQuery(channel,paste(t1,t2,t3) )   
# fetch Signal from GATC corresponding to the U95A chip for all samples in experiment. #
tmp <- apply(Ids,1,toupper)
probes <- data.frame(fetchSignal(experiment=tmp[1],channel, chip="hg_u95a")[,"PROBESET"])
test <- matrix(nrow=nrow(as.data.frame(probes)),ncol=nrow(Ids))
for(i in 1:nrow(as.data.frame(tmp))){
   test[,i] <- fetchSignal(experiment=tmp[i],channel, chip="hg_u95a")[,"SIGNAL"]
}
codes <- data.frame(apply(Ids,1,code<-function(x) substr(x,1,5)))
colnames(test) <- as.character(t(codes))
test <- test[,order(colnames(test))]
}
}

\author{Ported to R by Per Broberg. Original Oracle code by Petter Hallgren, with input from Petra Johansson.}
\keyword{database}
