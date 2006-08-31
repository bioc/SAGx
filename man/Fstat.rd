\name{Fstat}
\alias{Fstat}

\title{Calculation of F statistic by gene given a linear model}

\description{Calculates F statistic.}

\usage{Fstat(indata =  M, formula1 = ~as.factor(g), formula0 = "mean", design1 = NULL, design0 = NULL, B = NULL)}

\arguments{
\item{indata}{The data matrix}
\item{formula1}{a formula descibing the alternative linear model}
\item{formula0}{a formula describing the nullmodel. Use linear models syntax, except for one-way ANOVA ("mean")}
\item{design1}{the alternaive design matrix. If not NULL it overrides the formula argument}
\item{design0}{the null design matrix. If not NULL it overrides the formula argument}
\item{B}{the number of bootstrap replicates}
}

\author{Per Broberg}

\value{A list with the components
\item{Fstat}{the value of the F statistic}
\item{fnum}{the numerator degrees of freedom}
\item{fdenom}{the denominator degrees og freedom}
\item{design1}{the alternative design matrix}
\item{design0}{the null design matrix}
\item{SS1}{the sum of squares in the denominator of the F-statistic}
\item{SS0}{the sum of squares in the numerator  of the F-statistic} 
\item{pvalue}{the p-value for testing the alternative vs the null model}
}


\examples{
## Annette Dobson (1990) "An Introduction to Generalized Linear Models". 
## Page 9: Plant Weight Data. 
 ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14) 
 trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69) 
 group <- gl(2,10,20, labels=c("Ctl","Trt")) 
 weight <- c(ctl, trt) 
 anova(lm.D9 <- lm(weight ~ group)) 
# Analysis of Variance Table

# Response: weight
#          Df Sum Sq Mean Sq F value Pr(>F)
#group      1 0.6882  0.6882  1.4191  0.249
#Residuals 18 8.7292  0.4850               

 Fstat(indata = rbind(weight,weight),formula1=~group) # Fstat will need at least two genes to work with #
#$Fstat
#  weight   weight 
#1.419101 1.419101 

#$fnum
#[1] 18

#$fdenom
#[1] 1

#$design1
#   (Intercept) groupTrt
#1            1        0
#2            1        0
#3            1        0
#4            1        0
#5            1        0
#6            1        0
#7            1        0
#8            1        0
#9            1        0
#10           1        0
#11           1        1
#12           1        1
#13           1        1
#14           1        1
#15           1        1
#16           1        1
#17           1        1
#18           1        1
#19           1        1
#20           1        1
#attr(,"assign")
#[1] 0 1

# $design0
# NULL

# $SS1
# weight  weight 
#8.72925 8.72925 

#$SS0
#  weight   weight 
#0.688205 0.688205 


}

\keyword{models}
