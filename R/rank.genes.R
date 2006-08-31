# A PCA approach to ranking genes #
# It is assumed that genes come in rows and the criteria in columns #
# Furthermore, high values should be good #
rank.genes <- function(data = indats){
# data <- dats.min  #
   data <- apply(data, 2, rank)
   library(stats)
   tr.pca <- prcomp(data, retx=TRUE,center=FALSE,scale.= FALSE) # SVD DIRECTLY ON DATA MATRIX #
return(rank(tr.pca$x[,1]))
}
# The first PCA gives a sum up of the different ranking criteria #
