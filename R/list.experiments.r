##### Per Broberg 26JAN01. Modified 08FEB02 ##############

list.experiments <- function(channel, chip="HG_U95Av2"){
   chip <- toupper(chip)
t1 <- "select q3.name as EXPNAME from chip_design q1, physical_chip q2,"
t2 <-  paste("experiment q3 where upper(Q1.NAME) =","'",chip,"'",sep="")
t3 <- "and q2.design_id=q1.id and q3.physical_chip_id=q2.id"
t <- paste(t1,t2,t3)
sqlQuery(channel,t)
}


