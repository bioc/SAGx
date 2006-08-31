### Per Broberg 25JAN01 ###############################
# History : chip design criterion added 25SEP01;updated 09OCT03 #

fetchSignal <- function (experiment = "AZ33 ALI", channel, chip="HG_U95Av2") {
    chip <- toupper(chip)
    t1 <- "SELECT q2.NAME as FILENAME,q3.ITEM_NAME  as PROBESET,q1.SIGNAL as Signal,"
    t2 <- "q4.name as abs_call, q5.name as CHIPNAME FROM ABS_GENE_EXPR_RES_STAT Q1,BIOLOGICAL_ITEM Q3,"
    t3 <- "ANALYSIS Q2,ABS_GENE_EXPR_RESULT_TYPE Q4, CHIP_DESIGN Q5, ANALYSIS_SCHEME q6"
    t4 <- paste("WHERE Q3.ID = Q1.ITEM_ID  AND Q2.ID = Q1.ANALYSIS_ID  AND upper(Q5.NAME) =","'",chip,"'",sep="")
    t5 <- paste("AND Q4.ID = Q1.TYPE_ID AND upper(Q2.NAME) LIKE ", "'", experiment, "'", sep = "")
    t6 <- paste("and q5.id=q6.chip_design_id and q6.id=q2.scheme_id")
    test <- paste(t1, t2, t3, t4,t5,t6)
    return(sqlQuery(channel, paste(test)))
}
