
# Original code by Bill Venables ? #
mat2TeX <- function(mat, digits=4, rowNameTitle="", file="",
                    roundNum=NULL,
                    rowNameAlign="l",matAlign="r", prtHead=TRUE, prtEnd=TRUE,
                    extraTitle=NULL, rowNameCols=1, append=FALSE) {
  nr <- nrow(mat);      nc <- ncol(mat)
  rn <- rownames(mat);  cn <- colnames(mat)
  cmat <- array("",dim=dim(mat))
  if (!is.null(roundNum)) {
    if (length(roundNum)==nc){
      for(i in 1:nc) mat[,i] <- round(mat[,i],roundNum[i])
    }  else {
      errif(length(roundNum!=1))
      mat <- round(mat,roundNum)
    }
  }
  for(i in 1:nc){
    cmat[,i] <- format(mat[,i], digits=digits)
  }

  if(is.null(cn))   cn <- 1:nc
  if(is.null(rn))   rn <- 1:nr
  if(rowNameTitle == "") 
    rowNameTitle <- paste(rep(" ",rowNameCols), collapse="|")

  cat("\n ", file=file, append=append)  # write new if necessary 

  if (length(rowNameAlign)==1) rowNameAlign <- rep(rowNameAlign,rowNameCols)
  if (length(matAlign)==1) matAlign <- rep(matAlign,nc)

  if (prtHead){
    cat("\\begin{tabular}[c]{",
        paste(c("",rowNameAlign, matAlign,""), collapse="|"),
        "}\\hline\n", sep="",
        file=file, append=TRUE)
    if (!is.null(extraTitle))
      cat(extraTitle,"\\\\ \n",file=file, append=TRUE)
    cat(  paste(c(rowNameTitle,cn), collapse=" & "),"\\\\ \\hline\n",
        file=file, append=TRUE)
  }
  for(j in 1:nr){
    cat(paste(c(rn[j], cmat[j,]), collapse=" & "),
        "\\\\\n", file=file, append=TRUE)
  }
  cat("\\hline\n", file=file, append=TRUE)
  if (prtEnd) cat("\\end{tabular}\n", file=file, append=TRUE)
}
