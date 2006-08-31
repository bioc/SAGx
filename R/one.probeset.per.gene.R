# Author : PBG/17NOV05
# Description : selects the probeset with the highest present rate for each gene
# Detail : it is assumed that missing symbol is coded as "".
one.probeset.per.gene <-function(probeset = probeset, present = present, symbol = symbol){ 
   probeset.sel = probeset[order(symbol,-present)]
   symbol = symbol[order(symbol,-present)]
   probeset.sel = probeset.sel[!symbol==""]
   symbol = symbol[!symbol==""]
   indices = match(unique(symbol), symbol)
   probeset.sel = probeset.sel[indices]
   probeset.sel
}
