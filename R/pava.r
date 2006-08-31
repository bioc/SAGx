# Pool adjacent violators #
pava <- function (x, wt=rep(1,length(x)))
#  Compute the isotonic regression of numeric vector 'x', with 
#  weights 'wt', with respect to simple order.  The pool-adjacent-
#  violators algorithm is used.  Returns a vector of the same length
#  as 'x' containing the regression.

#  02 Sep 1994 / R.F. Raubertas
{
   n <- length(x)
   if (n <= 1) return (x)
   if (any(is.na(x)) || any(is.na(wt))) {
      stop ("Missing values in 'x' or 'wt' not allowed")
   }
   lvlsets <- (1:n)
   repeat {
      viol <- (as.vector(diff(x)) < 0)  # Find adjacent violators
      if (!(any(viol))) break

      i <- min( (1:(n-1))[viol])        # Pool first pair of violators
      lvl1 <- lvlsets[i]
      lvl2 <- lvlsets[i+1]
      ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
      x[ilvl] <- sum(x[ilvl]*wt[ilvl]) / sum(wt[ilvl])
      lvlsets[ilvl] <- lvl1
   }
   x
}
