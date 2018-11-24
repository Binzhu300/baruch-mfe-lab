# We construct a piecewise constant forward variance curve from variance swap data

xi.curve <- function(expiries,w) function (t) 
{
  n <- length(w)
  xi.vec <- c(w[1]/expiries[1], diff(w)/diff(expiries))
  exp.vec <- c(0,expiries)
  exp.vec[n+1] <- Inf # Long expiries are all in the last bucket
  pos <- findInterval(t, vec = exp.vec)
  res <- xi.vec[pos]
  return(res)
}