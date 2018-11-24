library(stinepack);


plotIvols <- function (ivolData, sviMatrix=NULL, modelVol=FALSE, slices = NULL) 
{
  bidVols <- as.numeric(ivolData$Bid)
  askVols <- as.numeric(ivolData$Ask)
  if(modelVol){modelVols <- as.numeric(ivolData$modelVol)}
  expDates <- unique(ivolData$Texp)
  nSlices <- length(expDates)
  if (!is.null(slices)) {
    nSlices <- length(slices)
  }
  else {
    slices <- 1:nSlices
  }
  colnum <- sqrt(nSlices * 2)
  rows <- round(colnum/2, 0)
  columns <- round(colnum, 0)
  while (rows * columns < nSlices) {
    rows <- rows + 1
  }
  atmVol <- numeric(nSlices)
  atmSkew <- numeric(nSlices)
  atmCurv <- numeric(nSlices)
  par(mfrow = c(rows, columns), mex = 0.5)
  for (slice in slices) {
    t <- expDates[slice]
    texp <- ivolData$Texp
    bidVol <- bidVols[texp == t]
    askVol <- askVols[texp == t]
    midVol <- (bidVol + askVol)/2
    if(modelVol){modVol <- modelVols[texp == t]}
    f <- (ivolData$Fwd[texp == t])[1]
    k <- log(ivolData$Strike[texp == t]/f)
    include <- !is.na(bidVol)
    kmin <- min(k[include])
    kmax <- max(k[include])
    ybottom <- 0.8 * min(bidVol[include])
    ytop <- 1.2 * max(askVol[include], na.rm = T)
    xrange <- c(kmin, kmax)
    yrange <- c(ybottom, ytop)
    plot(k, bidVol, col = "red", pch = 24, cex = 0.5, xlim = xrange, 
         ylim = yrange, main = paste("T =", format(t, digits = 2, 
                                                   nsmall = 2)), xlab = "Log-Strike", ylab = "Implied Vol.")
    points(k, askVol, col = "blue", pch = 25, cex = 0.5)
    if (modelVol) { lines(k, modVol, col = "green4",lwd=2) }
    if((!is.null(sviMatrix))){
      vol <- function(k){sqrt(svi(sviMatrix[slice,],k)/t)}
      curve(vol(x),from=kmin,to=kmax,col="orange",lwd=2,add=T);
    }
    kIn <- k[!is.na(midVol)]
    volIn <- midVol[!is.na(midVol)]
    volInterp <- function(xout) {
      stinterp(x = kIn, y = volIn, xout)$y
    }
    atmVol[slice] <- volInterp(0)
    sig <- atmVol[slice] * sqrt(t)
    atmSkew[slice] <- (volInterp(sig/10) - volInterp(-sig/10))/(2 * 
                                                                  sig/10)
    atmCurv[slice] <- (volInterp(sig/10) + volInterp(-sig/10) - 
                         2 * atmVol[slice])/(2 * (sig/10)^2)
  }
  par(mfrow = c(1, 1), mex = 1)
  par(new = F)
  return(list(expiries = expDates, atmVol = atmVol, atmSkew = atmSkew, 
              atmCurv = atmCurv))
}


