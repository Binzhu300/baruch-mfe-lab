# Computes var swaps for given vector of expirations
# This version has some error trapping that returns NA for the variance swap
# when the integration fails

###########################################
# Variance swap
###########################################
sviVarSwap <- function(sviMatrix,texp){
  
  nSlices <- length(texp)
  varRes <- numeric(nSlices)
  
  for (slice in 1:nSlices){
    
    t <- texp[slice]
    volBS <- function(k){sqrt(svi(sviMatrix[slice,],k)/t)}
    # cTilde <- function(k){BSFormula(1, exp(k), t, r=0, volBS(k))*exp(-k)}
    # pTilde <- function(k){BSFormulaPut(1, exp(k), t, r=0, volBS(k))*exp(-k)}
    cTilde <- function(y){BSFormula(1, 1/y, t, r=0, volBS(-log(y)))}
    pTilde <- function(y){BSFormulaPut(1, y, t, r=0, volBS(log(y)))/y^2}
    possibleError <- tryCatch(callIntegral <- integrate(cTilde,lower=0,upper=1)$value,error=function(e) e)
    possibleError <- tryCatch(putIntegral <- integrate(pTilde,lower=0,upper=1)$value,error=function(e) e)
    if(!inherits(possibleError,"error")){
      varRes[slice] <- 2*(callIntegral+putIntegral)/t
    } else {varRes[slice] <- NA}
  }
  return(varRes)  # Returns vector of variance swaps in variance (vol. squared) terms
}

###########################################
# Total variance E[w]
###########################################
svi.w <- function (sviMatrix) 
{
  nSlices <- dim(sviMatrix)[1]
  w.res <- numeric(nSlices)
  
  for (slice in 1:nSlices) {
    sigBS <- function(k) { sqrt(svi(sviMatrix[slice, ], k)) }
    
    call <- function(K){
      k <- log(K)
      BSFormula(1, K, 1, r = 0, sigBS(k))/K^2
    }
    
    put <- function(K){
      k <- log(K)
      BSFormulaPut(1, K, 1, r = 0, sigBS(k))/K^2
    }
    
    possibleError <- tryCatch(callIntegral <- integrate(call, 
                                                        lower = 1, upper = Inf)$value, error = function(e) e)
    possibleError <- tryCatch(putIntegral <- integrate(put, 
                                                       lower = 0, upper = 1)$value, error = function(e) e)
    
    if (!inherits(possibleError, "error")) {
      w.res[slice] <- 2 * (callIntegral + putIntegral)
    }
    else {
      w.res[slice] <- NA
    }
  }
  return(w.res)
}

###########################################
# Forward variance \xi_t(T)
###########################################
svi.xi <- function (sviMatrix,expiries) 
{
  w <- svi.w(sviMatrix)
  xi.vec <- c(w[1]/expiries[1],diff(w)/diff(expiries))
  xi.curve.raw <- function(tau){
    taup <- min(tau,max(expiries)) # Extrapolation at constant forward variance
    xi.vec[sum(expiries < taup)+1]
  }
  
  xi.curve <- function(tau){sapply(tau,xi.curve.raw)}
  return(xi.curve) # Returns a function of tau
}

###########################################
# Gamma swap 
###########################################
sviGammaSwap <- function(sviMatrix,texp){
  
  nSlices <- length(texp)
  varRes <- numeric(nSlices)
  
  for (slice in 1:nSlices){
    
    t <- texp[slice]
    volBS <- function(k){sqrt(svi(sviMatrix[slice,],k)/t)}
    cTilde <- function(y){BSFormula(1, 1/y, t, r=0, volBS(-log(y)))/y}
    pTilde <- function(y){BSFormulaPut(1, y, t, r=0, volBS(log(y)))/y}
    possibleError <- tryCatch(callIntegral <- integrate(cTilde,lower=0,upper=1)$value,error=function(e) e)
    possibleError <- tryCatch(putIntegral <- integrate(pTilde,lower=0,upper=1)$value,error=function(e) e)
    if(!inherits(possibleError,"error")){
      varRes[slice] <- 2*(callIntegral+putIntegral)/t
    } else {varRes[slice] <- NA}
  }
  return(varRes)  # Returns vector of variance swaps in variance (vol. squared) terms
}

###########################################
#  E[X_T^2]
###########################################
svi.EX2 <- function (sviMatrix) 
{
  nSlices <- dim(sviMatrix)[1]
  eX2.res <- numeric(nSlices)
  
  for (slice in 1:nSlices) {
    sigBS <- function(k) { sqrt(svi(sviMatrix[slice, ], k)) }
    
    call <- function(K){
      k <- log(K)
      BSFormula(1, K, 1, r = 0, sigBS(k))*(1-k)/K^2
    }
    
    put <- function(K){
      k <- log(K)
      BSFormulaPut(1, K, 1, r = 0, sigBS(k))*(1-k)/K^2
    }
    
    possibleError <- tryCatch(callIntegral <- integrate(call, 
                                                        lower = 1, upper = Inf)$value, error = function(e) e)
    possibleError <- tryCatch(putIntegral <- integrate(put, 
                                                       lower = 0, upper = 1)$value, error = function(e) e)
    
    if (!inherits(possibleError, "error")) {
      eX2.res[slice] <- 2 * (callIntegral + putIntegral)
    }
    else {
      eX2.res[slice] <- NA
    }
  }
  return(eX2.res)
}

  ###########################################
  #  E[X_T^2] + 2 E[X_T]
  ###########################################
  svi.EX2EX <- function (sviMatrix)
  {
    nSlices <- dim(sviMatrix)[1]
    eX2.res <- numeric(nSlices)

    for (slice in 1:nSlices) {
      sigBS <- function(k) { sqrt(svi(sviMatrix[slice, ], k)) }

      call <- function(K){
        k <- log(K)
         - BSFormula(1, K, 1, r = 0, sigBS(k))*k/K^2
      }

      put <- function(K){
        k <- log(K)
        - BSFormulaPut(1, K, 1, r = 0, sigBS(k))*k/K^2
      }

      possibleError <- tryCatch(callIntegral <- integrate(call,
                                                          lower = 1, upper = Inf)$value, error = function(e) e)
      possibleError <- tryCatch(putIntegral <- integrate(put,
                                                         lower = 0, upper = 1)$value, error = function(e) e)

      if (!inherits(possibleError, "error")) {
        eX2.res[slice] <- 2 * (callIntegral + putIntegral)
      }
      else {
        eX2.res[slice] <- NA
      }
    }
    return(eX2.res)
}

###########################################
#  Stochasticity \zeta_t(T)
###########################################
svi.zeta2 <- function(sviMatrix){

  return(svi.EX2(sviMatrix)
         -svi.w(sviMatrix)^2/4
         -svi.w(sviMatrix))
}

####################################################
#  Stlibaraochasticity \zeta_t(T) alternative computation
####################################################
svi.zeta <- function(sviMatrix){
  return( svi.EX2EX(sviMatrix)  - svi.w(sviMatrix)^2/4) }


# vols <- sqrt(sviVarSwap(sviMatrix,texp))
# plot(texp,vols,col="red",type="b")
# volsG <- sqrt(sviGammaSwap(sviMatrix,texp))
# points(texp,volsG,col="blue",type="b")
# plot(texp,vols-volsG,col="dark green",type="b")
# sZ <- sviZSwap(sviMatrix,texp)/texp
# plot(texp,sZ)
# varx <- sZ - (sviVarSwap(sviMatrix,texp)/2)^2
# plot(texp,varx)
# # Now subtract the variance swap 
# varx.qv <- varx - sviVarSwap(sviMatrix,texp)
# plot(texp,varx.qv,col="red",pch=20)
