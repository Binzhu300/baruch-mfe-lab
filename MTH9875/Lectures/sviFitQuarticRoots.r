# Fit SVI independently to each slice subject to crossing penalty

sviFitQR <- function(ivolData,sviGuess,penaltyFactor=100){

callVals <- as.numeric(ivolData$CallMid);
askVols <- as.numeric(ivolData$Ask);
bidVols <- as.numeric(ivolData$Bid);
midVols <- (bidVols+askVols)/2;
volSpreads <- (askVols-bidVols)/2;

expDates <- unique(ivolData$Texp);
nSlices <- length(expDates);
slices <- 1:nSlices;
sviMatrix <- sviGuess;

######################################
#for (slice in rev(slices)){  # Start from the latest slice

for (slice in slices){  # Start from the earliest slice

t <- expDates[slice];
texp <- ivolData$Texp;

midVal <- callVals[texp==t];

pick <- !is.na(midVal);
midVal <- midVal[pick];
midVols1 <- midVols[texp==t][pick];
volSpreads1 <- volSpreads[texp==t][pick];
f <- (ivolData$Fwd[texp==t])[1];
k <- log(ivolData$Strike[texp==t]/f)[pick]; # Plot vs log-strike

# Compute squared distance between actual mid-price and theoretical mid-price
sqDist <- function(sviparams){

    tmp1 <- as.list(sviparams);
    names(tmp1) <- c("a","b","sig","rho","m")# tmp1 is w-style
    sviVar <- svi(tmp1,k)/t; # Convert to variance-style

    outVal <- BSFormula(f, f*exp(k), t, 0, sqrt(abs(sviVar))) # Compute Black-Scholes call values

    tmp <- sum((midVal-outVal)^2,na.rm=T);
    
    return(tmp);
};

# Compute squared distance between mid-vol and theoretical mid-vol, weighted by spread
sqDistOld <- function(sviparams){

    tmp1 <- as.list(sviparams);
    names(tmp1) <- c("a","b","sig","rho","m")# tmp1 is w-style
    sviVar <- svi(tmp1,k)/t; # Convert to variance-style
    sviVol <- sqrt(abs(sviVar));
    
    tmp <- sum((midVols1-sviVol)^2/(volSpreads1^2),na.rm=T);
    
    return(tmp);
};

# Normalized squared distance
sqDistN <- function(sviparams){sqDist(sviparams)/sqDist(sviGuess[slice,])}

crossPenalty <- function(sviparams){
    
    tmp1 <- as.list(sviparams);
    names(tmp1) <- c("a","b","sig","rho","m")# tmp1 is w-style
    
    cPenalty <- 0;
    if (slice > 1) { # Compare with previous slice
    	slicePlusPrevious <- rbind(as.data.frame(tmp1),as.data.frame(sviMatrix[slice-1,]));
    	cPenalty <- ComputeSviRoots(slicePlusPrevious)$crossedness} 
    else { # Ensure no negative variances on first slice
        minVar <- tmp1$a+tmp1$b*tmp1$sig*sqrt(abs(1-tmp1$rho^2));
        negVarPenalty <- min(100,exp(-1/minVar));
        cPenalty <- negVarPenalty;
        }  
    if (slice < nSlices) { # Compare with next slice
    	slicePlusNext <- rbind(as.data.frame(tmp1),as.data.frame(sviMatrix[slice+1,]));
    	cPenalty <- cPenalty + ComputeSviRoots(slicePlusNext)$crossedness};
    
    return(cPenalty*1000);
};

# Compute objective function
obj <- function(sviparams){sqDistN(sviparams)+crossPenalty(sviparams)};

# Optimize
fit <- optim(sviGuess[slice,],obj,method = "L-BFGS-B", lower=c(-Inf,0,0,-1,-Inf));
sviMatrix[slice,] <- fit$par;
if(abs(fit$par["rho"])>1){ #Only send to L-BFGS-B if rho is outside permitted range
fit <- optim(sviGuess[slice,],obj,method="L-BFGS-B",lower=c(-1000,0,0.00000001,-.999,-10),upper=c(+1000,100,100,+.999,+10));
};


                    }# end of for loop
######################################
sviMatrix <- as.data.frame(sviMatrix);
colnames(sviMatrix) <- c("a","b","sig","rho","m");
return(sviMatrix);
}# end of function
