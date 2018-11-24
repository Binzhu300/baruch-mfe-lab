library(stinepack)

# Function to compute total variance given k and t
# texp is a vector of times to expiration
sviW <- function(sviMatrix,texp,k,t){

# Vector of SVI variance for a given strike k
sviWk <- function(k){
m <- dim(sviMatrix)[1];
sapply(1:m,function(i){svi(sviMatrix[i,],k)}); 
}

# Function to compute interpolated variance for any strike and expiration
wInterp <- function(k,t){stinterp(texp,sviWk(k),t)$y};

# Vectorized function that returns implied total variance for vectors k and t
return(sapply(k,function(k1){wInterp(k1,t)}));

} # End of sviW
