source("BlackScholes.R")

# SVI total variance
svi <- function(sviparams,k){
  a <- sviparams$a;
  b <- sviparams$b;
  sig <- sviparams$sig;
  rho <- sviparams$rho;
  m <- sviparams$m;
  return(a + b *(rho*(k-m)+ sqrt((k-m)*(k-m) + sig*sig)));
}

# SVI total variance skew
svi.skew <- function(sviparams,k){
  a <- sviparams$a;
  b <- sviparams$b;
  sig <- sviparams$sig;
  rho <- sviparams$rho;
  m <- sviparams$m;
  
  skew <- b * ((k-m)/sqrt((k-m)^2+sig^2) + rho )
  return(skew)
}


# SVI Density
psvi <- function(sviparams,k){sapply(k,function(k){
	
	b <- sviparams$b;
	sig <- sviparams$sig;
	rho <- sviparams$rho;
	m <- sviparams$m;
	dsqrt <- sqrt((k-m)^2+sig^2);
	w <- svi(sviparams,k);
	wk <- b*rho+b*(k-m)/dsqrt;
	wkk <- b*sig^2/dsqrt^3 ;
	tmp2 <- (1-k/(2*w)*wk)^2-1/4*(1/4+1/w)*wk^2+1/2*wkk;
	tmp1 <- exp(-(2*k - w)^2/(8*w))/(2*sqrt(2*pi)*sqrt(w));
	return(2*tmp1*tmp2);
	
	});
	
	}
	
# Params from TVS page 	
# sviparams <- list(a=0.04,b=0,rho=0,sig=0.01,m=0);
# curve(psvi(sviparams,x),from=-1,to=1);
# curve(dnorm(x,mean=0.02,sd=.2),col="red",lty=2,add=T);

