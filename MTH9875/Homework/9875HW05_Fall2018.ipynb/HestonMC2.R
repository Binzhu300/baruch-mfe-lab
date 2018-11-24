# Heston simulation code from Lecture 4

source("BlackScholes.R");

#-----------------------------------------------------------------------------
# A little function needed later

is.even <- function(j){as.logical((j+1) %% 2)};

#-----------------------------------------------------------------------------
# Alfonsi (2010) timestep: Full truncation

evolveAlfonsiF <- function(v,x,dt,W1,W2,L){
        
        eldt2 <- exp(-lambda*dt/2);
        
        #Variance process 
        vbarp <- vbar - eta^2/(4*lambda);
        psi <- (1-eldt2)/lambda;
        v1 <- v*eldt2+lambda*vbarp*psi;
        v2 <- (v1 > 0) * v1; # Take v2 = 0 if v1<0, else v2=v1
        par <- sqrt(v2) + eta/2 * sqrt(dt)*W2;
        vf <- eldt2*par*par +lambda*vbarp*psi + v1 - v2;
                                # Full truncation
        
        # Log-stock process (Andersen equation (33))
        vvf <- (v+vf > 0) * (v+vf);
        dw <- vvf/2*dt;
        x <- x  - dw/2 + rho2m1*sqrt(dw)*W1 +
        	rho/eta*(lambda*dw + vf-v -lambda*vbar*dt) ;
        # Impose martingale constraint
        x <- x - log(mean(exp(x)));
        v <- vf;
        return(cbind(x,v,dw));
}

#-----------------------------------------------------------------------------
# Heston Monte Carlo code from Lecture 4

HestonMC2 <- function(params){
  
  is.even <- function(j){as.logical((j+1) %% 2)}
  
  res <- function(S0, T, AK, N, m, evolve,exactVols=NULL)
  {
    
    lambda <<- params$lambda;
    rho <<- params$rho;
    eta <<- params$eta;
    vbar <<- params$vbar;
    v0 <<- params$v;
    
    n <- m*2; #n is number of timesteps = 2*m so we can use Romberg extrapolation
    sqrt2 <- sqrt(2);
    rho2m1 <<- sqrt(1-rho*rho);
    vbarp <<- vbar - eta^2/(4*lambda);
    
    negCount <- 0;
    
    # We use a vertical array, one element per M.C. path
    x <- rep(0,N); v <- rep(1,N)*v0;
    xm <- x; vm <- v;
    W1m <- rep(0,N); W2m <- rep(0,N); 
    
    # Loop for bias computation (N small, n big)
    for (i in 1:n)
    {
      # Two sets of correlated normal random vars.
      
      W1 <- rnorm(N);
      W2 <- rnorm(N);
      W1 <- W1 - mean(W1); W1 <- W1/sd(W1);
      W2 <- W2 - mean(W2); W2 <- W2/sd(W2);
      # Now W1 and W2 are forced to have mean=0 and sd=1
      
      W2p <- W2 - cor(W1,W2)*W1; # Eliminate actual correlation
      W2p <- W2p - mean(W2p); W2 <- W2p/sd(W2p); 
      # Now W1 and W2 have mean=0, sd=1 and correlation=0
      
      L <- rbinom(N, size=1, prob=1/2); # Bernoulli rv for NV step
      
      # Add code for subgrid
      W1m <- W1m + W1/sqrt2; W2m <- W2m + W2/sqrt2; # N(0,1) rv's for subgrid
      
      if (is.even(i)) {
        #print(c(i,mean(W1m),mean(W2m),sd(W1m),sd(W2m),cor(W1m,W2m)));
        resm <- evolve(vm,xm,T/m,W1m,W2m,L);
        xm <- resm[,1];
        vm <- resm[,2];
        W1m <- rep(0,N); W2m <- rep(0,N);
      }
      
      res <- evolve(v,x,T/n,W1,W2,L);
      x <- res[,1];
      v <- res[,2];
      negCount <- negCount +mean(v<0)/n; #Probability of negative variance per path per timestep
      
    }
    
    S <- S0*exp(x);
    Sm <- S0*exp(xm);
    
    # Now we have three vectors of final stock prices
    
    M <- length(AK);
    AV <- numeric(M); AVdev <- numeric(M);
    BSV <- numeric(M); BSVH <- numeric(M); BSVL <- numeric(M);
    iv2SD <- numeric(M); bias <- numeric(M);
    AVm <- numeric(M); AVmdev <- numeric(M);
    BSVm <- numeric(M); BSVHm <- numeric(M); BSVLm <- numeric(M);
    iv2SDm <- numeric(M);
    AV1 <- numeric(M); AV1dev <- numeric(M);
    BSV1 <- numeric(M); BSVH1 <- numeric(M); BSVL1 <- numeric(M);
    iv2SDrom <- numeric(M);biasRom <- numeric(M);
    
    # Evaluate mean call value for each path
    for (i in 1:M)
    {
      # 2*m timesteps
      K <- AK[i];
      V <- (S>K)*(S - K); # Boundary condition for European call
      AV[i] <- mean(V);
      AVdev[i] <- sqrt(var(V)/length(V));  
      
      BSV[i] <- BSImpliedVolCall(S0, K, T, 0, AV[i]);
      BSVL[i] <- BSImpliedVolCall(S0, K, T, 0, AV[i] - AVdev[i]);
      BSVH[i] <- BSImpliedVolCall(S0, K, T, 0, AV[i] + AVdev[i]);
      iv2SD[i] <- (BSVH[i]-BSVL[i]);
      
      # m timesteps
      Vm <- (Sm>K)*(Sm - K); # Boundary condition for European call
      AVm[i] <- mean(Vm);
      AVmdev[i] <- sd(Vm) / sqrt(N);
      BSVm[i] <- BSImpliedVolCall(S0, K, T, 0, AVm[i]);
      BSVLm[i] <- BSImpliedVolCall(S0, K, T, 0, AVm[i] - AVmdev[i]);
      BSVHm[i] <- BSImpliedVolCall(S0, K, T, 0, AVm[i] + AVmdev[i]);
      iv2SDm[i] <- (BSVH[i]-BSVL[i]);
      
      # Richardson extrapolation estimates 
      V1 <- 2*V - Vm;
      AV1[i] <- mean(V1);
      AV1dev[i] <- sd(V1) / sqrt(N);
      BSV1[i] <- BSImpliedVolCall(S0, K, T, 0, AV1[i]);
      BSVL1[i] <- BSImpliedVolCall(S0, K, T, 0, AV1[i] - AV1dev[i]);
      BSVH1[i] <- BSImpliedVolCall(S0, K, T, 0, AV1[i] + AV1dev[i]);
      iv2SDrom[i] <- (BSVH1[i]-BSVL1[i]);
      
      if(!is.null(exactVols)) {bias <- BSV-exactVols};
      if(!is.null(exactVols)) {biasRom <- BSV1-exactVols};
    }
    
    
    l.AK <- length(AK)      
    data.out <- data.frame(AK,rep(N,l.AK),rep(2*m,l.AK),BSV,bias,iv2SD,BSVm,BSV1,biasRom,iv2SDrom) 
    names(data.out) <- c("Strikes","Paths","Steps","ivol","bias","twoSd","ivolm", "ivolRichardson", "biasRichardson","twoSdRichardson") 
    return(data.out) 
    
    
  }
  return(res)
}

