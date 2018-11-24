# SVI quartic root finder

Svi <- function(svi.params, k){
# Computes the SVI total variance at the log-strike k.
  #
  # Args:
  #  svi.params: A list of SVI parameters labeled as a, b, sig, rho, m.
  #  k: The log-strike price.
  #
  # Returns:
  #   The SVI total variance.
  
  a <- svi.params$a
  b <- svi.params$b 
  sig <- svi.params$sig 
  rho <- svi.params$rho 
  m <- svi.params$m 
  return(a + b * (rho * (k - m) + sqrt((k - m) * (k - m) + sig * sig))) 
}


ComputeSviRoots <- function(svi.data){
  # Finds roots and computes crossedness of two SVI slices.
  #
  # Args:
  #  svi.data: A data.frame with two SVI total variance slices in the rows, the earlier 
  #  slice in row 1, the later slice in row 2.  The columns are {a,b,sigma,rho,m}.
  #
  # Returns:
  #   A list with the positions of the points where the lines cross and the crossedness of the two slices.
  
  a1 <- svi.data$a[1] 
  b1 <- svi.data$b[1] 
  s1 <- svi.data$sig[1] 
  r1 <- svi.data$rho[1] 
  m1 <- svi.data$m[1] 

  a2 <- svi.data$a[2] 
  b2 <- svi.data$b[2] 
  r2 <- svi.data$rho[2] 
  m2 <- svi.data$m[2] 
  s2 <- svi.data$sig[2] 

  # The standard form of the quartic is q4 x^4 + q3 x^3 +q2 x^2 +q1 x + q0 == 0
  # Note that multiplying all the coefficients qi by a factor should have no effect on the roots.

  q2 <- 1000000 * -2 * (-3 * b1 ^ 4 * m1 ^ 2 + b1 ^ 2 * b2 ^ 2 * m1 ^ 2 + 4 * b1 ^ 2 * b2 ^ 2 * m1 * m2 + 
            b1 ^ 2 * b2 ^ 2 * m2 ^ 2 - 3 * b2 ^ 4 * m2 ^ 2 + 6 * b1 ^ 4 * m1 ^ 2 * r1 ^ 2 + 
            b1 ^ 2 * b2 ^ 2 * m1 ^ 2 * r1 ^ 2 + 4 * b1 ^ 2 * b2 ^ 2 * m1 * m2 * r1 ^ 2 + 
            b1 ^ 2 * b2 ^ 2 * m2 ^ 2 * r1 ^ 2 - 3 * b1 ^ 4 * m1 ^ 2 * r1 ^ 4 - 6 * b1 ^ 3 * b2 * m1 ^ 2 * r1 * r2 - 
            6 * b1 ^ 3 * b2 * m1 * m2 * r1 * r2 - 6 * b1 * b2 ^ 3 * m1 * m2 * r1 * r2 - 
            6 * b1 * b2 ^ 3 * m2 ^ 2 * r1 * r2 + 6 * b1 ^ 3 * b2 * m1 ^ 2 * r1 ^ 3 * r2 + 
            6 * b1 ^ 3 * b2 * m1 * m2 * r1 ^ 3 * r2 + b1 ^ 2 * b2 ^ 2 * m1 ^ 2 * r2 ^ 2 + 
            4 * b1 ^ 2 * b2 ^ 2 * m1 * m2 * r2 ^ 2 + b1 ^ 2 * b2 ^ 2 * m2 ^ 2 * r2 ^ 2 + 6 * b2 ^ 4 * m2 ^ 2 * r2 ^ 2 - 
            3 * b1 ^ 2 * b2 ^ 2 * m1 ^ 2 * r1 ^ 2 * r2 ^ 2 - 12 * b1 ^ 2 * b2 ^ 2 * m1 * m2 * r1 ^ 2 * r2 ^ 2 - 
            3 * b1 ^ 2 * b2 ^ 2 * m2 ^ 2 * r1 ^ 2 * r2 ^ 2 + 6 * b1 * b2 ^ 3 * m1 * m2 * r1 * r2 ^ 3 + 
            6 * b1 * b2 ^ 3 * m2 ^ 2 * r1 * r2 ^ 3 - 3 * b2 ^ 4 * m2 ^ 2 * r2 ^ 4 - 
            a1 ^ 2 * (b1 ^ 2 * (-1 + 3 * r1 ^ 2) - 6 * b1 * b2 * r1 * r2 + b2 ^ 2 * (-1 + 3 * r2 ^ 2)) - 
            a2 ^ 2 * (b1 ^ 2 * (-1 + 3 * r1 ^ 2) - 6 * b1 * b2 * r1 * r2 + b2 ^ 2 * (-1 + 3 * r2 ^ 2)) - 
            2 * a2 * (3 * b1 ^ 3 * m1 * r1 * (-1 + r1 ^ 2) - b1 ^ 2 * b2 * (2 * m1 + m2) * (-1 + 
                3 * r1 ^ 2) * r2 - 3 * b2 ^ 3 * m2 * r2 * (-1 + r2 ^ 2) + b1 * b2 ^ 2 * (m1 + 2 * m2) * 
               r1 * (-1 + 3 * r2 ^ 2)) + 2 * a1 * (3 * b1 ^ 3 * m1 * r1 * (-1 + r1 ^ 2) - 
              b1 ^ 2 * b2 * (2 * m1 + m2) * (-1 + 3 * r1 ^ 2) * r2 - 3 * b2 ^ 3 * m2 * r2 * (-1 + 
                r2 ^ 2) + b1 * b2 ^ 2 * (m1 + 2 * m2) * r1 * (-1 + 3 * r2 ^ 2) + 
              a2 * (b1 ^ 2 * (-1 + 3 * r1 ^ 2) - 6 * b1 * b2 * r1 * r2 + b2 ^ 2 * (-1 + 3 * r2 ^ 2))) - 
            b1 ^ 4 * s1 ^ 2 + b1 ^ 2 * b2 ^ 2 * s1 ^ 2 + b1 ^ 4 * r1 ^ 2 * s1 ^ 2 - 2 * b1 ^ 3 * b2 * r1 * r2 * 
             s1 ^ 2 + b1 ^ 2 * b2 ^ 2 * r2 ^ 2 * s1 ^ 2 + b1 ^ 2 * b2 ^ 2 * s2 ^ 2 - b2 ^ 4 * s2 ^ 2 + 
            b1 ^ 2 * b2 ^ 2 * r1 ^ 2 * s2 ^ 2 - 2 * b1 * b2 ^ 3 * r1 * r2 * s2 ^ 2 + b2 ^ 4 * r2 ^ 2 * s2 ^ 2)

  q4 <- 1000000 * (b1 ^ 4 * (-1 + r1 ^ 2) ^ 2 - 4 * b1 ^ 3 * b2 * r1 * (-1 + r1 ^ 2) * r2 - 
            4 * b1 * b2 ^ 3 * r1 * r2 * (-1 + r2 ^ 2) + b2 ^ 4 * (-1 + r2 ^ 2) ^ 2 + 
            2 * b1 ^ 2 * b2 ^ 2 * (-1 - r2 ^ 2 + r1 ^ 2 * (-1 + 3 * r2 ^ 2)))

  q0 <- 1000000 * (a1 ^ 4 + a2 ^ 4 + b1 ^ 4 * m1 ^ 4 - 2 * b1 ^ 2 * b2 ^ 2 * m1 ^ 2 * m2 ^ 2 + b2 ^ 4 * m2 ^ 4 - 
            2 * b1 ^ 4 * m1 ^ 4 * r1 ^ 2 - 2 * b1 ^ 2 * b2 ^ 2 * m1 ^ 2 * m2 ^ 2 * r1 ^ 2 + b1 ^ 4 * m1 ^ 4 * r1 ^ 4 + 
            4 * b1 ^ 3 * b2 * m1 ^ 3 * m2 * r1 * r2 + 4 * b1 * b2 ^ 3 * m1 * m2 ^ 3 * r1 * r2 - 
            4 * b1 ^ 3 * b2 * m1 ^ 3 * m2 * r1 ^ 3 * r2 - 2 * b1 ^ 2 * b2 ^ 2 * m1 ^ 2 * m2 ^ 2 * r2 ^ 2 - 
            2 * b2 ^ 4 * m2 ^ 4 * r2 ^ 2 + 6 * b1 ^ 2 * b2 ^ 2 * m1 ^ 2 * m2 ^ 2 * r1 ^ 2 * r2 ^ 2 - 
            4 * b1 * b2 ^ 3 * m1 * m2 ^ 3 * r1 * r2 ^ 3 + b2 ^ 4 * m2 ^ 4 * r2 ^ 4 + 
            4 * a2 ^ 3 * (b1 * m1 * r1 - b2 * m2 * r2) - 4 * a1 ^ 3 * (a2 + b1 * m1 * r1 - 
              b2 * m2 * r2) + 2 * b1 ^ 4 * m1 ^ 2 * s1 ^ 2 - 2 * b1 ^ 2 * b2 ^ 2 * m2 ^ 2 * s1 ^ 2 - 
            2 * b1 ^ 4 * m1 ^ 2 * r1 ^ 2 * s1 ^ 2 + 4 * b1 ^ 3 * b2 * m1 * m2 * r1 * r2 * s1 ^ 2 - 
            2 * b1 ^ 2 * b2 ^ 2 * m2 ^ 2 * r2 ^ 2 * s1 ^ 2 + b1 ^ 4 * s1 ^ 4 - 2 * b1 ^ 2 * b2 ^ 2 * m1 ^ 2 * s2 ^ 2 + 
            2 * b2 ^ 4 * m2 ^ 2 * s2 ^ 2 - 2 * b1 ^ 2 * b2 ^ 2 * m1 ^ 2 * r1 ^ 2 * s2 ^ 2 + 
            4 * b1 * b2 ^ 3 * m1 * m2 * r1 * r2 * s2 ^ 2 - 2 * b2 ^ 4 * m2 ^ 2 * r2 ^ 2 * s2 ^ 2 - 
            2 * b1 ^ 2 * b2 ^ 2 * s1 ^ 2 * s2 ^ 2 + b2 ^ 4 * s2 ^ 4 + 4 * a2 * (b1 * m1 * r1 - b2 * m2 * r2) * 
             (-2 * b1 * b2 * m1 * m2 * r1 * r2 + b1 ^ 2 * (m1 ^ 2 * (-1 + r1 ^ 2) - s1 ^ 2) + 
              b2 ^ 2 * (m2 ^ 2 * (-1 + r2 ^ 2) - s2 ^ 2)) - 4 * a1 * (a2 + b1 * m1 * r1 - 
              b2 * m2 * r2) * (a2 ^ 2 - 2 * b1 * b2 * m1 * m2 * r1 * r2 + 2 * a2 * (b1 * m1 * r1 - 
                b2 * m2 * r2) + b1 ^ 2 * (m1 ^ 2 * (-1 + r1 ^ 2) - s1 ^ 2) + 
              b2 ^ 2 * (m2 ^ 2 * (-1 + r2 ^ 2) - s2 ^ 2)) + 2 * a2 ^ 2 * 
             (-6 * b1 * b2 * m1 * m2 * r1 * r2 + b1 ^ 2 * (m1 ^ 2 * (-1 + 3 * r1 ^ 2) - s1 ^ 2) + 
              b2 ^ 2 * (m2 ^ 2 * (-1 + 3 * r2 ^ 2) - s2 ^ 2)) + 
            2 * a1 ^ 2 * (3 * a2 ^ 2 - 6 * b1 * b2 * m1 * m2 * r1 * r2 + 6 * a2 * (b1 * m1 * r1 - 
                b2 * m2 * r2) + b1 ^ 2 * (m1 ^ 2 * (-1 + 3 * r1 ^ 2) - s1 ^ 2) + 
              b2 ^ 2 * (m2 ^ 2 * (-1 + 3 * r2 ^ 2) - s2 ^ 2)))

  q3 <- 1000000 * -4 * (b1 ^ 4 * m1 * (-1 + r1 ^ 2) ^ 2 - b1 ^ 3 * r1 * (-1 + r1 ^ 2) * 
             (a1 - a2 + b2 * (3 * m1 + m2) * r2) + b2 ^ 3 * (-1 + r2 ^ 2) * 
             ((a1 - a2) * r2 + b2 * m2 * (-1 + r2 ^ 2)) + b1 * b2 ^ 2 * r1 * 
             (a1 - 3 * a1 * r2 ^ 2 - b2 * (m1 + 3 * m2) * r2 * (-1 + r2 ^ 2) + 
              a2 * (-1 + 3 * r2 ^ 2)) + b1 ^ 2 * b2 * ((a1 - a2) * (-1 + 3 * r1 ^ 2) * r2 + 
              b2 * (m1 + m2) * (-1 - r2 ^ 2 + r1 ^ 2 * (-1 + 3 * r2 ^ 2))))

  q1 <- 1000000 * 4 * (-(b1 ^ 4 * m1 ^ 3) + b1 ^ 2 * b2 ^ 2 * m1 ^ 2 * m2 + b1 ^ 2 * b2 ^ 2 * m1 * m2 ^ 2 - 
            b2 ^ 4 * m2 ^ 3 + 2 * b1 ^ 4 * m1 ^ 3 * r1 ^ 2 + b1 ^ 2 * b2 ^ 2 * m1 ^ 2 * m2 * r1 ^ 2 + 
            b1 ^ 2 * b2 ^ 2 * m1 * m2 ^ 2 * r1 ^ 2 - b1 ^ 4 * m1 ^ 3 * r1 ^ 4 - b1 ^ 3 * b2 * m1 ^ 3 * r1 * r2 - 
            3 * b1 ^ 3 * b2 * m1 ^ 2 * m2 * r1 * r2 - 3 * b1 * b2 ^ 3 * m1 * m2 ^ 2 * r1 * r2 - b1 * b2 ^ 3 * m2 ^ 3 * r1 * r2 + b1 ^ 3 * b2 * m1 ^ 3 * r1 ^ 3 * r2 + 3 * b1 ^ 3 * b2 * m1 ^ 2 * m2 * 
             r1 ^ 3 * r2 + b1 ^ 2 * b2 ^ 2 * m1 ^ 2 * m2 * r2 ^ 2 + b1 ^ 2 * b2 ^ 2 * m1 * m2 ^ 2 * r2 ^ 2 + 
            2 * b2 ^ 4 * m2 ^ 3 * r2 ^ 2 - 3 * b1 ^ 2 * b2 ^ 2 * m1 ^ 2 * m2 * r1 ^ 2 * r2 ^ 2 - 
            3 * b1 ^ 2 * b2 ^ 2 * m1 * m2 ^ 2 * r1 ^ 2 * r2 ^ 2 + 3 * b1 * b2 ^ 3 * m1 * m2 ^ 2 * r1 * r2 ^ 3 + 
            b1 * b2 ^ 3 * m2 ^ 3 * r1 * r2 ^ 3 - b2 ^ 4 * m2 ^ 3 * r2 ^ 4 + a1 ^ 3 * (b1 * r1 - b2 * r2) + 
            a2 ^ 3 * (-(b1 * r1) + b2 * r2) + a2 ^ 2 * (b1 ^ 2 * (m1 - 3 * m1 * r1 ^ 2) + 3 * b1 * b2 * (m1 + m2) * r1 * r2 + b2 ^ 2 * m2 * (1 - 3 * r2 ^ 2)) + 
            a1 ^ 2 * (b1 ^ 2 * (m1 - 3 * m1 * r1 ^ 2) + 3 * b1 * r1 * (-a2 + b2 * (m1 + m2) * r2) + 
              b2 * (3 * a2 * r2 + b2 * (m2 - 3 * m2 * r2 ^ 2))) - b1 ^ 4 * m1 * s1 ^ 2 + 
            b1 ^ 2 * b2 ^ 2 * m2 * s1 ^ 2 + b1 ^ 4 * m1 * r1 ^ 2 * s1 ^ 2 - b1 ^ 3 * b2 * m1 * r1 * r2 * s1 ^ 2 - 
            b1 ^ 3 * b2 * m2 * r1 * r2 * s1 ^ 2 + b1 ^ 2 * b2 ^ 2 * m2 * r2 ^ 2 * s1 ^ 2 + 
            b1 ^ 2 * b2 ^ 2 * m1 * s2 ^ 2 - b2 ^ 4 * m2 * s2 ^ 2 + b1 ^ 2 * b2 ^ 2 * m1 * r1 ^ 2 * s2 ^ 2 - 
            b1 * b2 ^ 3 * m1 * r1 * r2 * s2 ^ 2 - b1 * b2 ^ 3 * m2 * r1 * r2 * s2 ^ 2 + 
            b2 ^ 4 * m2 * r2 ^ 2 * s2 ^ 2 + a2 * (b1 ^ 2 * b2 * r2 * (m1 ^ 2 * (-1 + 3 * r1 ^ 2) + 
                2 * m1 * m2 * (-1 + 3 * r1 ^ 2) - s1 ^ 2) + b1 ^ 3 * r1 * (-3 * m1 ^ 2 * 
                 (-1 + r1 ^ 2) + s1 ^ 2) + b2 ^ 3 * r2 * (3 * m2 ^ 2 * (-1 + r2 ^ 2) - s2 ^ 2) + 
              b1 * b2 ^ 2 * r1 * (m1 * m2 * (2 - 6 * r2 ^ 2) + m2 ^ 2 * (1 - 3 * r2 ^ 2) + s2 ^ 2)) + 
            a1 * (3 * a2 ^ 2 * (b1 * r1 - b2 * r2) + a2 * (2 * b1 ^ 2 * m1 * (-1 + 3 * r1 ^ 2) - 
                6 * b1 * b2 * (m1 + m2) * r1 * r2 + 2 * b2 ^ 2 * m2 * (-1 + 3 * r2 ^ 2)) + 
              b1 ^ 3 * r1 * (3 * m1 ^ 2 * (-1 + r1 ^ 2) - s1 ^ 2) + b1 ^ 2 * b2 * r2 * ( 
                m1 * m2 * (2 - 6 * r1 ^ 2) + m1 ^ 2 * (1 - 3 * r1 ^ 2) + s1 ^ 2) + 
              b1 * b2 ^ 2 * r1 * (2 * m1 * m2 * (-1 + 3 * r2 ^ 2) + m2 ^ 2 * (-1 + 3 * r2 ^ 2) - 
                s2 ^ 2) + b2 ^ 3 * r2 * (-3 * m2 ^ 2 * (-1 + r2 ^ 2) + s2 ^ 2)))

  term16 <- (2 * q2 ^ 3 + 27 * q3 ^ 2 * q0 - 72 * q4 * q2 * q0 - 9 * q3 * q2 * q1 + 27 * q4 * q1 ^ 2)

  term21 <- (q2 ^ 2 / 4 + 3 * q4 * q0 - 3 * q3 * q1 / 4)

  term1sq <- -256 * term21 ^ 3 + term16 ^ 2

  term1 <- sqrt(term1sq+0*1i)  # Note use of complex arithmetic in R

  term23 <- (term16 + term1) ^ (1/3)

  term22 <- 3 * q4 * term23 

  temp1 <- (4 * 2 ^ (1 / 3) * term21)
  temp2 <- (3 * 2 ^ (1 / 3) * q4)
  temp3 <- q3 ^ 2 / (4 * q4 ^ 2) - (2 * q2) / (3 * q4)
  temp4 <- temp1 / term22 + term23 / temp2 

  rr <- sqrt(temp3 + temp4) 

  temp5 <- q3 ^ 2 / (2 * q4 ^ 2) - (4 * q2) / (3 * q4)
  temp6 <- (-q3 ^ 3 / 4 + q4 * q3 * q2 - 2 * q4 ^ 2 * q1) / (q4 ^ 3)

  ee <- q3 ^ 2 / (2 * q4 ^ 2) - (4 * q2) / (3 * q4) - (4 * 2 ^ (1 / 3) * term21) / term22 - term23 / (3 * 2 ^ (1 / 3) * q4) -
          (-q3 ^ 3 / 4 + q4 * q3 * q2 - 2 * q4 ^ 2 * q1) / (q4 ^ 3 * rr) 

  dd <- q3 ^ 2 / (2 * q4 ^ 2) - (4 * q2) / (3 * q4) - (4 * 2 ^ (1 / 3) * term21) / term22 - term23 / (3 * 2 ^ (1 / 3) * q4) + 
          (-q3 ^ 3 / 4 + q4 * q3 * q2 - 2 * q4 ^ 2 * q1) / (q4 ^ 3 * rr) 

  temp7 <- -q3 / (4 * q4) 

  # Potential roots are given by
  roots <- c(
             -q3 / (4 * q4) +rr / 2 + sqrt(dd) / 2,
             -q3 / (4 * q4) +rr / 2 - sqrt(dd) / 2,
             -q3 / (4 * q4) -rr / 2 + sqrt(ee) / 2,
             -q3 / (4 * q4) -rr / 2 - sqrt(ee) / 2) 

  # Need to check these are really roots
  kr <- roots * (abs(Im(roots)) < 10 ^ (-10)) 
  test <- function(k){
                      (a1 + b1 * (r1 * (k - m1) + sqrt((k - m1) ^ 2 + s1 ^ 2))) - 
                      (a2 + b2 * (r2 * (k - m2) + sqrt((k - m2) ^ 2 + s2 ^ 2)))
  } 
            
  num <- which(abs(test(kr)) < 10 ^ (-10))  # Which potential root is actually a root?

  roots <- sort(Re(kr[num]))  # Roots in ascending order
  n.roots <- length(roots)

  crossedness <- 0

  if (n.roots>1) { 
    mid.points <- (roots[1:(n.roots-1)]+roots[2:n.roots])/2 
  } 
  else {
    mid.points <- c()
  }  

  if (n.roots>0) {
    sample.points <- c(roots[1] - 1, mid.points, roots[n.roots] + 1)  # Choose some sensible sampling points
	  svi.short <- Svi(svi.data[1, ], sample.points) 
	  svi.long <- Svi(svi.data[2, ], sample.points) 
	  crossedness <- max(c(0, max(svi.short - svi.long)))  # Maximal amount of crossing bounded below by zero
  } 

  return(list(roots=roots, crossedness=crossedness))   
}

# Example (Remove comments to run)
# slice1 <- c(1.8,.8,0,-.5,0)
# slice2 <- c(1,1,1,-.5,0)
# test.data <- as.data.frame(rbind(slice1, slice2))
# colnames(test.data)<- c("a","b","sig","rho","m")
# ComputeSviRoots(test.data)

