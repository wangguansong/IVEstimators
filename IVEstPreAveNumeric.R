IVEstPreAveNumeric <- function(steps,
                               FUN=function(x){pmin(x, 1-x)}){
  # Numerically compute the parameters (integrals of the weight 
  # function) for the Pre-Averaging approach.
  #
  # Reference:
  # Jacod, Li, Mykland, Podolskij, and Vetter (2009), Microstructure
  # Noise in the Continuous Case: The Pre-averaging Approach, Stochastic
  # Processes and Their Applications
  #
  # Args:
  #   steps: Bandwidth used.
  #   FUN: Weighting function.
  # Returns: A vector of parameter values, in order of
  #   (psi1, psi2, Phi11, Phi12, Phi22).
  ##################################################
  # Input Check
  if (steps==0) {
    return(c(1, 1/12, 1/6, 1/96, 151/80640))
  }
  grids <- 0:steps/steps
  FUNvalue <- FUN(grids)
  FUN1stdev <- c(FUNvalue[2]-FUNvalue[1],
                 (FUNvalue[3:(steps+1)]-FUNvalue[1:(steps-1)])/2,
                 FUNvalue[steps+1]-FUNvalue[steps]) * steps
  psi1 <- sum((FUN1stdev[1:steps]^2 + FUN1stdev[2:(steps+1)]^2) / 2) /
          steps
  psi2 <- sum((FUNvalue[1:steps]^2 + FUNvalue[2:(steps+1)]^2) / 2) /
          steps
  phi1 <- numeric(steps+1) 
  phi2 <- numeric(steps+1) 
  for (j in 0:steps) {
    phi1[j+1] <-
      sum(FUN1stdev[j:steps+1]*FUN1stdev[j:steps+1-j]) / steps
    phi2[j+1] <- sum(FUNvalue[j:steps+1]*FUNvalue[j:steps+1-j]) / steps
  }
  Phi11 <- sum((phi1[1:steps]^2 + phi1[2:(steps+1)]^2) / 2) / steps
  Phi22 <- sum((phi2[1:steps]^2 + phi2[2:(steps+1)]^2) / 2) / steps
  Phi12 <- sum(((phi1*phi2)[1:steps] + (phi1*phi2)[2:(steps+1)]) / 2)/
           steps
  return(c(psi1, psi2, Phi11, Phi12, Phi22))
}
