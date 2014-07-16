
IVEstRKFromGamma <- function(Gamma, weights=numeric(), scale.factor=1) {
  # Realized Kernel estimator, computed by given gamma's
  # Args:
  #   Gamma: a vector of gamma_0, gamma_1, ...
  #   weights:
  #   scale.factor: a vector, default 1
  # Returns:
  #   estimator value

  ##################################################
  # Input Check

  H <- min(length(weights), sum(!is.na(Gamma))-1)

  ##################################################

  # Estimator
  Gamma[-1] <- 2 * Gamma[-1]
  return(Gamma[1] + ifelse(H>0, sum(Gamma[1+1:H] * weights[1:H]), 0))
}
