IVEstLowerBound <- function(N, ma.parameter=1, sigma2=1, sigma2v=1) {
  # Parametric Lower Bound
  #
  # Arguments:
  #   N: sample size
  #   ma.parameter: beta, MA parameter of noise, with beta_0
  #   sigma2: integrated variation
  #   sigma2v: variance of white noise
  #
  # Returns:
  #   A scalar, parametric lower bound
  
  CN <- matrix(0, nrow=N, ncol=N)
  theta <- c(1, diff(c(ma.parameter, 0)))
  
  q <- length(theta) - 1
  phi <- numeric(q+1)
  
  for (i in 0:q) {
    phi[i+1] <- sum(theta[1:(q+1-i)] * theta[(1+i):(q+1)])
  }
  
  
  diag(CN) <- sigma2/N + phi[1]*sigma2v
  for (i in 1:q) {
    if (i+1<N) {
      diag(CN[(1+i):N, 1:(N-i)]) <- phi[i+1]*sigma2v
      diag(CN[1:(N-i), (1+i):N]) <- phi[i+1]*sigma2v
    } else if ((i+1)==N) {
      CN[(i+1), 1] <- phi[i+1]*sigma2v
      CN[1, (i+1)] <- phi[i+1]*sigma2v
    } else if ((i+1)>N) {
      break
    }
  }
  
  lambda <- eigen(CN, symmetric=T, only.values=T)[[1]]
  delta <- (sigma2v)^(1/4) * N^(-5/4) / lambda
  
  
  return(1 / (sum(delta^2)/2) * sqrt(sigma2v) )
}