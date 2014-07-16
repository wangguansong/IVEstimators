IVEstTS <- function(stock, K = NA) {
  # Time series
  #
  # Arguments:
  #   stock
  #   K
  #
  # Returns:
  #   ivest
  #
  # Arthor: Guansong Wang
  # Updates: 04/05/2013
  # 
  if (is.data.frame(stock)) {
    log.ret <- diff(log(stock[, 2]))
  } else {
    log.ret <- diff(log(stock))
  }
  N <- length(log.ret)
  if (N < 2) {
    cat("Warning from IVEstTS:\n Check class of log.returns")
  }
  Nbar <- (N-K+1) / K

  qv1 <- sum(log.ret^2)
  log.ret.K <- filter(x = log.ret, filter = rep(1,K),
                          method = "convolution", sides = 1)
  log.ret.K <- log.ret.K[!is.na(log.ret.K)]
  qvK <- sum(log.ret.K^2) / K

  iv.est <- (qvK - Nbar/N*qv1) / (1-Nbar/N)
  
  return(iv.est)
}
