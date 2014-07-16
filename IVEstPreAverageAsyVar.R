IVEstPreAverageAsyVar <- function(stock, bandwidth, type="logreturn",
                                  adj.flag=TRUE, numeric.par) {
  # Compute estimated asymptotic variance
  # Pre-Averaging approach, with simple weight function min(x, 1-x)
  #
  # Reference:
  # Jacod, Li, Mykland, Podolskij, and Vetter (2009), Microstructure
  # Noise in the Continuous Case: The Pre-averaging Approach, Stochastic
  # Processes and Their Applications
  #
  # Args:
  #   stock: A vector or a dataframe with the first column as time stamp
  #     and the second column as data.
  #   bandwidth: Preaveraging bandwith, (k_n in the paper)
  #   type: the type of data, "logreturn", "logprice", or "price".
  #   adj.flag: Boolean, adjusted or not.
  #   numeric.par: A vector of parameter values in finite cases in order
  #     of (psi1, psi2, Phi11, Phi12, Phi22). If not provided, they will
  #     be computed.
  # Returns:
  #   A number, estimator of asymptotic variance of PA estimator
  ##################################################
  # Input Check
  if (is.data.frame(stock)) {
    stock <- stock[, 2]
  }
  # change "stock" into log-return process
  if (type %in% c("logreturn", "log-return", "lr")) {
    log.ret <- stock
  } else if (type %in% c("logprice", "log-price", "lp")) {
    log.ret <- diff(stock)
  } else if (type %in% c("price", "p")) {
    log.ret <- diff(log(stock))
  }
  N <- length(log.ret)    # number of log-returns
  if (N<2 | N<bandwidth-1) {
    return(NA)
  }

  ##################################################
  # weight function
  g.fun <- function(x) {pmin(x, (1-x))}
  if (adj.flag & exists("numeric.par", mode="numeric")) {
    psi1 <- numeric.par[1]
    psi2 <- numeric.par[2]
    Phi11 <- numeric.par[3]
    Phi12 <- numeric.par[4]
    Phi22 <- numeric.par[5]
  } else if (adj.flag & !exists("numeric.par", mode="numeric")) {
    grids <- 0:bandwidth/ bandwidth
    g.value <- g.fun(grids)
    g.1stdev <- c(g.value[2]-g.value[1],
                  (g.value[3:(bandwidth+1)]-g.value[1:(bandwidth-1)])/2,
                  g.value[bandwidth+1]-g.value[bandwidth]) * bandwidth
    psi1 <- sum((g.1stdev[1:bandwidth]^2+g.1stdev[2:(bandwidth+1)]^2)/
                2) / bandwidth
    psi2 <- sum((g.value[1:bandwidth]^2 + g.value[2:(bandwidth+1)]^2)/
                2) / bandwidth
    phi1 <- numeric(bandwidth+1) 
    phi2 <- numeric(bandwidth+1) 
    for (j in 0:bandwidth) {
      phi1[j+1] <- sum(g.1stdev[j:bandwidth+1]*
                       g.1stdev[j:bandwidth+1-j]) / bandwidth
      phi2[j+1] <- sum(g.value[j:bandwidth+1]*
                       g.value[j:bandwidth+1-j]) / bandwidth
    }
    Phi11 <- sum((phi1[1:bandwidth]^2 + phi1[2:(bandwidth+1)]^2) / 2) /
              bandwidth
    Phi22 <- sum((phi2[1:bandwidth]^2 + phi2[2:(bandwidth+1)]^2) / 2) /
              bandwidth
    Phi12 <- sum(((phi1*phi2)[1:bandwidth]+(phi1*phi2)[2:(bandwidth+1)])
                 /2) / bandwidth
  } else {
    psi1 <- 1
    psi2 <- 1/12
    Phi11 <- 1/6
    Phi12 <- 1/96
    Phi22 <- 151/80640
  }
  theta <- bandwidth / sqrt(N)
  ##################################################
  # Create pre-averaged log-returns
  pre.avg.logret <- filter(x=log.ret,
                            filter=g.fun((bandwidth-1):1/bandwidth),
                            sides=1)
  pre.avg.logret <- as.numeric(pre.avg.logret)
  pre.avg.logret <- pre.avg.logret[!is.na(pre.avg.logret)]

  ##################################################
  # compute estimator
  part1 <- 4*Phi22/3/theta/psi2^4 * sum(pre.avg.logret^4)
  pre.avg.logret2 <- filter(x=log.ret[-(1:(bandwidth-1))]^2,
                            filter=rep(1, bandwidth),
                            sides=1)
  pre.avg.logret2 <- as.numeric(pre.avg.logret2)
  pre.avg.logret2 <- pre.avg.logret[!is.na(pre.avg.logret2)]
  part2 <- 4/N/theta^3 * (Phi12/psi2^3 - Phi22*psi1/psi2^4) *
           sum(pre.avg.logret[1:(N-2*bandwidth+2)]^2*pre.avg.logret2)
  part3 <- 1/N/theta^3 * (Phi11/psi2^2 - 2*Phi12*psi1/psi2^3 +
                          Phi22*psi1^2/psi2^4) *
           sum(log.ret[-c(N-1,N)]^2 * log.ret[-c(1,2)]^2)
  if (adj.flag) {
    return(1 / (1 - psi1/2/bandwidth^2/psi2)^2 *
           part1 * N / (N-bandwidth+2) +
           part2 * N / (N-bandwidth+2) +
           part3 * N / (N-2))
  } else {
    return(part1 + part2 + part3)
  }
}
