IVEstRKCore <- function (stock, weights=numeric(), flag.scale=T,
                         type="logreturn", flag.out.sample=F,
                         diagnose=F) {
  # Realized Kernel estimator
  # Args:
  #   stock: A vector or a dataframe with the first column as time stamp
  #     and the second column as data.
  #   weights: A vector of weights on AC_1, AC_2, ..., AC_H.
  #   flag.scale: T or F
  #   type: the type of data, "logreturn", "logprice", or "price".
  #   flag.out.sample: T or F
  #   diagnose: T or F
  # Returns:
  #   A number, Multi-Scale IV estimator

  ##################################################
  # Input Check

  if (is.data.frame(stock)) {
    stock <- stock[, 2]
    seconds <- stock[, 1]
  }

  # change "stock" into log-return process

  if (type %in% c("logreturn", "log-return", "lr")) {
    log.ret <- stock
  } else if (type %in% c("logprice", "log-price", "lp")) {
    log.ret <- diff(stock)
  } else if (type %in% c("price", "p")) {
    log.ret <- diff(log(stock))
  }

  N <- length(log.ret)
  if (exists("seconds", mode="numeric")) {
    time.interval <- diff(seconds)
  } else {
    time.interval <- rep(1/length(log.ret), N)
  }
  H <- length(weights)
  if (N<2 | N<H) {
    return(NA)
  }


  ##################################################
  # create autocovariations
  Gamma2 <- numeric(length=H+1)
  

  if (flag.out.sample | flag.scale) {
    for (i in 0:H) {
      Gamma2[i+1] <- sum(log.ret[(H+1):(N-H)] *
                        log.ret[(H+1+i):(N-H+i)]) +
                    sum(log.ret[(H+1):(N-H)] *
                        log.ret[(H+1-i):(N-H-i)])
    }
  } else {
    for (i in 0:H) {
      Gamma2[i+1] <- 2 * sum( log.ret[(1+i):N] * log.ret[1:(N-i)] )
    }
  }
  Gamma2[1] <- Gamma2[1] / 2
  
  # scale factor
  scale.factor <- 1
  if (flag.scale & H>0 & flag.out.sample==F) {
    scale.factor <- sum(time.interval) /
                      sum(time.interval[(H+1):(N-H)])
  }
  

  # Estimator
  value <- Gamma2[1] + ifelse(H>0, sum(Gamma2[-1] * weights), 0)
  value <- scale.factor * value

  ##################################################
  # return value
  if (diagnose) {
    estimator <- list()
    estimator$value <- value
    estimator$Gamma <- c(Gamma2[1], Gamma2[-1]/2)
    estimator$bandwidth <- H
    estimator$weights <- weights
    estimator$out.sample <- out.sample
    estimator$scale <- scale.factor
    return(estimator)
  } else
    return(value)
  

}
