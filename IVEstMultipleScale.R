IVEstMultipleScale <- function(stock, scale, type="logreturn") {
  # Multiple Scale estimator to integrated volatility
  # Args:
  #   stock: A vector or a dataframe with the first column as time stamp
  #     and the second column as data.
  #   scale: a number, indicating the second time scale K.
  #   type: the type of data, "logreturn", "logprice", or "price".
  # Returns:
  #   A number, Multi-Scale IV estimator
  
  
  if (is.data.frame(stock)) {
    stock <- stock[, 2]
    seconds <- stock[, 1]
  }
  # TODO: check "seconds"
  # change "stock" into log-price process
  if (type %in% c("logreturn", "log-return", "lr")) {
    stock <- c(0, cumsum(stock))
  } else if (type %in% c("logprice", "log-price", "lp")) {
  } else if (type %in% c("price", "p")) {
    stock <- log(stock)
  }
  N <- length(stock)
  if (N<2 | N<scale) {
    return(NA)
  }
  
  # compute rv and sub-sample rv:
  subsample.rv <- numeric(scale)
  for (i in 1:scale) {
    subsample.rv[i] <- sum((stock[(i+1):N]-stock[1:(N-i)])^2) / i
  }
  
  a <- 12 * (1:scale)/scale^2 * ((1:scale)/scale - 1/2 - 1/2/scale) /
       (1 - 1/scale^2)
  a[1] <- a[1] + 2/N
  a[2] <- a[2] - 2/N
  
  return(sum(a*subsample.rv))
  
}
