IVEstTwoScale <- function(stock, scale, type="logreturn",
                          bias.adjusted=TRUE) {
  # Two Scale estimator of integrated volatility
  #
  # Args:
  #   stock: A vector of prices or a dataframe with the first column as
  #     time stamp and the second column as data.
  #   scale: a number, indicating the second time scale K.
  #   type: the type of data, "logreturn", "logprice", or "price".
  #   bias.adjusted: boolean.
  # Returns:
  #   A number, two scale estimator
  #
  # Reference:
  # Zhang, Mykland, Ait-Sahalia (2005), A Tale of Two Time Scales: 
  # Determining Integrated Volatility With Noisy High-Frequency Data,
  # Journal of the American Statistical Association
  #
  # Note:
  #   When scale is 1, it returns a NaN.

  if (is.data.frame(stock)) {
    stock <- stock[, 2]
    seconds <- stock[, 1]
  }

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
  subsample.rv <- sum((stock[(scale+1):N]-stock[1:(N-scale)])^2) / scale
  rv <- sum(diff(stock)^2)

  # return (bias adjusted) estimator
  N.bar <- (N-scale+1) / scale
  ivest <- subsample.rv - N.bar/N*rv
  if (bias.adjusted) {
    return(ivest/(1-N.bar/N))
  } else {
    return(ivest)
  }


}
