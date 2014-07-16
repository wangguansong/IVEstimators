IVEstTwoScaleOpt <- function(stock, type="logreturn", noise.var=0,
                             pre.var=0, bias.adjusted=TRUE) {
  # Two Scale estimator of integrated volatility, with optimal K (the
  # second time scale)
  #
  # Args:
  #   stock: A vector or a dataframe with the first column as time stamp
  #     and the second column as data.
  #   type: The type of data, "logreturn", "logprice", or "price".
  #   noise.var: Variance of the additive noise. Default is 0 which is
  #     to be estimated using data.
  #   pre.var: Preliminary estimate of variation. Default is 0 which is
  #     to be estimmated using sub-sampled data, assuming 6.5 hours and
  #     in 15 minutes frequency.
  #   bias.adjusted: boolean.
  # Returns:
  #   A number, two scale estimator
  #
  # Reference:
  # Zhang, Mykland, Ait-Sahalia (2005), A Tale of Two Time Scales: 
  # Determining Integrated Volatility With Noisy High-Frequency Data,
  # Journal of the American Statistical Association
  # Optimal K is given as c*N^(2/3), where optimal c is given in
  # equation (63). Note: assuming constant volatility.
  #
  ########## Guansong Wang 10/30/2013 ##########
  ########## wang.guansong@hotmail.com ##########

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

  # compute optimal scale
  if (noise.var==0) {
    noise.var <- sum(diff(stock)^2) / 2 / (N-1)
  }
  if (pre.var==0) {
    if (!exists("seconds", mode="numeric")) {
      seconds <- seq(from=6.5*3600, to=16*3600, length.out=N)
    }
    grids.15min <- seq(from=6.5*3600, to=16*3600, by=15*60)
    index.15min <- findInterval(grids.15min, seconds)
    index.15min[which(index.15min==0)] <- 1
    pre.var <- sum(diff(stock[index.15min])^2)
  }
  scale <- (12 * noise.var^2 / pre.var^2)^(1/3) * (N-1)^(2/3)
  scale <- round(scale)
  if (scale==0) {
    scale <- 1
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
