IVEstSignaturePlot <- function(stock, FUN, ..., log.price=F, plot=T,
                               by.ticks=T, by.secs=F,
                               ticks=c(1:10, 20, 50, 100, 200, 500),
                               secs=c(1,2,5,10,30,60,120,300,600,1200)){

  # Draws an IV estimator signature plot of the given price process
  #
  # Args:
  #   stock: A data frame with column "TIME" (seconds) and "PRICE", or a
  #     numeric vector of price assuming equal distance sampling.
  #   FUN: A function object of an IV estimator.
  #   ...: Parameters besides stock into the estimator function
  #   by.ticks, by.secs: Logical values. Default is by.ticks.
  #   log.price: A logical value, if the stock price is log price.
  #     Default yes.
  #   ticks: (Optional) A grid of ticks to plot.
  #   secs: (Optional) A grid of secs to plot.
  #
  # Date: 2013/04/04
  
  # Check input
  if (is.numeric(stock)) {
    N <- length(stock) - 1
    stock <- data.frame(TIME=(0:N)+34200, PRICE=stock)
  }
  if (log.price) {
    stock$PRICE <- exp(stock$PRICE)
  }
  
  if (by.ticks) {

    ivest <- numeric(length(ticks))
    for (i in 1:length(ivest)) {
      ivest[i] <-
        FUN(stock[seq(from=1, to=nrow(stock), by=ticks[i]), ], ...)
    }
    if (plot) {
      X11()
      plot(ticks, ivest, type="l")
    }
  }
  if (by.secs) {
    ivest <- numeric(length(secs))
    for (i in 1:length(ivest)) {
      time.grid <- seq(from=stock$TIME[1], to=stock$TIME[nrow(stock)],
                       by=secs[i])
      ivest[i] <- FUN(stock[findInterval(time.grid, stock$TIME), ], ...)
    }
    if (plot) {
      X11()
      plot(secs, ivest, type="l")
    }
  }
  return(invisible(ivest))
}
