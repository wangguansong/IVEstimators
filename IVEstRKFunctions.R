##################################################

# Equal Unit Kernel
# Bartlett Kernel
# Second Order Kernel
# Epanechnikov Kernel
# Cubic Kernel, asymptotically equivalent to Multi-Scale estimator
# Parzen Kernel
# Tukey Hanning p Kernel
# n-th Order Kernel, n=5~8
# Generic Kernel function

# Returns:
#   A list containing dk0, dk1 (derivative of the kernel function at
#   0 and 1), k00, k11, k22 (characteristics of the kernel), cstar
#   (the optimal value of c), and rate (convergence rate).

##################################################

# Equal Unit Kernel
RKF_Unit <- function (x) {
  kernel <- list()
  kernel$weights <- rep(1, length(x))
  kernel$dk0 <- 0
  kernel$dk1 <- 0
  kernel$k00 <- 1
  kernel$k11 <- 0
  kernel$k22 <- 0
  kernel$cstar <- NA
  kernel$rate <- NA
  return(kernel)
}

# Bartlett Kernel
RKF_Bartlett <- function (x) {
  kernel <- list()
  kernel$weights <- (1-x) * (x>=0 & x<=1)
  kernel$dk0 <- -1
  kernel$dk1 <- -1
  kernel$k00 <- 1/3
  kernel$k11 <- 1
  kernel$k22 <- 0
  kernel$cstar <- 2.28
  kernel$rate <- 1/6
  return(kernel)
}

# Second Order Kernel
RKF_2ndOrder <- function (x) {
  kernel <- list()
  kernel$weights <- (1 - 2*x + x^2) * (x>=0 & x<=1)
  kernel$dk0 <- -2
  kernel$dk1 <- 0 
  kernel$k00 <- 1/5
  kernel$k11 <- 4/3
  kernel$k22 <- 4
  kernel$cstar <- 3.42
  kernel$rate <- 1/6
  return(kernel)
}

# Epanechnikov Kernel
RKF_Epanechnikov <- function (x) {
  kernel <- list()
  kernel$weights <- (1 - x^2) * (x>=0 & x<=1)
  kernel$dk0 <- 0
  kernel$dk1 <- -2
  kernel$k00 <- 8/15
  kernel$k11 <- 4/3
  kernel$k22 <- 4
  kernel$cstar <- 2.46
  kernel$rate <- 1/6
  return(kernel)
}

# Cubic Kernel, asymptotically equivalent to Multi-Scale estimator
RKF_Cubic <- function (x) {
  kernel <- list()
  kernel$weights <- (1 - 3*x^2 + 2*x^3) * (x>=0 & x<=1)
  kernel$dk0 <- 0
  kernel$dk1 <- 0
  kernel$k00 <- 0.371
  kernel$k11 <- 1.20
  kernel$k22 <- 12.0
  kernel$cstar <- 3.68
  kernel$rate <- 1/4
  return(kernel)
}

# Parzen Kernel
RKF_Parzen <- function (x) {
  kernel <- list()
  kernel$weights <- (1 - 6*x^2 + 6*x^3) * (x>=0 & x<=1/2) +
                    (2 * (1-x)^3) * (x>1/2 & x<=1)
  kernel$dk0 <- 0
  kernel$dk1 <- 0
  kernel$k00 <- 0.269
  kernel$k11 <- 1.50
  kernel$k22 <- 24.0
  kernel$cstar <- 4.77
  kernel$rate <- 1/4
  return(kernel)
}

# Tukey Hanning p Kernel
RKF_TukeyHanning <- function (x, p=2) {
  if ( ! p %in% c(1, 2, 5, 16))
    p <- 2
  kernel <- list()
  kernel$weights <- sin(pi/2 * (1-x)^p)^2 * (x>=0 & x<=1)
  kernel$dk0 <- 0
  kernel$dk1 <- 0
  switch(p,
    "1" = {
      kernel$k00 <- 0.375
      kernel$k11 <- 1.23
      kernel$k22 <- 12.1
      kernel$cstar <- 3.70
    }, 
    "2" = {
      kernel$k00 <- 0.219
      kernel$k11 <- 1.71
      kernel$k22 <- 41.7
      kernel$cstar <- 5.74
    },
    "5" = {
      kernel$k00 <- 0.097
      kernel$k11 <- 3.50
      kernel$k22 <- 489.0
      kernel$cstar <- 8.07
    },
    "16" = {
      kernel$k00 <- 0.032
      kernel$k11 <- 10.26
      kernel$k22 <- 14374.0
      kernel$cstar <- 39.16
    })
  kernel$rate <- 1/4
  return(kernel)
}

# 5th Order Kernel
RKF_5thOrder <- function (x) {
  kernel <- list()
  kernel$weights <- (1 - 10*x^3 + 15*x^4 - 6*x^5) * (x>=0 & x<=1) 
  kernel$dk0 <- 0
  kernel$dk1 <- 0
  kernel$k00 <- 0.391
  kernel$k11 <- 1.42
  kernel$k22 <- 17.1
  kernel$cstar <- 3.70
  kernel$rate <- 1/4
  return(kernel)
}
RKF_6thOrder <- function (x) {
  kernel <- list()
  kernel$weights <- (1 - 15*x^4 + 24*x^5 - 10*x^6) * (x>=0 & x<=1) 
  kernel$dk0 <- 0
  kernel$dk1 <- 0
  kernel$k00 <- 0.471
  kernel$k11 <- 1.55
  kernel$k22 <- 22.8
  kernel$cstar <- 3.97
  kernel$rate <- 1/4
  return(kernel)
}
RKF_7thOrder <- function (x) {
  kernel <- list()
  kernel$weights <- (1 - 21*x^5 + 35*x^6 - 15*x^7) * (x>=0 & x<=1)
  kernel$dk0 <- 0
  kernel$dk1 <- 0
  kernel$k00 <- 0.533
  kernel$k11 <- 1.71
  kernel$k22 <- 31.8
  kernel$cstar <- 4.11
  kernel$rate <- 1/4
  return(kernel)
}
RKF_8thOrder <- function (x) {
  kernel <- list()
  kernel$weights <- (1 - 28*x^6 + 48*x^7 - 21*x^8) * (x>=0 & x<=1)
  kernel$dk0 <- 0
  kernel$dk1 <- 0
  kernel$k00 <- 0.582
  kernel$k11 <- 1.87
  kernel$k22 <- 43.8
  kernel$cstar <- 4.31
  kernel$rate <- 1/4
  return(kernel)
}

# Generic Kernel function
RKF_Generic <- function (f, ..., H=100, rate=0, tol=1/H, rho=1) {
  #
  # Args:
  #   f: It can be a kernel function object, a list object returned by
  #     a kernel function, or a vector of weights.
  #   ...: extra arguments to the kernel function.
  #   H: discretization level when computing integration of the kernel.
  #     Default is 100.
  #   rate: the convergence rate of the kernel. Default is 0, meaning
  #     it is determined by dk0 and dk1.
  #   tol: tolerance to zero. (dk0^2+dk1^2<=tol means zero, i.e. rate is
  #     1/4.
  #   rho: IV/sqrt(QV). Default is 1.

  if (is.function(f)) {
    weights <- f( seq(from=0, to=H)/H, ... )
    if (is.list(weights)) {
      weights <- weights$weights
    }
  } else if (is.numeric(f)) {
    H <- length(f)
    weights <- c(f, 0)
  }

  dk0 <- (weights[2]-weights[1]) * H
  dk1 <- (weights[H+1]-weights[H]) * H
  first <- c(dk0, diff(weights, lag=2) * H/2, dk1)
  second <- c((first[2]-first[1]) * H,
              diff(first, lag=2) * H/2,
              (first[H+1]-first[H]) * H)
  
  k00 <- sum(((weights[-1] + weights[-(H+1)]) / 2)^2) / H
  k11 <- sum(((first[-1] + first[-(H+1)]) / 2)^2) / H
  k22 <- sum(((second[-1] + second[-(H+1)]) / 2)^2) / H
  if (rate==0) 
    rate <- ifelse((kernel$dk0^2 + kernel$dk1^2 <= tol), 1/4, 1/6)
  
  if (rate==1/6)
    cstar <- (2 * (dk0^2 + dk1^2) / k00)^(1/3)
  else if (rate==1/4)
    cstar <- sqrt( rho*k11/k00 * (1 + sqrt(1+3*k00*k22/k11^2/rho)) )
  
  kernel$weights <- weights
  kernel$dk0 <- dk0
  kernel$dk1 <- dk1
  kernel$k00 <- k00
  kernel$k11 <- k11
  kernel$k22 <- k22
  kernel$cstar <- cstar
  kernel$rate <- rate
  
  return(kernel)
}


# not tested!!
RKF_name <- function (kfname) {
  kfname <- tolower(kfname)

# Equal Unit Kernel
  if (kfname %in% c("equal", "unit", "eq"))
    kernelfunction <- RKF_Unit 
# Bartlett Kernel
  else if (kfname %in% c("bartlett", "bart", "tsrv", "twoscale"))
    kernelfunction <- RKF_Bartlett
# Second Order Kernel
  else if (kfname %in% c("2nd", "secondorder", "2ndorder"))
    kernelfunction <- RKF_2ndOrder
# Epanechnikov Kernel
  else if (kfname %in% c("epanechnikov", "epan"))
    kernelfunction <- RKF_Epanechnikov
# Cubic Kernel, asymptotically equivalent to Multi-Scale estimator
  else if (kfname %in% c("cubic", "3rd", "ms", "multiscale"))
    kernelfunction <- RKF_Cubic
# Parzen Kernel
  else if (kfname %in% c("parzen", "par"))
    kernelfunction <- RKF_Parzen
# Tukey Hanning p Kernel
  else if (kfname %in% c("th", "tukeyhanning", "tukey-hanning"))
    kernelfunction <- RKF_TukeyHanning
  else if (kfname %in% c("5th", "5thOrder"))
# n-th Order Kernel, n=5~8
    kernelfunction <- RKF_5thOrder
  else if (kfname %in% c("6th", "6thOrder"))
    kernelfunction <- RKF_6thOrder
  else if (kfname %in% c("7th", "7thOrder"))
    kernelfunction <- RKF_7thOrder
  else if (kfname %in% c("8th", "8thOrder"))
    kernelfunction <- RKF_8thOrder
  else
# Generic Kernel function
    kernelfunction <- RKF_Generic

}
