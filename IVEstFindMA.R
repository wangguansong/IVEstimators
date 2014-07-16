IVEstFindMA <- function(yacf, q=length(yacf), ret="theta") {
  # Find MA parameters from Gamma
  
  yacf <- yacf[1:q]
  
  f <- function(x, acf1) {
    q <- length(acf1)
    
    theta <- c(1, x, -(1+sum(x)))
    phi <- numeric(q+1)
    for (i in 0:q) {
      phi[i+1] <- sum(theta[0:(q-i)+1+i] * theta[0:(q-i)+1])
    }
    acf2 <- phi[-1] / phi[1]
    return(sum((acf1 - acf2)^2))
  }
  x<-nlm(f, c(rep(-1/q, q-1)), yacf, gradtol=1e-10)
  if (ret=="theta") {
    return(c(1, x$estimate, -1-sum(x$estimate)))
  } else if (ret=="nlm") {
    return(x)
  } else if (ret=="beta") {
    return(cumsum(c(1, x$estimate)))
  }
  return(x)
  

}