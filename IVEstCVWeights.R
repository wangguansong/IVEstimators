IVEstCVWeights <- function(stock, q, H, adjust=T, diagnose=F) {
  
  if (H==0) {
    return(numeric())
  }
  
  if (is.data.frame(stock)) {
    seconds <- stock[, 1]
    prices <- stock[, 2]
    N <- nrow(stock) - 1
  } else if (is.vector(stock, mode="numeric")) {
    prices <- stock
    N <- length(stock) - 1
    seconds <- (0:N)/N * 23400 + 9.5*3600
  }
  
  # preliminary estimator
  time.grid <- seq(from=seconds[1], to=seconds[N+1],
                   by=20*60)
  pre.ivest <- sum(diff(log(prices[findInterval(time.grid, seconds)]))^2)
  
  log.ret <- diff(log(prices))
  sigma2.phi <- numeric(length=q+1)
  
  for (i in 0:q) {
    sigma2.phi[i+1] <- mean( log.ret[(1+i):N] * log.ret[1:(N-i)] )
  }
  sigma2.phi[1] <- sigma2.phi[1] - pre.ivest/N
  
  # matrices
  K <- 1 + q + H
  # A
  matrixA <- diag(2, nrow=K, ncol=K)
  matrixA[1] <- 1
  
  # major diagonal of B and C
  sigma2.phi.long <- c(rev(sigma2.phi[-1]), sigma2.phi)
  matrixB <- diag(sigma2.phi[1]*2, nrow=K, ncol=K)
  matrixB[1] <- matrixB[1] / 2
  matrixC <- diag(sum(sigma2.phi.long^2), nrow=K, ncol=K)

  for (i in 1:q) {
    psi <- sum(sigma2.phi.long[(i+1):(1+2*q)] * sigma2.phi.long[1:(1+2*q-i)])
    if ((i+1)<K) {
      #B
      diag(matrixB[(i+1):K, 1:(K-i)]) <- sigma2.phi[i+1]*2
      diag(matrixB[1:(K-i), (i+1):K]) <- sigma2.phi[i+1]*2
      
      #C
      
      diag(matrixC[(i+1):K, 1:(K-i)]) <- psi
      diag(matrixC[1:(K-i), (i+1):K]) <- psi
    } else if ((i+1)==K) {
      #B
      matrixB[i+1, 1] <- sigma2.phi[i+1]*2
      matrixB[1, i+1] <- sigma2.phi[i+1]*2
      #C
      matrixC[i+1, 1] <- psi
      matrixC[1, i+1] <- psi
    }
  }
  

  for (i in (q+1):min(2*q, q+H)) {
    psi <- sum(sigma2.phi.long[(i+1):(1+2*q)] * sigma2.phi.long[1:(1+2*q-i)])
    if ((i+1)<K) {

    diag(matrixC[(i+1):K, 1:(K-i)]) <- psi
    diag(matrixC[1:(K-i), (i+1):K]) <- psi
    } else if ((i+1)==K) {
      matrixC[i+1, 1] <- psi
      matrixC[1, i+1] <- psi
    }
  }
  
  var.matrix <- 2*pre.ivest^2/N*matrixA + 4*pre.ivest*matrixB + 4*N*matrixC

  if (adjust) {
    adj.vector <- N /  (N- 0:q)
  } else {
    adj.vector <- rep(1, q+1)
  }
  

  
  if (diagnose) {
    return(list(var=var.matrix,
                A=matrixA,
                B=matrixB,
                C=matrixC))
  } else {
    return(as.numeric(
      -solve(var.matrix[(2+q):K, (2+q):K, drop=F], 
             t(t(var.matrix[(2+q):K, 1:(1+q), drop=F])*adj.vector)
             %*% matrix(rep(1, q+1), ncol=1))))
  }

}
