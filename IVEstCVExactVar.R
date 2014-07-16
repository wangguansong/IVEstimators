# memory demanding, crash
IVEstCVExactVar <- function(K, beta, sigma2.path, sigma2v, var.v2=3*sigma2v, 
                            diagnose=F) {
  N <- length(sigma2.path)
  dt <- 1/N
  q <- length(beta) + 1
  theta <- c(1, diff(c(1, beta, 0)))
  
  # IV, iQ
  iv <- sum(sigma2.path * dt)
  iq <- sum(sigma2.path^2 * dt)
  
  matrixT <- matrix(0, nrow=N, ncol=N+q)
  for (i in 1:N) {
    matrixT[i, 0:q+i] <- theta
  }
  
  matrixA <- matrix(0, nrow=K, ncol=K)
  matrixB <- matrix(0, nrow=K, ncol=K)
  matrixNCD <- matrix(0, nrow=K, ncol=K)
  matrixNEF <- matrix(0, nrow=K, ncol=K)
  
  for (nr in 1:K) {
    for (nc in 1:K) {
      k <- nr - 1
      l <- nc - 1
      if (abs(k-l)>2*q) {
        next
      } else {
        matrixJk <- matrix(0, nrow=N, ncol=N)
        if (k+1<N) {
          diag(matrixJk[(1+k):N, 1:(N-k)]) <- 1
          diag(matrixJk[1:(N-k), (1+k):N]) <- 1
        } else if ((k+1)==N) {
          matrixJk[k+1, 1] <- 1
          matrixJk[1, k+1] <- 1
        }
        matrixJl <- matrix(0, nrow=N, ncol=N)
        if (l+1<N) {
          diag(matrixJl[(1+l):N, 1:(N-l)]) <- 1
          diag(matrixJl[1:(N-l), (1+l):N]) <- 1
        } else if ((l+1)==N) {
          matrixJl[l+1, 1] <- 1
          matrixJl[1, l+1] <- 1
        }
        
        # matrix A
        if (k==0 & l==0) {
          matrixA[nr, nc] <- 1
        } else if (k==l) {
          matrixA[nr, nc] <- 2 / iq * dt *
            sum(sigma2.path[1:(N-k)]*sigma2.path[(k+1):N]) 
        }
        
        # matrix B
        matrixB[nr, nc] <- 1 / iv / 4 *
          sum(sigma2.path * diag(matrixJk%*%matrixT%*%t(matrixT)%*%matrixJl)) 

        # matrix N*E+F
        matrixNEF[nr, nc] <- sum(diag(t(matrixT)%*%matrixJk%*%matrixT) *
                                   diag(t(matrixT)%*%matrixJl%*%matrixT))
        
        # matrix N*C+D
        matrixNCD[nr, nc] <- (sum((t(matrixT)%*%matrixJk%*%matrixT)*
                                    (t(matrixT)%*%matrixJl%*%matrixT)) -
                                matrixNEF[nr, nc]) / 2

      }
    }
  }
  if (diagnose) {
    return(list(var.matrix=2*iq/N * matrixA + 4*sigma2v*iv * matrixB +
                  4*sigma2v^2 * matrixNCD + var.v2 * matrixNEF,
                A=matrixA,
                B=matrixB,
                NCD=matrixNCD,
                NEF=matrixNEF))
                
  } else {
    return(2*iq/N * matrixA + 4*sigma2v*iv * matrixB +
             4*sigma2v^2 * matrixNCD + var.v2 * matrixNEF)
  }
  
  
}