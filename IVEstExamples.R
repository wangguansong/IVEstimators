freq <- 23400
Tend <- 1
mu <- 0
V <- .2/sqrt(252)
p0 <- 1   # log price
leverage <- 0
dWc <- 0
sigma.seed <- 101

ma.parameter <- .5
ma.parameter <- 1
var.noise <- 4.2e-8
alpha <- -1/40
beta1 <- 1/8
beta0 <- beta1^2/(2*alpha)


# sv1f
sigma <- GenSigmaSV1F(freq=freq, V=V, Tend=Tend, p0=p0, seed=sigma.seed, 
                      alpha=alpha, beta1=beta1, leverage=leverage, diagnose=T)
x.seed <- round((as.numeric(Sys.time()) - floor(as.numeric(Sys.time()))) * 1e5)
y <- exp(AddLogNoise(GenLatLogPrice(freq=freq, mu=mu, sigma=sigma$sigma,
                                    Tend=Tend, p0=p0, dWc=sigma$dWc,
                                    leverage=sigma$leverage, seed=x.seed),
                     ma.parameter=ma.parameter, var.noise=var.noise))

cat("Real Value \n")
print(sum(sigma$sigma^2/freq))
plot(y, type="l")

cat("Preliminary IV Estimate, 20 minutes \n")
N <- length(y) - 1
seconds <- (0:N)/N * 23400 + 9.5*3600
time.grid <- seq(from=seconds[1], to=seconds[N+1],
                 by=20*60)
pre.ivest <- sum(diff(log(y[findInterval(time.grid, seconds)]))^2)
print(pre.ivest)


cat("Realized Variation \n")
print(IVEstRKCore(y))

cat("Magnitude of noise \n")
omega2 <- (IVEstRKCore(y)-pre.ivest) / 2 / N
print(omega2)

cat("Realized Variation, AC1 \n")
print(IVEstRKCore(y, weights=c(1)))

cat("Realized Variation, AC2 \n")
print(IVEstRKCore(y, weights=c(1,1)))

cat("Realized Variation, AC3 \n")
print(IVEstRKCore(y, weights=c(1,1,1)))

cat("Realized Kernel, TH2 \n")

xi <- sqrt(omega2 / pre.ivest)
th2weights <- RKF_TukeyHanning(0, 2)
Hstar <- round(th2weights$cstar * xi * sqrt(N))
th2weights <- RKF_TukeyHanning( (0:(Hstar-1) / Hstar), p=2)$weights
print(IVEstRKCore(y, th2weights, scale=T))


cat("Realized Kernel, Parzen 1/5 \n")

xi <- sqrt(omega2 / pre.ivest)
Hstar <- round(3.5134 * xi^(4/5) * N^(3/5))
parzen.weights <- RKF_Parzen( (0:(Hstar-1) / Hstar))$weights
print(IVEstRKCore(y, parzen.weights, scale=T))


cat("Realized Kernel, CV \n")
H <- round(Hstar)
q <- 1
cvweights <- IVEstCVWeights2(beta=1, H=Hstar, ivest=pre.ivest, sigma2v=omega2,
                             N=N)
print(IVEstRKCore(y, c(rep(1, q), cvweights), scale=F))

cat("Realized Kernel, bias corrected, TH2 \n")
q <- length(ma.parameter) + 1
xi <- sqrt(omega2 / pre.ivest)
th2weights <- RKF_TukeyHanning(0, 2)
Hstar <- round(th2weights$cstar * xi * sqrt(N))
th2weights <- RKF_TukeyHanning( (0:(Hstar-1) / Hstar), p=2)$weights
print(IVEstRKCore(y, c(rep(1,q), th2weights), scale=T))

cat("Realized Kernel, bias corrected, CV \n")
H <- Hstar
q <- length(ma.parameter) + 1
cvweights <- IVEstCVWeights(y, q, H, T)
print(IVEstRKCore(y, c(rep(1, q), cvweights), scale=F))
