# Surplus production to estimate reference point under non-equilibrium assumption for Schaefer using time series data

rm(list = ls())

### Create fake data from known parameters
r <- 0.2
K <- 1000
q <- 0.00025
Bo <- K
maxEffort <- 100
nYears <- 20
effort <- c(seq(1,500, length.out = nYears/2), rev(seq(1,500, length.out = nYears / 2)))

B <- CPUE <- C <- rep(NA, nYears)
procError <- 0.05
catchError <- 0.05
for (i in 1:nYears) {
  if (i == 1) B[i] <- Bo
  if (i>1) B[i] <- B[i-1] + rlnorm(1, meanlog = -procError^2/2, sdlog = procError) * B[i-1] * r * (1 - B[i-1]/K) - C[i-1]
  C[i] <- q * effort[i] * B[i] * rlnorm(1, meanlog = -catchError^2/2, sdlog = catchError)
  CPUE[i] <- C[i] / effort[i]
}  

df <- data.frame(catch = C,
                 effort = effort)

MSY <- r*K/4
Emsy <- r/2*q
inpars <- c(MSY, Bo, Emsy)

### create Schaefer function and plot for data fitting
RPpar.S <- function(inpars, df){
  MSY <- exp(inpars[1])
  B0 <- exp(inpars[2])
  Emsy <- exp(inpars[3])
  
  r <- (MSY*4)/k
  k <- (MSY*4)/r
  q <- r/(2*Emsy)

  CPUE <- df$catch/df$effort
  EstCatch <- EstCPUE <- EstBt <- vector(length=nrow(df))
  
  for (i in 1:nYears) {
    if (i == 1) EstBt[i] <- B0
    if (i>1) EstBt[i] <- EstBt[i-1] +  EstBt[i-1] * r * (1 - EstBt[i-1]/k) - df$catch[i-1]
    EstCatch[i] <- EstBt[i] * q * df$effort[i]
  }
  
  EstCPUE <- EstCatch / df$effort
  
  plot(CPUE, xlab="Year", ylab="CPUE", type="b", col="blue")
  lines(EstCPUE, type="b", col="red")
  legend("topright", legend=c("Observation", "Estimation"), lty=c(1, 1), col=c("blue", "red"), box.lty=1, cex=0.7)
  
  return(data.frame(CPUE=CPUE, EstBt=EstBt, EstCatch=EstCatch, EstCPUE=EstCPUE))
}

### create minimizing function
RPpar.S_opt <- function(inpars, df){
  MSY <- exp(inpars[1])
  B0 <- exp(inpars[2])
  Emsy <- exp(inpars[3])
  sigma <- exp(inpars[4])
  
  r <- (MSY*4)/k
  k <- (MSY*4)/r
  q <- r/(2*Emsy)
  
  CPUE <- df$catch/df$effort
  EstCatch <- EstCPUE <- EstBt <- vector(length=nrow(df))
  
  for (i in 1:nYears) {
    if (i == 1) EstBt[i] <- B0
    if (i>1) EstBt[i] <- max(c(0.01,
                               EstBt[i-1] +  EstBt[i-1] * r * (1 - EstBt[i-1]/k) - df$catch[i-1]
    )
    )
  }
  EstCPUE <-  EstBt * q
  
  nll <- -sum(dlnorm(x= na.omit(CPUE), meanlog = log(na.omit(EstCPUE)), sdlog = sigma, log = TRUE))
  return(nll)
}

### Estimate parameters using optim
startPars <- c(log(MSY), log(Bo), log(Emsy), log(0.1))

fit <- optim(par=startPars, 
             fn=RPpar.S_opt, 
             df=df, 
             method="Nelder-Mead")

fitted_pars <- exp(fit$par)
cbind(fitted_pars, c(MSY, Bo, Emsy, NA))

### plot fit
predicted <- RPpar.S(inpars = fitted_pars[1:4], df)


