# Surplus production under non-equilibrium assumption for Fox using time series data

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
  if (i == 1) B[i] <- log(Bo)
  if (i>1) B[i] <- B[i-1] + rlnorm(1, meanlog = -procError^2/2, sdlog = procError) * B[i-1] * r * (1 - log(B[i-1])/log(K)) - C[i-1]
  C[i] <- q * effort[i] * B[i] * rlnorm(1, meanlog = -catchError^2/2, sdlog = catchError)
  CPUE[i] <- C[i] / effort[i]
}  

df <- data.frame(catch = C,
                 effort = effort)

inpars <- c(K, Bo, r, q)

### create Fox function and plot for data fitting
SPpar.F <- function(inpars, df){
  k <- inpars[1]
  B0 <- inpars[2]
  r <- inpars[3]
  q <- inpars[4]
  
  CPUE <- df$catch/df$effort
  EstCatch <- EstCPUE <- EstBt <- vector(length=nrow(df))
  
  for (i in 1:nYears) {
    if (i == 1) EstBt[i] <- log(B0)
    if (i>1) EstBt[i] <- EstBt[i-1] +  EstBt[i-1] * r * (1 - log(EstBt[i-1])/log(k)) - df$catch[i-1]
    EstCatch[i] <- EstBt[i] * q * df$effort[i]
  }
  
  EstCPUE <- EstCatch / df$effort
  
  plot(CPUE, xlab="Year", ylab="CPUE", type="b", col="blue")
  lines(EstCPUE, type="b", col="red")
  legend("topright", legend=c("Observation", "Estimation"), lty=c(1, 1), col=c("blue", "red"), box.lty=1, cex=0.7)
  
  return(data.frame(CPUE=CPUE, EstBt=EstBt, EstCatch=EstCatch, EstCPUE=EstCPUE))
}

### create minimizing function
SPpar.F_opt <- function(inpars, df){
  k <- exp(inpars[1])
  B0 <- exp(inpars[2])
  r <- exp(inpars[3])
  q <- exp(inpars[4])
  sigma <- exp(inpars[5])
  
  CPUE <- df$catch/df$effort
  EstCatch <- EstCPUE <- EstBt <- vector(length=nrow(df))
  
  for (i in 1:nYears) {
    if (i == 1) EstBt[i] <- log(B0)
    if (i>1) EstBt[i] <- max(c(0.01,
                               EstBt[i-1] +  EstBt[i-1] * r * (1 - log(EstBt[i-1])/log(k)) - df$catch[i-1]
    )
    )
  }
  EstCPUE <-  EstBt * q
  
  nll <- -sum(dlnorm(x= na.omit(CPUE), meanlog = log(na.omit(EstCPUE)), sdlog = sigma, log = TRUE))
  
  return(nll)
}

### Estimate parameters using optim
startPars <- c(log(K), log(Bo), log(r), log(q), log(0.1))

fit <- optim(par=startPars, 
             fn=SPpar.F_opt, 
             df=df, 
             method="Nelder-Mead")

fitted_pars <- exp(fit$par)
cbind(fitted_pars, c(K, Bo, r, q, NA))

### plot fit
predicted <- SPpar.F(inpars = fitted_pars[1:4], df)


