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
plot(CPUE, type="l")

Bmsy <- K/2
MSY <- (r*K)/4
Emsy <- r/(2*q)
inpars <- c(Bmsy, MSY, Emsy, Bo)

### create Schaefer function and plot for data fitting
RPpar.S <- function(inpars, df){
  Bmsy <- inpars[1]
  MSY <- inpars[2]
  Emsy <- inpars[3]
  B0 <- inpars[4]
  
  k <- 2*Bmsy
  r <- (MSY*4)/k
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
  Bmsy <- exp(inpars[1])
  MSY <- exp(inpars[2])
  Emsy <- exp(inpars[3])
  B0 <- exp(inpars[4])
  sigma <- exp(inpars[5])
  
  k <- 2*Bmsy
  r <- (MSY*4)/k
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
  
  nll <- -sum(dlnorm(x= CPUE, meanlog = log(EstCPUE), sdlog = sigma, log = TRUE))
  return(nll)
}

### Estimate parameters using optim
startPars <- c(log(Bmsy), log(MSY), log(Emsy), log(Bo), log(0.1))

fit <- optim(par=startPars, 
             fn=RPpar.S_opt, 
             df=df, 
             method="Nelder-Mead",
             hessian=F)

fitted_pars <- exp(fit$par)

cbind(fitted_pars, c(Bmsy, MSY, Emsy, Bo, NA))

### plot fit
predicted <- RPpar.S(inpars = fitted_pars[1:4], df)
predicted

### create likelihood profile for Bmsy, MSY, Emsy
Bmsy_Sprofiler <- function(inpars, inpBmsy, df){
  Bmsy <- exp(inpBmsy)
  MSY <- exp(inpars[2])
  Emsy <- exp(inpars[3])
  B0 <- exp(inpars[4])
  sigma <- exp(inpars[5])
  
  k <- 2*Bmsy
  r <- (MSY*4)/k
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
  
  nll <- -sum(dlnorm(x= CPUE, meanlog = log(EstCPUE), sdlog = sigma, log = TRUE))
  return(nll)
}

MSY_Sprofiler <- function(inpars, inpMSY, df){
  Bmsy <- exp(inpars[1])
  MSY <- exp(inpMSY)
  Emsy <- exp(inpars[3])
  B0 <- exp(inpars[4])
  sigma <- exp(inpars[5])
  
  k <- 2*Bmsy
  r <- (MSY*4)/k
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
  
  nll <- -sum(dlnorm(x= CPUE, meanlog = log(EstCPUE), sdlog = sigma, log = TRUE))
  return(nll)
}

Emsy_Sprofiler <- function(inpars, inpEmsy, df){
  Bmsy <- exp(inpars[1])
  MSY <- exp(inpars[1])
  Emsy <- exp(inpEmsy)
  B0 <- exp(inpars[4])
  sigma <- exp(inpars[5])
  
  k <- 2*Bmsy
  r <- (MSY*4)/k
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
  
  nll <- -sum(dlnorm(x= CPUE, meanlog = log(EstCPUE), sdlog = sigma, log = TRUE))
  return(nll)
}

Bmsy.vec <- log(seq(from=300, to=700, by=1))
Bmsy.NLL <- vector(length=length(Bmsy.vec))
for (i in 1:length(Bmsy.vec)){
  tryCatch({
    x <- optim(par=startPars, fn=Bmsy_Sprofiler, inpBmsy=Bmsy.vec[i],
               df=df, method="Nelder-Mead")
    Bmsy.NLL[i] <- x$value
  }, error=function(e){cat("Error :", conditionMessage(e), "\n")}
  )
}

MSY.vec <- log(seq(from=30, to=70, by=0.01))
MSY.NLL <- vector(length=length(MSY.vec))
for (i in 1:length(MSY.vec)){
  tryCatch({
    x <- optim(par=startPars, fn=MSY_Sprofiler, inpMSY=MSY.vec[i],
             df=df, method="Nelder-Mead")
  MSY.NLL[i] <- x$value
  }, error=function(e){cat("Error :", conditionMessage(e), "\n")}
  )
}

Emsy.vec <- log(seq(from=100, to=800, by=1)) # unfinished code
Emsy.NLL <- vector(length=length(Emsy.vec))
for (i in 1:length(Emsy.vec)){
  tryCatch({
    x <- optim(par=startPars, fn=Emsy_Sprofiler, inpEmsy=Emsy.vec[i],
               df=df, method="Nelder-Mead")
    Emsy.NLL[i] <- x$value
  }, error=function(e){cat("Error :", conditionMessage(e), "\n")}
  )
}


cbind(exp(Bmsy.vec), Bmsy.NLL)
cbind(exp(MSY.vec), MSY.NLL)
cbind(exp(Emsy.vec), Emsy.NLL)

par(mfrow=c(3,1))
plot(Bmsy.NLL, type="l")
plot(MSY.NLL, type="l")
plot(Emsy.NLL, type="l")
