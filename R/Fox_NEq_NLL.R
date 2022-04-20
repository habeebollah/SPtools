function(inpars, df){
  k <- exp(inpars[1])
  B0 <- exp(inpars[2])
  r <- exp(inpars[3])
  q <- exp(inpars[4])
  sigma <- exp(inpars[5])
  
  CPUE <- df$catch/df$effort
  EstCatch <- EstCPUE <- EstBt <- vector(length=nrow(df))
  
  for (i in 1:nYears) {
    if (i == 1) EstBt[i] <- B0
    if (i>1) EstBt[i] <- max(c(0.01,
                               EstBt[i-1] +  EstBt[i-1] * r * (1 - log(EstBt[i-1])/log(k)) - df$catch[i-1]
    )
    )
  }
  EstCPUE <-  EstBt * q
  
  nll <- -sum(dlnorm(x= CPUE, meanlog = log(EstCPUE), sdlog = sigma, log = TRUE))
  
  return(nll)
}