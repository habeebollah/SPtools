function(inpars, df){
  MSY <- exp(inpars[1])
  B0 <- exp(inpars[2])
  Emsy <- exp(inpars[3])
  sigma <- exp(inpars[4])

  r*k/4 <- MSY
  r/2*q <- Emsy

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
  #nll <- sum(0.5*log(2*pi)+log(sigma)+log(CPUE)+(log(CPUE)-log(EstCPUE))^2/(2*sigma^2))
  return(nll)
}
