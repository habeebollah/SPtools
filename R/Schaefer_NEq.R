function(inpars, df){
  k <- inpars[1]
  B0 <- inpars[2]
  r <- inpars[3]
  q <- inpars[4]
  
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