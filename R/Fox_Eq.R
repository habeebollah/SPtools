function(data){
  CPUE <- log(data$catch/data$effort)
  ObsLR <- lm(CPUE ~ data$effort)
  mod.sum <- summary(ObsLR)
  rsquared <- c(mod.sum$r.squared, mod.sum$adj.r.squared)
  
  CI_a0.250 <- confint(ObsLR)[1] # lower CI for intercept
  CI_a0.975 <- confint(ObsLR)[3] # upper CI for intercept
  
  a <- mod.sum$coefficients[1] # intercept
  b <- mod.sum$coefficients[2] # slope
  
  MSY <- -(1/b)*exp(a-1)
  fMSY <- -1/b

  MSY_CI <- c(-(1/b)*exp(CI_a0.250-1), -(1/b)*exp(CI_a0.975-1))
  fMSY_CI <- c(-1/b, -1/b)

return(list(a=a, b=b, rsquared=rsquared, MSY=MSY, fMSY=fMSY, MSY_CI=MSY_CI, fMSY_CI=fMSY_CI))
}

#Schaefer_Eq(njava)
#Fox_Eq(njava)
