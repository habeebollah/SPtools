function(data){
  CPUE <- data$catch/data$effort
  ObsLR <- lm(CPUE ~ data$effort)
  mod.sum <- summary(ObsLR)
  rsquared <- c(mod.sum$r.squared, mod.sum$adj.r.squared)
  
  CI_a0.250 <- confint(ObsLR)[1] # lower CI for intercept
  CI_a0.975 <- confint(ObsLR)[3] # upper CI for intercept
  
  a <- mod.sum$coefficients[1] # intercept
  b <- mod.sum$coefficients[2] # slope
  
  MSY <- -0.25 * a^2/b
  fMSY <- -0.5 * a/b
  
  MSY_CI <- c(-0.25 * CI_a0.250^2/b, -0.25 * CI_a0.975^2/b)
  fMSY_CI <- c(-0.5 * CI_a0.250/b, -0.5 * CI_a0.975/b)
  
  res <- list()
  #return(list(a=a, b=b, rsquared=rsquared, MSY=MSY, fMSY=fMSY, MSY_CI=MSY_CI, fMSY_CI=fMSY_CI))
}

#prod_mod(data = trawl_fishery_Java, plot = TRUE)
