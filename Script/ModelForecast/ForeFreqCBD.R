#' Forecast CBD-model frequentist style
#' Just assume that both kappa-terms are running a random walk-process
#' 

horizon <- 30

# k1 freq fit
k1_freq_lm <- lm(diff(dat_EW_mort$k1_CBD[which(dat_EW_mort$Age==l_age)]) ~ 1)
# k2 freq fit
k2_freq_lm <- lm(diff(dat_EW_mort$k2_CBD[which(dat_EW_mort$Age==l_age)]) ~ 1)

# forecast k1
forecast_k1_freq <- data.frame(E_k1 = dat_EW_mort$k1_CBD[nrow(dat_EW_mort)]+(1:horizon)*coef(k1_freq_lm))
forecast_k1_freq$l_k1 <- forecast_k1_freq$E_k1 + sqrt(1:horizon)*sd(resid(k1_freq_lm))*qnorm(0.025)
forecast_k1_freq$u_k1 <- forecast_k1_freq$E_k1 + sqrt(1:horizon)*sd(resid(k1_freq_lm))*qnorm(0.975)

# forecast k2
forecast_k2_freq <- data.frame(E_k2 = dat_EW_mort$k2_CBD[nrow(dat_EW_mort)]+(1:horizon)*coef(k2_freq_lm))
forecast_k2_freq$l_k2 <- forecast_k2_freq$E_k2 + sqrt(1:horizon)*sd(resid(k2_freq_lm))*qnorm(0.025)
forecast_k2_freq$u_k2 <- forecast_k2_freq$E_k2 + sqrt(1:horizon)*sd(resid(k2_freq_lm))*qnorm(0.975)


####---- Insert into list ----####

# list with projections
CBD.forecast_mortality = list(matrix(nrow = length(1:iAge),ncol = horizon[1]))
CBD.forecast_mortality_UPPER = list(matrix(nrow = length(1:iAge),ncol = horizon[1]))
CBD.forecast_mortality_LOWER = list(matrix(nrow = length(1:iAge),ncol = horizon[1]))

# expected mortality
CBD.forecast_mortality <- 
  t(matrix(forecast_k1_freq$E_k1,ncol=iAge,nrow=horizon[1]))+
  (unique(dat_EW_mort$Age)-mean(dat_EW_mort$Age))%*%t(as.matrix(forecast_k2_freq$E_k2))

CBD.forecast_mortality_UPPER <- 
  t(matrix(forecast_k1_freq$u_k1,ncol=iAge,nrow=horizon[1]))+
  (unique(dat_EW_mort$Age)-mean(dat_EW_mort$Age))%*%t(as.matrix(forecast_k2_freq$u_k2))

CBD.forecast_mortality_LOWER <- 
  t(matrix(forecast_k1_freq$l_k1,ncol=iAge,nrow=horizon[1]))+
  (unique(dat_EW_mort$Age)-mean(dat_EW_mort$Age))%*%t(as.matrix(forecast_k2_freq$l_k2))

# plot(diag(exp(CBD.forecast_mortality)),type="b")
# lines(diag(exp(CBD.forecast_mortality_UPPER)))
# lines(diag(exp(CBD.forecast_mortality_LOWER)))