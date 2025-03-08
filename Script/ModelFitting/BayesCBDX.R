#' Script til at fitte en Bayesian udgave af CBD-X modellen, som de bruger i "Affluence index"-artiklen
#' Skal tillægge iAge - kolonner i Kronecker produktet.
#' 


CBD.bayesX <- stan_glm(formula=as.formula(dat_EW_mort$mD ~ -1 + XX),
                      family = poisson,
                      offset = log(dat_EW_mort$mE), iter=2000, chains = 2,sparse=TRUE,seed = 1991) # har maximum 4 cores

# setting priors to ML-estimates
CBD.bayesX <- stan_glm(formula=as.formula(dat_EW_mort$mD ~ -1 + XX),
                       family = poisson,
                       offset = log(dat_EW_mort$mE), iter=1000, chains = 2,sparse=TRUE,seed = 1991,
                       prior = normal(location=coef(CBD.freqX),scale=c(rep(0.001,iAge),rep(1,2*iYear))),
                       control = list(adapt_delta=0.9999)) # har maximum 4 cores

CBD.bayesX <- stan_glm(formula=as.formula(as.integer(dat_full_pop_men$sD) ~ -1 + CBDX),
                       family = poisson,
                       offset = log(dat_full_pop_men$sE), iter=1000, chains = 2,sparse=TRUE,seed = 1991,
                       prior = normal(location=c(log(dat_full_pop_men$sD/dat_full_pop_men$sE)[1:iAge],
                                                 seq(0,-0.5,length.out = iYear),
                                                 seq(0,0.015,length.out = iYear)),scale=c(rep(0.05,iAge),rep(1,2*iYear))),
                       control = list(adapt_delta=0.9999)) # har maximum 4 cores

pairs(CBD.bayesX,pars=c("XXBx_62","XXBx_63"))

summary(CBD.bayesX)

prior_summary(CBD.bayesX)



# shiny stan
launch_shinystan(CBD.bayesX)

# sampling diagnostics
available_mcmc(pattern = "_nuts_")
log_posterior(CBD.bayesX)
nuts_params(CBD.bayesX)

mcmc_parcoord(as.array(CBD.bayesX),nuts_params(CBD.bayesX))
parameters(CBD.bayesX)

mcmc_nuts_divergence(nuts_params(CBD.bayesX), log_posterior(CBD.bayesX))

# trace of plot of 
mcmc_trace( CBD.bayesX, pars = "CBDXk2_2014",np=nuts_params(CBD.bayesX))

mcmc_nuts_energy(nuts_params(CBD.bayesX))

# effective sample sizes
neff_ratio(CBD.bayesX)

# rhat
rhat(CBD.bayesX)

pp_check(CBD.bayesX)

mcmc_rhat(rhat(CBD.bayesX)[(iAge+iYear+1):(iAge+2*iYear)])+ yaxis_text(hjust = 1)
mcmc_rhat(rhat(CBD.bayesX)[1:iAge])+ yaxis_text(hjust = 1)

####---- Comparing ML with Bayes ----####
round(coef(summary(CBD.freqX)), 3)
cbind(Median = coef(CBD.bayesX), MAD_SD = se(CBD.bayesX))

predictive_interval(CBD.bayesX,prob=0.95,pars="XXk2_2008")

# k1
plot(l_time:u_time,coef(CBD.bayesX)[(iAge+1):(iAge+iYear)],type="l")
plot(l_time:u_time,coef(CBD.freqX)[(iAge+1):(iAge+iYear)],type="l")


plot(l_time:u_time,coef(CBD.freqX_men)[(iAge+1):(iAge+iYear)],type="l")

# k2
plot(l_time:u_time,coef(CBD.bayesX)[(iAge+iYear+1):(iAge+2*iYear)],type="l")
plot(l_time:u_time,coef(CBD.freqX)[(iAge+iYear+1):(iAge+2*iYear)],type="l")

####---- Posterior predictive distribution ----####
# use method as used in Modelling Mortality with Actuarial Applications by MacDonald
X1_new = kronecker(diag(iYear+30),rep(1,iAge))
X2_new = kronecker(diag(iYear+30),unique(list_men[[1]]$Age)-mean(l_age:u_age))

  # Column names
colnames(X1_new) <- paste(rep("k1_",iYear+30),l_time:(u_time+30),sep="")
colnames(X2_new) <- paste(rep("k2_",iYear+30),l_time:(u_time+30),sep="")

X_new = cbind(X1_new,X2_new)

# Add age to coloumns
X_age_new <-  kronecker(rep(1,iYear+30),diag(iAge))

colnames(X_age_new) <- paste(rep("Bx_",iAge),l_age:(u_age),sep="")

CBDX_new <- cbind(X_age_new,X_new)

post_pred_new <- posterior_predict(CBD.bayesX, newdata = as.data.frame(CBDX_new))
