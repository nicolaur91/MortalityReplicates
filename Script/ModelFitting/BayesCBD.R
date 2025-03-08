#' Fit the Bayesian version of the CBD-model

CBD.bayes <- stan_glm(formula=as.formula(dat_EW_mort$mD ~ -1 + X),
                family = poisson,
                offset = log(dat_EW_mort$mE),iter=10000, chains = 2, sparse=TRUE,seed = 1991)

summary(CBD.bayes)
iYear*2

print(CBD.bayes)
covmat(CBD.bayes)
CBD.bayes$covmat

# Sampling diagnostics
# effective sample sizes
neff_ratio(CBD.bayes)

# rhat
rhat(CBD.bayes)
mcmc_nuts_divergence(nuts_params(CBD.bayes), log_posterior(CBD.bayes))

mcmc_parcoord(as.array(CBD.bayes),nuts_params(CBD.bayes))

# trace plot of different variables
mcmc_trace( CBD.bayes, pars = "Xk1_1965",np=nuts_params(CBD.bayes))

# posterior predictive checks
yrep_CBD <- posterior_predict(CBD.bayes, draws = 500)

ppc_hist(dat_EW_mort$mD, yrep_CBD[1:5, ])
ppc_dens_overlay(dat_EW_mort$mD, yrep_CBD[1:50, ])

pp_check(CBD.bayes)

# test statistics - number of zeros
mean(dat_EW_mort$mD==0)
mean(yrep_CBD==0)
summary(rowMeans(yrep_CBD==0)) # to many zeros. Why is this so?


# Fitted values
exp(matrix(coef(CBD.bayes)[1:iYear],nrow=iAge,ncol=iYear)+as.matrix(unique(dat_EW_mort$Age)-mean(l_age:u_age)) %*% t(coef(CBD.bayes)[(iYear+1):(2*iYear)]))

# predicted interval
pred_inter_CBD <- predictive_interval(CBD.bayes,prob=0.95,pars="k2_1997")


####---- Comparing ML with Bayes ----####
round(coef(summary(CBD.male)), 3)
cbind(Median = coef(CBD.bayes), MAD_SD = se(CBD.bayes))

# conclusion: ML-estimates have low and zero SD, while Bayesian have just low SD.


####---- insert bayes-CBD values into data set ----####
dat_EW_mort$k1_B_CBD <- rep(coef(CBD.bayes)[1:length(unique(dat_EW_mort$Year))],each=length(unique(dat_EW_mort$Age)))
dat_EW_mort$k2_B_CBD <- rep(coef(CBD.bayes)[(length(unique(dat_EW_mort$Year))+1):length(coef(CBD.bayes))],
                            each=length(unique(dat_EW_mort$Age)))