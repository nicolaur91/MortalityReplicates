#' Bayesian forecast and simulation
#' If I use Kappa_1_new then I get 99% R^2, but the hedge weights are of.. 

####---- Kappa-processes ----####

# # below only needed if need to rerun model
stan_data <- list(T = iYear,
                  Kappa_1 = as.numeric(dat_EW_mort$k1_CBD[which(dat_EW_mort$Age==l_age)]),
                  T_new = horizon+1)
# 
fitted.model <- stan(file = 'Kappa_1_predict.stan',data=stan_data,chains=4,iter = 10000,control = list(adapt_delta = 0.9999,max_treedepth=15)) # måske er beta centraliseret ved 0.5 det bedste.

summary(fitted.model)$summary

# for kappa 2 process
stan_data <- list(T = iYear,
                  Kappa_1 = as.numeric(dat_EW_mort$k2_CBD[which(dat_EW_mort$Age==l_age)]),
                  T_new = horizon+1)
# 
fitted.model_2 <- stan(file = 'Kappa_1_predict.stan',data=stan_data,chains=4,iter = 10000,control = list(adapt_delta = 0.9999,max_treedepth=15)) # måske er beta centraliseret ved 0.5 det bedste.

summary(fitted.model_2)$summary

####---- Put in simulations in expected ----####
# full population list men
CBD.bayes.forecast_mortality = list(matrix(nrow = length(1:iAge),ncol = horizon[1]))
CBD.bayes.forecast_mortality_UPPER = list(matrix(nrow = length(1:iAge),ncol = horizon[1]))
CBD.bayes.forecast_mortality_LOWER = list(matrix(nrow = length(1:iAge),ncol = horizon[1]))

# with kappa_new
# kt_temp_CBD <- matrix(cbind(
#   (summary(fitted.model)$summary[which(substr(row.names(as.data.frame(summary(fitted.model)$summary)),
#                                                  1,11)=="Kappa_1_new"),
#                                     1])[2:(horizon+1)],
#   (summary(fitted.model_2)$summary[which(substr(row.names(as.data.frame(summary(fitted.model_2)$summary)),
#                                                  1,11)=="Kappa_1_new"),
#                                     1])[2:(horizon+1)]),
#   nrow=horizon,ncol=2)

# with kappa_pred
kt_temp_CBD <- matrix(cbind(
  (summary(fitted.model)$summary[which(substr(row.names(as.data.frame(summary(fitted.model)$summary)),
                                              1,12)=="Kappa_1_pred"),
                                 1])[1:horizon],
  (summary(fitted.model_2)$summary[which(substr(row.names(as.data.frame(summary(fitted.model_2)$summary)),
                                                1,12)=="Kappa_1_pred"),
                                   1])[1:horizon]),
  nrow=horizon,ncol=2) 

CBD.bayes.forecast_mortality <- t(matrix(kt_temp_CBD[,1],ncol=iAge,nrow=horizon[1]))+
  (l_age:u_age-mean(l_age:u_age))%*%t(as.matrix(kt_temp_CBD[,2]))

# lower 2.5%
kt_temp_CBD <- matrix(cbind(
  (summary(fitted.model)$summary[which(substr(row.names(as.data.frame(summary(fitted.model)$summary)),
                                              1,12)=="Kappa_1_pred"),
                                 4])[2:(horizon+1)],
  (summary(fitted.model_2)$summary[which(substr(row.names(as.data.frame(summary(fitted.model_2)$summary)),
                                                1,12)=="Kappa_1_pred"),
                                   4])[2:(horizon+1)]),
  nrow=horizon,ncol=2) 

CBD.bayes.forecast_mortality_UPPER <- t(matrix(kt_temp_CBD[,1],ncol=iAge,nrow=horizon[1]))+
                            (l_age:u_age-mean(l_age:u_age))%*%t(as.matrix(kt_temp_CBD[,2]))


# upper 97.5%
kt_temp_CBD <- matrix(cbind(
  (summary(fitted.model)$summary[which(substr(row.names(as.data.frame(summary(fitted.model)$summary)),
                                              1,12)=="Kappa_1_pred"),
                                 8])[2:(horizon+1)],
  (summary(fitted.model_2)$summary[which(substr(row.names(as.data.frame(summary(fitted.model_2)$summary)),
                                                1,12)=="Kappa_1_pred"),
                                  8])[2:(horizon+1)]),
  nrow=horizon,ncol=2) 

CBD.bayes.forecast_mortality_LOWER <- t(matrix(kt_temp_CBD[,1],ncol=iAge,nrow=horizon[1]))+
  (l_age:u_age-mean(l_age:u_age))%*%t(as.matrix(kt_temp_CBD[,2]))

####---- Hedge effectiveness ----####

# extract simulations from Stan-object
extr_pred_CBD_1 <- rstan::extract(fitted.model)
extr_pred_CBD_2 <- rstan::extract(fitted.model_2)

CBD_bayes_sim_mortality = list()

for (i in 1:nsim) {
  
  kt_temp_CBD <- matrix(cbind(as.data.frame(extr_pred_CBD_1[which(names(extr_pred_CBD_1)=="Kappa_1_pred")])[i,,],
                        as.data.frame(extr_pred_CBD_2[which(names(extr_pred_CBD_2)=="Kappa_1_pred")])[i,,]
                        ),nrow=horizon,ncol=2)[1:max_horizon,]
  
  # # with full population as reference group
  # kt_temp_Pop_CBDX <- matrix(as.data.frame(extr_pred_CBDX_men[which(names(extr_pred_CBDX_men)=="Kappa_full_pred")])[i,,],
  #                            nrow=30,ncol=2)[1:max_horizon,]
  
  CBD_bayes_sim_mortality[[i]] = matrix(nrow = iAge,ncol = max_horizon)
  
  CBD_bayes_sim_mortality[[i]] = exp(
      t(matrix(unlist(kt_temp_CBD[,1]),ncol=iAge,nrow=max_horizon))+
      (unique(dat_EW_mort$Age)-mean(dat_EW_mort$Age))%*%matrix(unlist(kt_temp_CBD[,2]),nrow=1,ncol=max_horizon)
  ) 
  # # with full pop as reference group
  # CBDX_sim_Pop[[i]] = exp(
  #                           matrix(dat_full_pop_men$cbdx_bx[1:iAge],nrow=iAge,ncol=max_horizon)+
  #                           t(matrix(unlist(kt_temp_Pop_CBDX[,1]),ncol=iAge,nrow=max_horizon))+
  #                           (l_age:u_age-mean(l_age:u_age))%*%matrix(unlist(kt_temp_Pop_CBDX[,2]),nrow=1,ncol=max_horizon)
  #   
  # )
  
  
  # transform to end-of-year mortality rates
  CBD_bayes_sim_mortality[[i]] <- CBD_bayes_sim_mortality[[i]]/(1+0.5*CBD_bayes_sim_mortality[[i]])

}

# simulated kappa
sim_mat_bayes_k <- matrix(nrow=nsim,ncol=max_horizon)

for (i in 1:nsim) {
  
  sim_mat_bayes_k[i,] <- diag((CBD_bayes_sim_mortality[[i]])[x:iAge,1:max_horizon])
  
}

####---- Calculate hedge ----####

#calculate survival rates and expected survival rates
SLL_bayes = QLL_bayes = list()


for (j in 1:nsim) {
  SLL_bayes[[j]] = QLL_bayes[[j]] = vector(length = max_horizon)
  
  # S1[[j]] = cumprod(1-diag(lcsim_Pop[[j]][x:length(lcsim_Pop[[j]][,1]),]))
  SLL_bayes[[j]] = cumprod(1-diag(CBD_bayes_sim_mortality[[j]][x:length(CBD_bayes_sim_mortality[[j]][,1]),]))
  
  for (i in 1:max_horizon) {
    
    # Q1[[j]][i]  = lcsim_Pop[[j]][x+i-1,i]
    QLL_bayes[[j]][i] = CBD_bayes_sim_mortality[[j]][x+i-1,i]
  
  }
}

# R^2 94.86 med my normal(0,1), sigma ~ cauchy(0,1)% - R^2 94.86 med my normal(0,1), sigma ~ inv_gamma(11,0.1) 
# R^2 93.9 med my normal(0,2), sigma ~ inv_gamma(11,0.1)
# R^2 94.14 med my normal(0,10), sigma ~ inv_gamma(11,0.1)
# R^2 95 med informal prior.
summary(lm(((1+r)^-(1:max_horizon) %*% t(do.call(rbind,SLL_bayes)))[1,] ~ as.matrix((1+r)^-5*do.call(rbind,QLL_bayes)[,5])[,1]+
             as.matrix((1+r)^-10*do.call(rbind,QLL_bayes)[,10])[,1]+as.matrix((1+r)^-15*do.call(rbind,QLL_bayes)[,15])[,1]+
             as.matrix((1+r)^-20*do.call(rbind,QLL_bayes)[,20])[,1]+as.matrix((1+r)^-25*do.call(rbind,QLL_bayes)[,25])[,1]))
