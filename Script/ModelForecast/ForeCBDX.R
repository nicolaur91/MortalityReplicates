#' CBD-X through STAN - both deaths and dynamic processes
#' 

# first long to wide

temp_df_D <- reshape(cbind.data.frame(list_men[[i]][,1:2],D=as.integer(list_men[[i]][,8])), idvar = "Age", timevar = "Year", direction = "wide")
temp_df_D <- data.frame(matrix(unlist(temp_df_D), nrow=length(temp_df_D), byrow=TRUE),stringsAsFactors=FALSE)
temp_df_E <- as.data.frame(reshape(list_men[[i]][,c(1,2,7)], idvar = "Age", timevar = "Year", direction = "wide"))
temp_df_E <- data.frame(matrix(unlist(temp_df_E), nrow=length(temp_df_E), byrow=TRUE),stringsAsFactors=FALSE)



# run STAN code
setwd("C:/Users/nsl.fi/Dropbox/PhD/Backup data/Mortality replicate/Bayesian CBD/Script/ModelForecast")

# males
temp_d <- rbind(as.matrix(list_men[[1]]$D),
                as.matrix(list_men[[2]]$D),
                as.matrix(list_men[[3]]$D),
                as.matrix(list_men[[4]]$D),
                as.matrix(list_men[[5]]$D))

temp_e <- rbind(as.matrix(list_men[[1]]$E),
                as.matrix(list_men[[2]]$E),
                as.matrix(list_men[[3]]$E),
                as.matrix(list_men[[4]]$E),
                as.matrix(list_men[[5]]$E))

temp_mort <- t(cbind(as.matrix(log(list_men[[1]]$D/list_men[[1]]$E))[which(list_men[[1]]$Year==1995)],
                   as.matrix(log(list_men[[2]]$D/list_men[[2]]$E))[which(list_men[[2]]$Year==1995)],
                   as.matrix(log(list_men[[3]]$D/list_men[[3]]$E))[which(list_men[[3]]$Year==1995)],
                   as.matrix(log(list_men[[4]]$D/list_men[[4]]$E))[which(list_men[[4]]$Year==1995)],
                   as.matrix(log(list_men[[5]]$D/list_men[[5]]$E))[which(list_men[[5]]$Year==1995)]))

# females
temp_d <- rbind(as.matrix(list_women[[1]]$D),
                as.matrix(list_women[[2]]$D),
                as.matrix(list_women[[3]]$D),
                as.matrix(list_women[[4]]$D),
                as.matrix(list_women[[5]]$D))

temp_e <- rbind(as.matrix(list_women[[1]]$E),
                as.matrix(list_women[[2]]$E),
                as.matrix(list_women[[3]]$E),
                as.matrix(list_women[[4]]$E),
                as.matrix(list_women[[5]]$E))

temp_mort <- t(cbind(as.matrix(log(list_women[[1]]$D/list_women[[1]]$E))[which(list_women[[1]]$Year==1995)],
                     as.matrix(log(list_women[[2]]$D/list_women[[2]]$E))[which(list_women[[2]]$Year==1995)],
                     as.matrix(log(list_women[[3]]$D/list_women[[3]]$E))[which(list_women[[3]]$Year==1995)],
                     as.matrix(log(list_women[[4]]$D/list_women[[4]]$E))[which(list_women[[4]]$Year==1995)],
                     as.matrix(log(list_women[[5]]$D/list_women[[5]]$E))[which(list_women[[5]]$Year==1995)]))


# # below only needed if need to rerun model
stan_data <- list(J = as.integer(iAge),
                  T = as.integer(iYear),
                  H = as.integer(iSubPop),
                  #H = as.integer(1),
                  Tfor = as.integer(30),
                  d = as.integer(as.vector(temp_d)),
                  e = as.integer(as.vector(temp_e)),
                  #k_bar = as.vector(seq(1,iYear)),
                  #k2_bar = as.vector(seq(1,iYear)),
                  # d = as.integer(as.vector(list_men[[1]]$sD)),
                  # e = as.integer(as.vector(list_men[[1]]$sE)),
                  #ones = as.vector(rep(1,length=iSubPop)),
                  mort = temp_mort,
                  age = unique(list_men[[1]]$Age),
                  family = 0)

stan_data <- list(J = as.integer(iAge),
                  T = as.integer(iYear),
                  H = 1,
                  #H = as.integer(1),
                  Tfor = as.integer(30),
                  d = as.integer(as.vector( as.matrix(list_men[[1]]$D)
                                            )),
                  e = as.integer(as.vector(as.matrix(list_men[[1]]$E)
                                                )),
                  #k_bar = as.vector(seq(1,iYear)),
                  #k2_bar = as.vector(seq(1,iYear)),
                  # d = as.integer(as.vector(list_men[[1]]$sD)),
                  # e = as.integer(as.vector(list_men[[1]]$sE)),
                  #ones = as.vector(rep(1,length=iSubPop)),
                  mort =as.matrix(log(list_men[[1]]$D/list_men[[1]]$E))[which(list_men[[1]]$Year==1995)],
                  age = unique(list_men[[1]]$Age),
                  family = 0)

stan_data <- list(J = as.integer(iAge),
                  T = as.integer(iYear),
                  #H = 1,
                  #H = as.integer(1),
                  Tfor = as.integer(30),
                  d = as.integer(as.vector( as.matrix(list_men[[1]]$D))),
                  e = as.integer(as.vector(as.matrix(list_men[[1]]$E))),
                  sqrt_d_inv = sqrt(1/as.vector( as.matrix(list_men[[1]]$D))),
                  #k_bar = as.vector(seq(1,iYear)),
                  #k2_bar = as.vector(seq(1,iYear)),
                  # d = as.integer(as.vector(list_men[[1]]$sD)),
                  # e = as.integer(as.vector(list_men[[1]]$sE)),
                  #ones = as.vector(rep(1,length=iSubPop)),
                  #mort =as.matrix(log(list_men[[1]]$D/list_men[[1]]$E))[which(list_men[[1]]$Year==1995)],
                  age = unique(list_men[[1]]$Age),
                  family = 0)

fitted.model <- stan(file = 'CBDX.stan',data=stan_data,chains=4,thin = 1,iter = 4000,warmup = 1500,cores=4,save_warmup = FALSE,
                     control = list( max_treedepth =15),
                     seed = 2020) # måske er beta centraliseret ved 0.5 det bedste.
fitted.model <- stan(file = 'CBDX_RP.stan',data=stan_data,chains=4,thin = 1,iter = 6000,warmup = 5000,cores=4,save_warmup = FALSE,
                     control = list( max_treedepth =15),
                     seed = 2020) # måske er beta centraliseret ved 0.5 det bedste.

# high autocorrelation try thinning - not neccessarily the right way to go.

#save(fitted.model,file = "CBX_jointtest_men.RData")
# load("CBX_jointtest_men.RData")

# launch_shinystan(CBD.BayesX_men)
# check fit
fit2 <- rstan::extract(fitted.model)
nsim <- dim(fit2$k)[1] 
fit_mort <- data.frame(m_1 = rep(NA,length=iYear*iSubPop))
# fitted mortality:
for (i in 1:iSubPop) {
  # fit_mort[((i-1)*iYear*iAge+1):(i*iYear*iAge),1] <- exp(as.vector(matrix(colMeans(fit2$bx)[i,],nrow=iAge,ncol=iYear)+
  #                                                            t(matrix(colMeans(fit2$k)[i,],nrow=iYear,ncol=iAge))+
  #                                                            (unique(list_men[[1]]$Age)-mean(unique(list_men[[1]]$Age)))%*%t(colMeans(fit2$k2)[i,])))
  
  fit_mort[((i-1)*iYear*iAge+1):(i*iYear*iAge),1] <- exp(as.vector(matrix(colMeans(fit2$bx)[i,],nrow=iAge,ncol=iYear)+
                                                                     t(matrix(colMeans(fit2$k)[i,],nrow=iYear,ncol=iAge))+
                                                                     (unique(list_men[[1]]$Age[which(list_men[[1]]$Age>=60)])-
                                                                        mean(unique(list_men[[1]]$Age[which(list_men[[1]]$Age>=60)])))%*%t(colMeans(fit2$k2)[i,])))
  
  # fit_mort[((i-1)*iYear*iAge+1):(i*iYear*iAge),1] <- exp(as.vector(matrix(t(temp_bx)[i,],nrow=iAge,ncol=iYear)+
  #                                                                    t(matrix(colMeans(fit2$k)[i,],nrow=iYear,ncol=iAge))+
  #                                                                    (unique(list_men[[1]]$Age)-mean(unique(list_men[[1]]$Age)))%*%t(colMeans(fit2$k2)[i,])))
}


# fit_mort[1:(iYear*iAge)] <- exp(as.vector(matrix(colMeans(fit2$bx)[,1],nrow=iAge,ncol=iYear)+
#                                         t(matrix(colMeans(fit2$k)[,1],nrow=iYear,ncol=iAge))+
#                                     (unique(list_men[[1]]$Age)-mean(unique(list_men[[1]]$Age)))%*%t(colMeans(fit2$k2)[,1])))

# # standardized residuals - from Cairns, MKL (2019), pp. 37
sr <- (temp_d/temp_e - fit_mort[,1])/sqrt(fit_mort[,1]/temp_e)
summary(sr) # ought to be approx. 0
var(sr) # ought to be around 1
# plot(fit_mort[,1],sr)
# abline(h=c(-2,0,2),lty=c(2,1,2))
# # estimated overdispersion
# 1/(length(sr)-550)*sum(sr^2)
# 1-pchisq(sum(sr^2),length(sr)-550) # if sum(sr^2) = length(sr)-550 then we have approximately equi-distributed. The null is equi-distributed.

posterior_fit <- as.array(fitted.model)

mcmc_scatter(posterior_fit,
             pars=c("k[1,2]","k_bar[2]"))

# with StanMoMo
fitCBD=cbd_stan(death = list_men[[1]]$D,exposure=list_men[[1]]$E, J=iAge,T=iYear,age=temp_df_D[,1], forecast = 1, family = "poisson")

####---- Make for kappa-terms ----####
stan_data <- list(T = iYear,
                  H = iSubPop,
                  P = 1,
                  Kappa_1 = t(cbind(as.numeric(coef(get(paste("CBX_grp_men",1,sep="_")))[(iAge+1):(iAge+iYear)]),
                                  as.numeric(coef(get(paste("CBX_grp_men",2,sep="_")))[(iAge+1):(iAge+iYear)]),
                                  as.numeric(coef(get(paste("CBX_grp_men",3,sep="_")))[(iAge+1):(iAge+iYear)]),
                                  as.numeric(coef(get(paste("CBX_grp_men",4,sep="_")))[(iAge+1):(iAge+iYear)]),
                                  as.numeric(coef(get(paste("CBX_grp_men",5,sep="_")))[(iAge+1):(iAge+iYear)]))),
                  Kappa_2 = t(cbind(as.numeric(coef(get(paste("CBX_grp_men",1,sep="_")))[(iAge+iYear+1):(iAge+2*iYear)]),
                                  as.numeric(coef(get(paste("CBX_grp_men",2,sep="_")))[(iAge+iYear+1):(iAge+2*iYear)]),
                                  as.numeric(coef(get(paste("CBX_grp_men",3,sep="_")))[(iAge+iYear+1):(iAge+2*iYear)]),
                                  as.numeric(coef(get(paste("CBX_grp_men",4,sep="_")))[(iAge+iYear+1):(iAge+2*iYear)]),
                                  as.numeric(coef(get(paste("CBX_grp_men",5,sep="_")))[(iAge+iYear+1):(iAge+2*iYear)]))),
                  Kappa_1_bar = rowMeans(matrix(c(as.numeric(coef(CBX_grp_men_1)[(iAge+1):(iAge+iYear)]),
                                                      as.numeric(coef(CBX_grp_men_2)[(iAge+1):(iAge+iYear)]),
                                                      as.numeric(coef(CBX_grp_men_3)[(iAge+1):(iAge+iYear)]),
                                                      as.numeric(coef(CBX_grp_men_4)[(iAge+1):(iAge+iYear)]),
                                                      as.numeric(coef(CBX_grp_men_5)[(iAge+1):(iAge+iYear)])),nrow=iYear,ncol=iSubPop)),
                  Kappa_2_bar = rowMeans(matrix(c(as.numeric(coef(CBX_grp_men_1)[(iAge+iYear+1):(iAge+2*iYear)]),
                                                      as.numeric(coef(CBX_grp_men_2)[(iAge+iYear+1):(iAge+2*iYear)]),
                                                      as.numeric(coef(CBX_grp_men_3)[(iAge+iYear+1):(iAge+2*iYear)]),
                                                      as.numeric(coef(CBX_grp_men_4)[(iAge+iYear+1):(iAge+2*iYear)]),
                                                      as.numeric(coef(CBX_grp_men_5)[(iAge+iYear+1):(iAge+2*iYear)])),nrow=iYear,ncol=iSubPop)), # last term is a TxM kappa
                  T_new = horizon+1)
# 
# fitted.model <- stan(file = 'Kappa_p_predict.stan',data=stan_data,chains=2,iter = 10000,control = list(adapt_delta = 0.9999,max_treedepth=15)) # måske er beta centraliseret ved 0.5 det bedste.
fitted.model <- stan(file = 'Kappa_p_predict_single.stan',data=stan_data,chains=4,cores=4,iter = 5000,control = list(adapt_delta = 0.9999,max_treedepth=15,stepsize=2)) # måske er beta centraliseret ved 0.5 det bedste.
# 10 min.

# set to CBDX
CBX_grp_men_1_5 <- fitted.model

# simulations
pens_group <- 1
ref_group <- 3

#number of simulations
nsim = dim(extr_pred_CBDX_men$Kappa_1_pred)[1]

nsim = 5000

#CBD-X simulations - pension-cohort
CBDX_sim_mortality_PH = list()
CBDX_sim_Pop = list()

 

for (i in 1:nsim) {
  
  kt_temp_PH_CBDX <- cbind(matrix(as.data.frame(extr_pred_CBDX_men[which(names(extr_pred_CBDX_men)=="Kappa_1_pred")])[i,,],nrow = 31,ncol=5)[1:max_horizon,pens_group],
                           matrix(as.data.frame(extr_pred_CBDX_men[which(names(extr_pred_CBDX_men)=="Kappa_2_pred")])[i,,],nrow = 31,ncol=5)[1:max_horizon,pens_group])
  
  # with group 3 as reference group
  kt_temp_Pop_CBDX <- cbind(matrix(as.data.frame(extr_pred_CBDX_men[which(names(extr_pred_CBDX_men)=="Kappa_1_pred")])[i,,],nrow = 31,ncol=5)[1:max_horizon,ref_group],
                           matrix(as.data.frame(extr_pred_CBDX_men[which(names(extr_pred_CBDX_men)=="Kappa_2_pred")])[i,,],nrow = 31,ncol=5)[1:max_horizon,ref_group])
  
  CBDX_sim_mortality_PH[[i]] = matrix(nrow = iAge,ncol = max_horizon)
  
  CBDX_sim_Pop[[i]] = matrix(nrow = iAge,ncol = max_horizon)
  
  # CBDX
  CBDX_sim_mortality_PH[[i]] = exp(
    # matrix(list_men[[pens_group]]$cbdx_bx[1:iAge],nrow=iAge,ncol=max_horizon)+
    matrix(as.numeric(coef(get(paste("CBX_grp_men",pens_group,sep="_")))[1:iAge]),nrow=iAge,ncol=max_horizon)+
      t(matrix(unlist(kt_temp_PH_CBDX[,1]),ncol=iAge,nrow=max_horizon))+
      (unique(list_men[[pens_group]]$Age)-mean(list_men[[pens_group]]$Age))%*%matrix(unlist(kt_temp_PH_CBDX[,2]),nrow=1,ncol=max_horizon)
  )
  
  # # with full pop as reference group
  # CBDX_sim_Pop[[i]] = exp(
  #                           matrix(dat_full_pop_men$cbdx_bx[1:iAge],nrow=iAge,ncol=max_horizon)+
  #                           t(matrix(unlist(kt_temp_Pop_CBDX[,1]),ncol=iAge,nrow=max_horizon))+
  #                           (l_age:u_age-mean(l_age:u_age))%*%matrix(unlist(kt_temp_Pop_CBDX[,2]),nrow=1,ncol=max_horizon)
  #   
  # )
  
  # CBDX - with group 3 as reference group  
  CBDX_sim_Pop[[i]] = exp(
    # matrix(list_men[[ref_group]]$cbdx_bx[1:iAge],nrow=iAge,ncol=max_horizon)+
    matrix(as.numeric(coef(get(paste("CBX_grp_men",ref_group,sep="_")))[1:iAge]),nrow=iAge,ncol=max_horizon)+
      t(matrix(unlist(kt_temp_Pop_CBDX[,1]),ncol=iAge,nrow=max_horizon))+
      (unique(list_men[[ref_group]]$Age)-mean(list_men[[ref_group]]$Age))%*%matrix(unlist(kt_temp_Pop_CBDX[,2]),nrow=1,ncol=max_horizon)
    
  )
  
  # transform to end-of-year mortality rates
  CBDX_sim_mortality_PH[[i]] <- CBDX_sim_mortality_PH[[i]]/(1+0.5*CBDX_sim_mortality_PH[[i]])
  CBDX_sim_Pop[[i]] <- CBDX_sim_Pop[[i]]/(1+0.5*CBDX_sim_Pop[[i]])
}


####---- Calculate hedge ----####

#calculate survival rates and expected survival rates
S1 = SLL = Q = Q1 = QLL = list()
EQ1 = ES1 = ESLL = EQLL = vector(length = max_horizon)

for (j in 1:nsim) {
  S1[[j]] = SLL[[j]] = QLL[[j]] = Q1[[j]] = vector(length = max_horizon)
  
  S1[[j]] = cumprod(1-diag(CBDX_sim_Pop[[j]][x:length(CBDX_sim_Pop[[j]][,1]),]))
  SLL[[j]] = cumprod(1-diag(CBDX_sim_mortality_PH[[j]][x:length(CBDX_sim_mortality_PH[[j]][,1]),]))
  
  for (i in 1:max_horizon) {
    
    Q1[[j]][i]  = CBDX_sim_Pop[[j]][x+i-1,i]
    QLL[[j]][i] = CBDX_sim_mortality_PH[[j]][x+i-1,i]
    ES1[i] = ES1[i] + S1[[j]][i] 
    ESLL[i] = ESLL[i] + SLL[[j]][i] 
    EQ1[i] = EQ1[i] + Q1[[j]][i]
    EQLL[i] = EQLL[i] + QLL[[j]][i]
  }
}

# with basis risk
summary(lm(((1+r)^-(1:max_horizon) %*% t(do.call(rbind,SLL)))[1,] ~ as.matrix((1+r)^-5*do.call(rbind,Q1)[,5])[,1]+
                as.matrix((1+r)^-10*do.call(rbind,Q1)[,10])[,1]+as.matrix((1+r)^-15*do.call(rbind,Q1)[,15])[,1]+
                as.matrix((1+r)^-20*do.call(rbind,Q1)[,20])[,1]+as.matrix((1+r)^-25*do.call(rbind,Q1)[,25])[,1]))

# without basis risk
summary(lm(((1+r)^-(1:max_horizon) %*% t(do.call(rbind,SLL)))[1,] ~ as.matrix((1+r)^-5*do.call(rbind,QLL)[,5])[,1]+
                as.matrix((1+r)^-10*do.call(rbind,QLL)[,10])[,1]+as.matrix((1+r)^-15*do.call(rbind,QLL)[,15])[,1]+
                as.matrix((1+r)^-20*do.call(rbind,QLL)[,20])[,1]+as.matrix((1+r)^-25*do.call(rbind,QLL)[,25])[,1]))


####---- With joint ----####
# men
# set to CBDX
CBX_grp_men_1_5 <- fitted.model

# mufor has vectored size: iSubPop*30*iAge - starting with sub-population, next time period and last age. Iterates through age for iSubPop=1, Tfor=1, etc.

extr_pred_CBDX_men <- extract(CBX_grp_men_1_5)

# number of simulations
nsim <- dim(extr_pred_CBDX_men$k)[1]

# simulations
pens_group <- 5
ref_group <- 3

#CBD-X simulations - pension-cohort
CBDX_sim_mortality_PH = list()
CBDX_sim_Pop = list()

for (i in 1:nsim) {
  CBDX_sim_mortality_PH[[i]] = matrix(nrow = iAge,ncol = max_horizon)
  
  CBDX_sim_Pop[[i]] = matrix(nrow = iAge,ncol = max_horizon)
  
  CBDX_sim_mortality_PH[[i]] <- matrix((extr_pred_CBDX_men$mufor[i,])[((pens_group-1)*horizon*iAge+1):((pens_group-1)*horizon*iAge+iAge*horizon)],nrow=iAge,ncol=30)
  
  CBDX_sim_Pop[[i]] <- matrix((extr_pred_CBDX_men$mufor[i,])[((ref_group-1)*horizon*iAge+1):((ref_group-1)*horizon*iAge+iAge*horizon)],nrow=iAge,ncol=30)
  
  # end-of-year mortality rates
  CBDX_sim_mortality_PH[[i]] <- CBDX_sim_mortality_PH[[i]]/(1+0.5*CBDX_sim_mortality_PH[[i]])
  CBDX_sim_Pop[[i]] <- CBDX_sim_Pop[[i]]/(1+0.5*CBDX_sim_Pop[[i]])
  
}

####---- Calculate hedge ----####

#calculate survival rates and expected survival rates
S1 = SLL = Q = Q1 = QLL = list()
EQ1 = ES1 = ESLL = EQLL = vector(length = max_horizon)

for (j in 1:nsim) {
  S1[[j]] = SLL[[j]] = QLL[[j]] = Q1[[j]] = vector(length = max_horizon)
  
  S1[[j]] = cumprod(1-diag(CBDX_sim_Pop[[j]][x:length(CBDX_sim_Pop[[j]][,1]),]))
  SLL[[j]] = cumprod(1-diag(CBDX_sim_mortality_PH[[j]][x:length(CBDX_sim_mortality_PH[[j]][,1]),]))
  
  for (i in 1:max_horizon) {
    
    Q1[[j]][i]  = CBDX_sim_Pop[[j]][x+i-1,i]
    QLL[[j]][i] = CBDX_sim_mortality_PH[[j]][x+i-1,i]
    ES1[i] = ES1[i] + S1[[j]][i] 
    ESLL[i] = ESLL[i] + SLL[[j]][i] 
    EQ1[i] = EQ1[i] + Q1[[j]][i]
    EQLL[i] = EQLL[i] + QLL[[j]][i]
  }
}
# survival and mortality matrices from simulations
QLL_men_list[[pens_group]] <- do.call(rbind,QLL)
SLL_men_list[[pens_group]] <- do.call(rbind,SLL)

# with basis risk
summary(lm(((1+r)^-(1:max_horizon) %*% t(do.call(rbind,SLL)))[1,] ~ as.matrix((1+r)^-5*do.call(rbind,Q1)[,5])[,1]+
             as.matrix((1+r)^-10*do.call(rbind,Q1)[,10])[,1]+as.matrix((1+r)^-15*do.call(rbind,Q1)[,15])[,1]+
             as.matrix((1+r)^-20*do.call(rbind,Q1)[,20])[,1]+as.matrix((1+r)^-25*do.call(rbind,Q1)[,25])[,1]))

# without basis risk
summary(lm(((1+r)^-(1:max_horizon) %*% t(do.call(rbind,SLL)))[1,] ~ as.matrix((1+r)^-5*do.call(rbind,QLL)[,5])[,1]+
             as.matrix((1+r)^-10*do.call(rbind,QLL)[,10])[,1]+as.matrix((1+r)^-15*do.call(rbind,QLL)[,15])[,1]+
             as.matrix((1+r)^-20*do.call(rbind,QLL)[,20])[,1]+as.matrix((1+r)^-25*do.call(rbind,QLL)[,25])[,1]))


sim_mat_k <- matrix(nrow=nsim,ncol=max_horizon)
sim_mat_Pop <- matrix(nrow=nsim,ncol=max_horizon)

# premiums
for (i in 1:nsim) {
  
  sim_mat_k[i,] <- diag((CBDX_sim_mortality_PH[[i]])[x:iAge,1:max_horizon])
  sim_mat_Pop[i,] <- diag((CBDX_sim_Pop[[i]])[x:iAge,1:max_horizon])
  
}

prem_men[pens_group,] <- (exp(colMeans(log(sim_mat_k))+RP*sqrt(1:max_horizon)*colSds(log(sim_mat_k))+
                                0.5*colVars(log(sim_mat_k)))/colMeans(sim_mat_k)-1)*100

# absolute premium
prem_men_abs[pens_group,] <- RP*(1:max_horizon)*colVars(log(QLL_men_list[[pens_group]]))

####---- women ----####
# set to CBDX
CBX_grp_women_1_5 <- fitted.model

# mufor has vectored size: iSubPop*30*iAge - starting with sub-population, next time period and last age. Iterates through age for iSubPop=1, Tfor=1, etc.

extr_pred_CBDX_women <- extract(CBX_grp_women_1_5)

# simulations
pens_group <- 5
ref_group <- 3

#CBD-X simulations - pension-cohort
CBDX_sim_mortality_PH = list()
CBDX_sim_Pop = list()

for (i in 1:nsim) {
  CBDX_sim_mortality_PH[[i]] = matrix(nrow = iAge,ncol = max_horizon)
  
  CBDX_sim_Pop[[i]] = matrix(nrow = iAge,ncol = max_horizon)
  
  CBDX_sim_mortality_PH[[i]] <- matrix((extr_pred_CBDX_women$mufor[i,])[((pens_group-1)*horizon*iAge+1):((pens_group-1)*horizon*iAge+iAge*horizon)],nrow=iAge,ncol=30)
  
  CBDX_sim_Pop[[i]] <- matrix((extr_pred_CBDX_women$mufor[i,])[((ref_group-1)*horizon*iAge+1):((ref_group-1)*horizon*iAge+iAge*horizon)],nrow=iAge,ncol=30)
  
  # end-of-year mortality rates
  CBDX_sim_mortality_PH[[i]] <- CBDX_sim_mortality_PH[[i]]/(1+0.5*CBDX_sim_mortality_PH[[i]])
  CBDX_sim_Pop[[i]] <- CBDX_sim_Pop[[i]]/(1+0.5*CBDX_sim_Pop[[i]])
  
}

####---- Calculate hedge ----####

#calculate survival rates and expected survival rates
S1 = SLL = Q = Q1 = QLL = list()
EQ1 = ES1 = ESLL = EQLL = vector(length = max_horizon)

for (j in 1:nsim) {
  S1[[j]] = SLL[[j]] = QLL[[j]] = Q1[[j]] = vector(length = max_horizon)
  
  S1[[j]] = cumprod(1-diag(CBDX_sim_Pop[[j]][x:length(CBDX_sim_Pop[[j]][,1]),]))
  SLL[[j]] = cumprod(1-diag(CBDX_sim_mortality_PH[[j]][x:length(CBDX_sim_mortality_PH[[j]][,1]),]))
  
  for (i in 1:max_horizon) {
    
    Q1[[j]][i]  = CBDX_sim_Pop[[j]][x+i-1,i]
    QLL[[j]][i] = CBDX_sim_mortality_PH[[j]][x+i-1,i]
    ES1[i] = ES1[i] + S1[[j]][i] 
    ESLL[i] = ESLL[i] + SLL[[j]][i] 
    EQ1[i] = EQ1[i] + Q1[[j]][i]
    EQLL[i] = EQLL[i] + QLL[[j]][i]
  }
}

# with basis risk
summary(lm(((1+r)^-(1:max_horizon) %*% t(do.call(rbind,SLL)))[1,] ~ as.matrix((1+r)^-5*do.call(rbind,Q1)[,5])[,1]+
             as.matrix((1+r)^-10*do.call(rbind,Q1)[,10])[,1]+as.matrix((1+r)^-15*do.call(rbind,Q1)[,15])[,1]+
             as.matrix((1+r)^-20*do.call(rbind,Q1)[,20])[,1]+as.matrix((1+r)^-25*do.call(rbind,Q1)[,25])[,1]))

# without basis risk
summary(lm(((1+r)^-(1:max_horizon) %*% t(do.call(rbind,SLL)))[1,] ~ as.matrix((1+r)^-5*do.call(rbind,QLL)[,5])[,1]+
             as.matrix((1+r)^-10*do.call(rbind,QLL)[,10])[,1]+as.matrix((1+r)^-15*do.call(rbind,QLL)[,15])[,1]+
             as.matrix((1+r)^-20*do.call(rbind,QLL)[,20])[,1]+as.matrix((1+r)^-25*do.call(rbind,QLL)[,25])[,1]))


