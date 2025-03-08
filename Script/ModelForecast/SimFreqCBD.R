#' Simulating the CBD-model 
#' 

require(matrixStats)
require(tsDyn)
require(expm)

####---- Initial values ----####
#Initial parameter values
t0 = u_time # Time 0
n = 5 # number of q-forward contracts
CurrentAge = l_age # Age at Time 1 at pension
x_u = u_age #upper age limit
ye_betw = round(length(CurrentAge:x_u)/n)-1 # years between each equi-distant contract
Ages = (CurrentAge:x_u)[seq(1:n)*ye_betw]  # equi-distant ages
n = length(Ages)
r = 0.02  #interest rate
RP = 0.2 # yearly risk premium for buying q-forwards
delta = 0.001 #small change in mortality rate used for calculation of key q-duration

Tj = vector(length = n) #vector containing relevant years
for (i in 1:n) {
  Tj[i] = t0 + Ages[i]-CurrentAge
}
x = CurrentAge - min(dat_EW_mort$Age) +1

# max horizon
max_horizon = min(horizon[1],x_u-CurrentAge+1)

#number of simulations
nsim = 10000 

####---- Simulate ----####

####---- Hedge effectiveness ----####
set.seed(1991)
#Need to generate mortality scenarios in order to calculate hedge effectiveness.

#Unexpected cashflows
#Without hedge:   X = L(q) - L(E(q))
#With hedge:      X* = X - H(q) + H(E(q))

#Li-Lee simulations - pension-cohort
cbd_sim_mortality = list()

for (i in 1:nsim) {
  
  sim_temp_k1 <- forecast_k1_freq$E_k1[1:max_horizon]+cumsum(rnorm(n=max_horizon,mean = 0,
                                sd=sd(resid(k1_freq_lm))))
  sim_temp_k2 <- forecast_k2_freq$E_k2[1:max_horizon]+cumsum(rnorm(n=max_horizon,mean = 0,
                                sd=sd(resid(k2_freq_lm))))
  
  cbd_sim_mortality[[i]] = matrix(nrow = iAge,ncol = max_horizon)
  
  cbd_sim_mortality[[i]] = exp(t(matrix(sim_temp_k1,ncol=iAge,nrow=max_horizon))+
                                  (unique(dat_EW_mort$Age)-mean(dat_EW_mort$Age))%*%t(as.matrix(sim_temp_k2))
                                
  )
  
  # transform to end-of-year deaths
  cbd_sim_mortality[[i]] <- cbd_sim_mortality[[i]]/(1+0.5*cbd_sim_mortality[[i]])
}

# simulated kappa
sim_mat_k <- matrix(nrow=nsim,ncol=max_horizon)

for (i in 1:nsim) {
  
  sim_mat_k[i,] <- diag((cbd_sim_mortality[[i]])[x:iAge,1:max_horizon])

}

####---- Calculate hedge ----####

#calculate survival rates and expected survival rates
S1 = SLL = Q = Q1 = QLL = list()
EQ1 = ES1 = ESLL = EQLL = vector(length = max_horizon)

for (j in 1:nsim) {
  S1[[j]] = SLL[[j]] = QLL[[j]] = Q1[[j]] = vector(length = max_horizon)
  
  # S1[[j]] = cumprod(1-diag(lcsim_Pop[[j]][x:length(lcsim_Pop[[j]][,1]),]))
  SLL[[j]] = cumprod(1-diag(cbd_sim_mortality[[j]][x:length(cbd_sim_mortality[[j]][,1]),]))
  
  for (i in 1:max_horizon) {
    
    # Q1[[j]][i]  = lcsim_Pop[[j]][x+i-1,i]
    QLL[[j]][i] = cbd_sim_mortality[[j]][x+i-1,i]
    # ES1[i] = ES1[i] + S1[[j]][i] 
    ESLL[i] = ESLL[i] + SLL[[j]][i] 
    # EQ1[i] = EQ1[i] + Q1[[j]][i]
    EQLL[i] = EQLL[i] + QLL[[j]][i]
  }
}

####---- Hedge effectiveness ----####

# R^2 98 % 
summary(lm(((1+r)^-(1:max_horizon) %*% t(do.call(rbind,SLL)))[1,] ~ as.matrix((1+r)^-5*do.call(rbind,QLL)[,5])[,1]+
                        as.matrix((1+r)^-10*do.call(rbind,QLL)[,10])[,1]+as.matrix((1+r)^-15*do.call(rbind,QLL)[,15])[,1]+
                        as.matrix((1+r)^-20*do.call(rbind,QLL)[,20])[,1]+as.matrix((1+r)^-25*do.call(rbind,QLL)[,25])[,1]))
