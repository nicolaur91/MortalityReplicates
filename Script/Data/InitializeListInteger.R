#' Initialize important integers and matrices
#' 

# upper and lower limit of cohorts - perhaps subtract -4 in this?
l_Coh <- min(dat_EW_mort$mC)
u_Coh <- max(dat_EW_mort$mC)
# Number of ages and years
iAge <- length(l_age:u_age)
iYear <- length(l_time:u_time)
iCoh <- length(l_Coh:u_Coh)

# Set weights equal to zero for cohorts with less than 5 observation and for 1886.
# since those models are modelled using glm it make perhaps more sense to do it in  
dat_EW_mort$w <- 1

excl_coh <- c("1886",names(table(dat_EW_mort$mC)[which(table(dat_EW_mort$mC)<5)]))

dat_EW_mort$w[which(dat_EW_mort$mC %in% excl_coh)] <- 0

# mean-age to use in CBD-models
# set in mean age deviation to use in CBD-model
# dat_EW_mort$mAge <- dat_EW_mort$Age-mean(l_age:u_age)

# forecast horizon
horizon <- 30

# number of simulations
nsim <- 10^4

# Risk premium/Sharpe ratio
RP <- 0.2

# premium
prem_q <- matrix(data=NA,nrow=7,ncol=25)