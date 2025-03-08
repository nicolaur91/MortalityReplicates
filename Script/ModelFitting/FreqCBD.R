#' M5 - generalized CBD two-factor mortality model
#' 

# use method as used in Modelling Mortality with Actuarial Applications by MacDonald
X1 = kronecker(diag(iYear),rep(1,iAge))
X2 = kronecker(diag(iYear),unique(dat_EW_mort$Age)-mean(l_age:u_age))

# Column names
colnames(X1) <- paste(rep("k1_",iYear),l_time:u_time,sep="")
colnames(X2) <- paste(rep("k2_",iYear),l_time:u_time,sep="")

X = cbind(X1,X2)

# important to set age and years to factor
set.seed(1991)
CBD.male <- glm(formula=as.formula(dat_EW_mort$mD ~ -1 + X),
                family = poisson,
                offset = log(dat_EW_mort$mE))


# Extracting coefficients
### Extract coefficients - first 1:length(unique(Age)) is a_x   
#### second 
k1_CBD <- coef(CBD.male)[1:iYear]

k2_CBD <- coef(CBD.male)[(iYear+1):length(coef(CBD.male))]

# plot
plot(l_time:u_time,k1_CBD,type="l")
plot(l_time:u_time,k2_CBD,type="l")


# put in data.frame in list
# add common term bx and kt from Li-Lee fit
subpop_k1_CBD <- data.frame(cbind(k1_CBD,l_time:u_time))
names(subpop_k1_CBD) <- c("k1_CBD","Year")
subpop_k2_CBD <- data.frame(cbind(k2_CBD,(l_time:u_time)))
names(subpop_k2_CBD) <- c("k2_CBD","Year")

# add parameters to dataset
temp_dat3 <- dat_EW_mort
temp_dat4 <- merge(temp_dat3,subpop_k1_CBD,by="Year")
temp_dat5 <- merge(temp_dat4,subpop_k2_CBD,by="Year")

# sort by age and then by year:
temp_dat7 <- temp_dat5[order(temp_dat5$Year,temp_dat5$Age),]

#temp_dat7$aic <- LL.male$aic[1]

#temp_dat7$res <- LL.male$residuals

dat_EW_mort <- temp_dat7 


