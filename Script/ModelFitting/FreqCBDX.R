#' script for ML-estimates of CBDX-model

# Add age to coloumns
X_age <-  kronecker(rep(1,iYear),diag(iAge))

colnames(X_age) <- paste(rep("Bx_",iAge),l_age:u_age,sep="")

XX <- cbind(X_age,X)

# run model
CBD.freqX <- glm(formula=as.formula(dat_EW_mort$mD ~ -1 + XX),
                family = poisson,
                offset = log(dat_EW_mort$mE))

# retrieve ML-estimates

# Extracting coefficients
### Extract coefficients - first 1:length(unique(Age)) is a_x   
#### second 
Bx_CBDX <- coef(CBD.freqX)[1:iAge]
k1_CBDX <- coef(CBD.freqX)[(iAge+1):(iAge+iYear)]

k2_CBDX <- coef(CBD.freqX)[(iAge+iYear+1):length(coef(CBD.freqX))]

# plot
plot(l_age:u_age,Bx_CBDX,type="l")
plot(l_time:u_time,k1_CBDX,type="l")
plot(l_time:u_time,k2_CBDX,type="l")