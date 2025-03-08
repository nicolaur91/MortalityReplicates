#' Make 1 dataset for exposures and deaths for EW males population
#' 

# Combine data from HMD
# Set fixed ages and time horizon
# l_time = 1961
l_time = 1977
u_time = 2008
# l_age = 60
# u_age = 84

l_age = 67
u_age = 95


#####--------- Subset Death and exposure matrix for EW males -------#######
dat_mat_d <- dat_EW_d[which(dat_EW_d$Year %in% ((l_time-1)+seq(l_time:u_time)) &
                              dat_EW_d$Age %in% ((l_age-1)+seq(l_age:u_age))),c(1,2,4)]

# rename coloumns
names(dat_mat_d)[3] <- "D"

#####--------- Central exposure matrix for each country -------#######
dat_mat_expc <- dat_EW_expc[which(dat_EW_expc$Year %in% ((l_time-1)+seq(l_time:u_time)) &
                                    dat_EW_expc$Age %in% ((l_age-1)+seq(l_age:u_age))),c(1,2,4)]


# rename coloumns
names(dat_mat_expc)[3] <- "E"


# Total Dataset to gnm model
dat_EW_mort <- data.frame(Age = dat_mat_expc$Age,
                         Year = dat_mat_expc$Year,
                         mD = dat_mat_d$D,
                         mE = dat_mat_expc$E,
                         mC = dat_mat_expc$Year-dat_mat_expc$Age)


# Number of ages and years
iAge <- length(l_age:u_age)
iYear <- length(l_time:u_time)

#####----- leave everything out exc. matrix and average -----######
rm(list=setdiff(ls(), c("dir_wd","dat_EW_mort","l_time","u_time",
                        "l_age","u_age", "iAge","iYear")))