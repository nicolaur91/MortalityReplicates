#' Script til at estimere en CBD-model vha. først GLM, dernæst vha. Bayesian ved brug af stan_glm
#' Download data fra HMD - UK population males.

# Clean working environment
rm(list=ls())

# Set working directory
setwd("C:/Users/nsl.fi/Desktop/Backup data/Mortality replicate/Bayesian CBD/Script")
dir_wd <- getwd()

# Load packages
source("LoadPackages.R")

#### Download data - From Human Mortality Database ####
# Load data - only use crude mortality rates, i.e. no smoothing
setwd(paste(dir_wd,"Data",sep="/"))
source("HMD.R")

#### Combine data - ####
#### input: matrix: N different vec ####
#### output: matrix: dat_mat with N coloumns ####
source("CombineDataHMD.R")


source("InitializeListInteger.R")

#### Fitting model - ####
#### i.e. Lee-Carter for whole population ####
setwd(paste(dir_wd,"ModelFitting",sep="/"))
source("FreqCBD.R")

# source("BayesCBD.R") # Bayesian udgave af CBD-model - giver samme point estimates som frek.

# CBD-X models
# source("FreqCBDX.R") # ML-estimates of CBD
# source("BayesCBDX.R") # Bayesian udgave af CBD-X - model

####---- Forecasting model ----####
setwd(paste(dir_wd,"ModelForecast",sep="/"))
source("ForeFreqCBD.R")
source("SimFreqCBD.R")
source("ForeSimBayesCBD.R")