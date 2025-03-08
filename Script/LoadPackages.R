#' Script for loading packages to use in all of the script
#' 


list.of.packages <- c("HMDHFDplus","ggplot2")

require("HMDHFDplus")
require("StMoMo")
require("expm")
require("matlib")
require("vars")
require(matrixStats) # to do matrix operations and colsds
require("tsDyn") # to simulate VAR-processes
require("rstan") # to fit bayesian model
require("rstanarm") # to fit bayesian model
require("bayesplot") # to do diagnostics


# Install packages
# for (i in 1:length(list.of.packages)) {
#   ifelse(!require(list.of.packages[i]),
#          install.packages(list.of.packages[i]),
#          require(list.of.packages[i]))
# }

