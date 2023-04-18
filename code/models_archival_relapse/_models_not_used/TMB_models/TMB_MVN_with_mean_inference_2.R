rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(MASS)
library(mvtnorm)
library(TMB)
library(matrixcalc)
source("../../../GlobalDA/code/2_inference_TMB/helper_TMB.R")

TMB::compile("tmb_RE/tmb_MVN_with_mean_2.cpp")
dyn.load(dynlib("tmb_RE/tmb_MVN_with_mean_2"))

## Simulate data with a random correlation matrix
d <- 2
A <- matrix(runif(d^2)*2-1, ncol=d)
Sigma <- t(A) %*% A
Sigma_without_logs = Sigma
matrixcalc::is.positive.definite(Sigma)
mu = runif(d)
Y = mvtnorm::rmvnorm(n = 1000, mean = mu, sigma = Sigma)



## Create TMB objects
TMB_data = list(Y = Y,
                x = matrix(rep(1, nrow(Y))))
TMB_params = list(logs_sd_RE=runif(n = d, min = 0, max = 2),
                  cov_RE= runif(n = ((d)*(d)-(d))/2, min = 0.1, max = 0.2),
                  beta = matrix(runif(d, min = -4, max = 4), nrow = 1))

## Inference
obj <- MakeADFun(data = TMB_data, parameters = TMB_params, DLL="tmb_MVN_with_mean_2")
obj$hessian <- TRUE
opt <- do.call("optim", obj)
opt
opt$hessian ## <-- FD hessian from optim
report = sdreport(obj)
report
