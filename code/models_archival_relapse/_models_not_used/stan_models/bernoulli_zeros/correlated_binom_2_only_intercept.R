rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source("header.R")
library(MASS)
library(mvtnorm)

TMB::compile("tmb_RE/tmb_correlated_multinom_only_intercepts.cpp", "-std=gnu++17")
dyn.load(dynlib("tmb_RE/tmb_correlated_multinom_only_intercepts"))

# exposures_zeros = exposures_zeros[,-1]
# exposures_zeros = as.matrix(exposures_zeros[,3])

exposures_transformed = t(sapply(as.vector(exposures_zeros), function(i) table(factor(i, levels=c(0,1)))))
# exposures_zeros = matrix(sample(x = c(0,1), size = length(exposures_zeros), replace = T)) ## one-dimensional
exposures_zeros = cbind(sample(x = c(0,1), size = nrow(exposures_zeros), replace = T),
                        sample(x = c(0,1), size = nrow(exposures_zeros), replace = T)) ## two-dimensional
d = ncol(exposures_zeros)

TMB_data = list(Y = exposures_zeros,
                num_individuals = num_indiv,
                num_sigs = ncol(exposures),
                x = matrix(rep(1, nrow(exposures))),
                z = give_z_matrix_from_labels(patient.meta$PATIENT_ID))


give_res_report = function(TMB_data, d){
  TMB_params = list(beta = matrix(runif(d, min = -4, max = 4), nrow = 1),
                    # u_large = matrix(rep(1, (d)*(TMB_data$num_individuals)), nrow=TMB_data$num_individuals),
                    logs_sd_RE=runif(n = d, min = 0, max = 2),
                    cov_RE= runif(n = ((d)*(d)-(d))/2, min = 0.1, max = 0.2)
  )
  
  
  obj <- MakeADFun(data = TMB_data, parameters = TMB_params, DLL="tmb_correlated_multinom_only_intercepts")
  obj$hessian <- TRUE
  opt <- do.call("optim", obj)
  opt
  opt$hessian ## <-- FD hessian from optim
  report = sdreport(obj)
  
  return(list(TMB_params=TMB_params, report=report))
  
}

read_cached=F
image(t(TMB_data$Y))
if(read_cached){
}else{
  results = give_res_report(TMB_data, d)
}
results

## now with replicates
## Do the same with several observations
num_replicates = 4
beta_true = matrix(c(1,1,0.2, 2), nrow=2, byrow = T)
x_with_replicates = do.call('rbind', replicate(num_replicates, cbind(1, as.numeric(as.factor(patient.meta$group))-1), simplify = F))
z_with_replicates = do.call('rbind', replicate(num_replicates, give_z_matrix_from_labels(patient.meta$PATIENT_ID), simplify = F))
u_true = mvtnorm::rmvnorm(n = num_indiv, mean = c(0,0), sigma = diag(c(1,1)))
## treat each value j as a log-ratio, individually
exposures_with_replicates = apply(x_with_replicates %*% beta_true + z_with_replicates %*% u_true, 2, function(i){
  sapply(i, function(j) {prob_j = exp(j)/(1+exp(j)); sample(x = c(0,1), size = 1, prob = c(1-prob_j, prob_j))} )
})

TMB_data_with_replicates = list(Y = exposures_with_replicates,
                                num_individuals = num_indiv,
                                num_sigs = ncol(exposures_with_replicates),
                                x = matrix(x_with_replicates[,1]),
                                z = z_with_replicates)

results_with_replicates = give_res_report(TMB_data_with_replicates, ncol(exposures_with_replicates))
results_with_replicates ### same problem if we have multiple observations for the same patient


### drawing from a multivariate normal
library(mvtnorm)
library(matrixcalc)
without_logs = diag(3)
without_logs[2,1] = without_logs[1,2] = 3
with_logs = without_logs
diag(with_logs) = c(1,3,2)
mvtnorm::rmvnorm(n = 100, mean = c(0,0,0), sigma = without_logs)
is.positive.definite(without_logs)
is.positive.definite(with_logs)

##
n <- 4  
A <- matrix(runif(n^2)*2-1, ncol=n) 
Sigma <- t(A) %*% A
Sigma_without_logs = Sigma
diag(Sigma_without_logs) = 1
is.positive.definite(Sigma)
is.positive.definite(Sigma_without_logs)
