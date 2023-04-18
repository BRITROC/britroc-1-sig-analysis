rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(MASS)
library(mvtnorm)
library(TMB)
library(matrixcalc)
source("../../../../../../GlobalDA/code/2_inference_TMB/helper_TMB.R")

TMB::compile("../../tmb_RE/tmb_MVN.cpp", "-std=gnu++17")
dyn.load(dynlib("../../tmb_RE/tmb_MVN"))

# ## Simulate data with a random correlation matrix
# d <- 2
# A <- matrix(runif(d^2)*2-1, ncol=d)
# Sigma <- t(A) %*% A
# Sigma_without_logs = Sigma
# matrixcalc::is.positive.definite(Sigma)
# Y = mvtnorm::rmvnorm(n = 1000, mean = rep(0,d), sigma = Sigma)
# 
# ## Create TMB objects
# TMB_data = list(Y = Y)
# TMB_params = list(logs_sd_RE=runif(n = d, min = 0, max = 2),
#                   cov_RE= runif(n = ((d)*(d)-(d))/2, min = 0.1, max = 0.2))
# 
# ## Inference
# obj <- MakeADFun(data = TMB_data, parameters = TMB_params, DLL="tmb_MVN")
# obj$hessian <- TRUE
# opt <- do.call("optim", obj)
# opt
# opt$hessian ## <-- FD hessian from optim
# report = sdreport(obj)
# report
# 
# ## using the covariances as they are, and exponentiating logs_sd_RE for the diagonal
# cov_mat_est = give_UNSTRUCTURED_CORR_t_matrix(python_like_select_name(report$par.fixed, "cov_RE"), dim_mat = d)
# diag(cov_mat_est) = exp(python_like_select_name(report$par.fixed, "logs_sd_RE"))
# 
# ## using the covariances as they are, and the logs_sd_RE as they are
# cov_mat_est_2 = give_UNSTRUCTURED_CORR_t_matrix(python_like_select_name(report$par.fixed, "cov_RE"), dim_mat = d)
# diag(cov_mat_est_2) = (python_like_select_name(report$par.fixed, "logs_sd_RE"))
# 
# ## trying to scale column-wise, but that gives a non-symmetrical matrix
# cov_mat_est_2 = give_UNSTRUCTURED_CORR_t_matrix(python_like_select_name(report$par.fixed, "cov_RE"), dim_mat = d)
# cov_mat_est_2= sweep(cov_mat_est_2, 2, exp(python_like_select_name(report$par.fixed, "logs_sd_RE")), '*')
# 
# par(mfrow=c(1,5))
# image(Sigma, main='True')
# image(cov_mat_est, main='Est + log')
# image(cov_mat_est_2, main='Est + unlog')
# image(cov_mat_est_2, main='Est with sweep + log')
# 
# Y_sim = mvtnorm::rmvnorm(1000, mean = rep(0,d), sigma = cov_mat_est)
# 
# par(mfrow=c(1,1))
# plot(Y)
# points(Y_sim, col='blue', pch=2, cex=0.2)
# 
# ## check if the estimated covariance matrices are positive semidefinite
# mvtnorm::rmvnorm(1, mean = rep(0,d), sigma = cov_mat_est) ## PSD
# mvtnorm::rmvnorm(1, mean = rep(0,d), sigma = cov_mat_est_2) ## not PSD
# mvtnorm::rmvnorm(1, mean = rep(0,d), sigma = cov_mat_est_2) ## not even symmetric

## Now doing the same, but with higher-dimensional matrices
d <- 4
A5 <- matrix(runif(d^2)*2-1, ncol=d)
Sigma5 <- t(A5) %*% A5
Sigma_without_logs5 = Sigma5
matrixcalc::is.positive.definite(Sigma5)
Y = mvtnorm::rmvnorm(n = 5000, mean = rep(0,d), sigma = Sigma5)

## Create TMB objects
TMB_data5 = list(Y = Y)
TMB_params5 = list(logs_sd_RE=runif(n = d, min = 0, max = 2),
                  cov_RE= runif(n = ((d)*(d)-(d))/2, min = 0.1, max = 0.2))

obj5 <- MakeADFun(data = TMB_data5, parameters = TMB_params5, DLL="tmb_MVN")
obj5$hessian <- TRUE
opt5 <- do.call("optim", obj5)
opt5
opt5$hessian ## <-- FD hessian from optim
report5 = sdreport(obj5)
report5
report5$env$parameters ## what are these parameters??? they always seem to be positive...
log(report5$env$parameters$cov_RE)

## using the covariances as they are, and exponentiating logs_sd_RE for the diagonal.
cov_mat_est = give_UNSTRUCTURED_CORR_t_matrix(python_like_select_name(report5$par.fixed, "cov_RE"), dim_mat = d)
diag(cov_mat_est) = exp(python_like_select_name(report5$par.fixed, "logs_sd_RE"))

## using the covariances as they are, and the logs_sd_RE as they are
cov_mat_est_2 = give_UNSTRUCTURED_CORR_t_matrix(python_like_select_name(report5$par.fixed, "cov_RE"), dim_mat = d)
diag(cov_mat_est_2) = (python_like_select_name(report5$par.fixed, "logs_sd_RE"))

#' using the covariances as they are, and exponentiating logs_sd_RE for the diagonal,
#' and then taking its power, to go from sd to variance. THIS IS AFAIC THE CORRECT VERSION
cov_mat_est_3 = give_UNSTRUCTURED_CORR_t_matrix(python_like_select_name(report5$par.fixed, "cov_RE"), dim_mat = d)
diag(cov_mat_est_3) = exp(python_like_select_name(report5$par.fixed, "logs_sd_RE"))**2

## trying to scale column-wise, but that gives a non-symmetrical matrix
# cov_mat_est_2 = give_UNSTRUCTURED_CORR_t_matrix(python_like_select_name(report5$par.fixed, "cov_RE"), dim_mat = d)
# cov_mat_est_2= sweep(cov_mat_est_2, 2, exp(python_like_select_name(report5$par.fixed, "logs_sd_RE")), '*')

## using env$parameters, which is not clear what they are
cov_mat_est_4 <- fill_covariance_matrix(d,
                       (exp(report5$env$parameters$logs_sd_RE))**2,
                       report5$env$parameters$cov_RE, v=F)

## Rescale the values (correlations) to get covariances
## get the values of the correlation matrix
cov_mat_est_5 = give_UNSTRUCTURED_CORR_t_matrix(python_like_select_name(report5$par.fixed, "cov_RE"), dim_mat = d)
## for each element, get the covariance from the correlation, by multiplying by the standard deviations
.stdevs <- python_like_select_name(report5$par.fixed, "logs_sd_RE")
for(i in 1:d){
  for(j in i:d){
    cov_mat_est_5[i,j] = cov_mat_est_5[j,i] =cov_mat_est_5[i,j]*.stdevs[i]*.stdevs[j]
  }
}
## numerically it might give a matrix which is not symmetric
## we don't need to add the variances in the diagonal because we already haev them:
all(diag(cov_mat_est_5) == .stdevs**2)

##----------
## THE CORRECT WAY
L = give_UNSTRUCTURED_CORR_t_matrix(python_like_select_name(report5$par.fixed, "cov_RE"), dim_mat = d)
D = diag(L %*% t(L))
Dmat2 <- matrix(0, nrow = d, ncol = d)
diag(Dmat2) <- D
sqrt(solve(Dmat2))
Dmat2 <- solve(sqrt(Dmat2))
Sigma <- (Dmat2) %*% L %*% t(L) %*% (Dmat2)
Sigma
diag(Sigma) <- exp(python_like_select_name(report5$par.fixed, "logs_sd_RE"))**2
Sigma
mvtnorm::rmvnorm(1, mean = rep(0,d), sigma = Sigma) ## PSD

##

## check if the estimated covariance matrices are positive semidefinite
mvtnorm::rmvnorm(1, mean = rep(0,d), sigma = cov_mat_est) ## not PSD
mvtnorm::rmvnorm(1, mean = rep(0,d), sigma = cov_mat_est_2) ## not PSD
mvtnorm::rmvnorm(1, mean = rep(0,d), sigma = cov_mat_est_3) ## PSD in general, but not always!
#' Warning message:
#' In mvtnorm::rmvnorm(1, mean = rep(0, d), sigma = cov_mat_est_3) :
#'   sigma is numerically not positive semidefinite
mvtnorm::rmvnorm(1, mean = rep(0,d), sigma = cov_mat_est_5)  ## not PSD

## making a clearly non PSD matrix
cov_mat_est_3_mod <- cov_mat_est_3
cov_mat_est_3_mod[2,1] = cov_mat_est_3_mod[1,2] = -400
mvtnorm::rmvnorm(1, mean = rep(0,d), sigma = cov_mat_est_3_mod) ## PSD in general, but not always!

par(mfrow=c(1,6))
image(Sigma5, main='True')
image(cov_mat_est, main='Est + log')
image(cov_mat_est_2, main='Est + unlog')
image(cov_mat_est_3, main='Est + log + squared')
image(cov_mat_est_4, main='from env params')
image(cov_mat_est_5, main='from env params')


report5$pdHess
opt5$convergence

matrixcalc::is.positive.definite(Sigma5)
matrixcalc::is.positive.definite(cov_mat_est)
matrixcalc::is.positive.definite(cov_mat_est_2)
matrixcalc::is.positive.definite(cov_mat_est_3)
matrixcalc::is.positive.definite(cov_mat_est_5)

my_line <- function(x,y,...){
  points(x,y,...)
  abline(a = 0,b = 1,...)
}
pairs((cbind.data.frame(true=as.vector(Sigma5), as.vector(cov_mat_est),
                         as.vector(cov_mat_est_2), as.vector(cov_mat_est_3),
                        as.vector(cov_mat_est_5))),
      lower.panel = my_line, upper.panel = my_line)


Matrix::nearPD(cov_mat_est_3)$mat
Sigma5

report5$par.fixed
report5$env$parameters ## ??? are these the reported values? they seem to be different.
## re-run without reporting and see that happens

report5$env$random

summary_tmb <- summary(report5)
report5$par.fixed

opt5$convergence
opt5$hessian
report5$par.fixed
