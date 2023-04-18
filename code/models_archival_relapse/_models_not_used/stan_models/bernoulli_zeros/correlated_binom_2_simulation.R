rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

### Depending on dataset it works or not
### When using ubsets of the signatures, depending on the RUN
## (i.e. the initial values) it converges or not (!!!!)

source("header.R")
library(MASS)

# TMB::compile("tmb_RE/tmb_correlated_multinom_2_allFE_b.cpp")
# dyn.load(dynlib("tmb_RE/tmb_correlated_multinom_2_allFE_b"))

# exposures_zeros = exposures_zeros[,-1]
# exposures_zeros = as.matrix(exposures_zeros[,3])

# exposures_transformed = t(sapply(as.vector(exposures_zeros), function(i) table(factor(i, levels=c(0,1)))))
simulation = F
if(simulation){
  exposures_zeros = cbind(sample(x = c(0,1), size = nrow(exposures_zeros), replace = T),
                          sample(x = c(0,1), size = nrow(exposures_zeros), replace = T)) ## two-dimensional
}

exposures_zeros = exposures_zeros[,1:2]

d = ncol(exposures_zeros)

TMB_data = list(Y = exposures_zeros,
                num_individuals = num_indiv,
                num_sigs = ncol(exposures),
                x = cbind(1, as.numeric(as.factor(patient.meta$group))-1),
                z = give_z_matrix_from_labels(patient.meta$PATIENT_ID))


give_res_report = function(TMB_data, d){
  TMB_params = list(beta = (matrix(runif(d*2, min = -4, max = 4),
                                   nrow = 2, byrow=TRUE)),
                    u_large = matrix(rep(1, (d)*(TMB_data$num_individuals)), nrow=TMB_data$num_individuals),
                    logs_sd_RE=runif(n = d, min = 0, max = 2),
                    cov_RE= runif(n = ((d)*(d)-(d))/2, min = 0.1, max = 0.2)
  )
  
  
  obj <- MakeADFun(data = TMB_data, parameters = TMB_params, DLL="tmb_correlated_multinom_2_allFE", random = "u_large")
  # obj <- MakeADFun(data = TMB_data, parameters = TMB_params, DLL="tmb_correlated_multinom_2_allFE_b")
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


## Do the same with several observations
num_replicates = 4
beta_true = matrix(c(1,1,0.2, 2), nrow=2, byrow = T)
x_with_replicates = do.call('rbind', replicate(num_replicates, cbind(1, as.numeric(as.factor(patient.meta$group))-1), simplify = F))
z_with_replicates = do.call('rbind', replicate(num_replicates, give_z_matrix_from_labels(patient.meta$PATIENT_ID), simplify = F))
u_true = rmvnorm(n = num_indiv, mean = c(0,0), sigma = diag(c(1,1)))
## treat each value j as a log-ratio, individually
exposures_with_replicates = apply(x_with_replicates %*% beta_true + z_with_replicates %*% u_true, 2, function(i){
  sapply(i, function(j) {prob_j = exp(j)/(1+exp(j)); sample(x = c(0,1), size = 1, prob = c(1-prob_j, prob_j))} )
  })

TMB_data_with_replicates = list(Y = exposures_with_replicates,
                num_individuals = num_indiv,
                num_sigs = ncol(exposures_with_replicates),
                x = x_with_replicates,
                z = z_with_replicates)

results_with_replicates = give_res_report(TMB_data_with_replicates, ncol(exposures_with_replicates))
results_with_replicates ### same problem if we have multiple observations for the same patient

# report = results$report
# TMB_params = results$TMB_params

# 
# fitted_logR = simulate_from_correlated_binom(tmb_fit_object = report, full_RE = T, x_matrix=TMB_data$x, z_matrix=TMB_data$z, return_logratios=T)
# fitted_probs = simulate_from_correlated_binom(tmb_fit_object = report, full_RE = T, x_matrix=TMB_data$x, z_matrix=TMB_data$z, return_logratios=F)
# colnames(fitted_probs) = colnames(TMB_data$Y)
# 
# 
# ggplot(melt(fitted_probs))+
#   geom_raster( aes( x = Var2, y = Var1, fill = 1-value )  )+
#   scale_fill_gradient(low="blue", high="red")+
#   geom_point(data=melt((TMB_data$Y)) %>% filter(value == 0),
#              aes( x = Var2, y = Var1 ), col='white', shape=18, size=.3)+
#   labs(x='Signature', y='Sample')+
#   ggtitle('Probabilities of zero exposure and observed zero exposures')
# # ggsave("../results/zeros_modelling/correlated_binomial_2_v2.pdf", width = 4, height = 4)
# # ggsave("../results/zeros_modelling/correlated_binomial_2.pdf", width = 6, height = 6)
# 
# 
# ggplot(cbind(fitted=melt(fitted_probs), true=melt((TMB_data$Y))), aes(x=fitted.value, y=true.value))+geom_point(aes(col=fitted.value>0.5))+
#   facet_wrap(.~fitted.Var2)
# # ggsave("../results/zeros_modelling/correlated_binomial_2_scatter.pdf", width = 6, height = 6)
# 
# ##' replicate with different starting points. I had a case in which it converged when using independent
# ##' random effects but then all others didn't!
# # multiple_results = replicate(n = 30, give_res_report())
# 
# # sapply(multiple_results[2,], function(j) j$pdHess)
# 
# image(give_z_matrix(n_times_2 = TMB_data$num_individuals))
# length(report$par.random)
# matrix(python_like_select_name(report$par.fixed, 'beta'), nrow=2) ## important fixed effects
# re_vector_to_matrix(report$par.random, 7) ## negligible RE
# (python_like_select_name(report$par.fixed, 'logs_sd_RE'))
# fitted_cov_mat = give_UNSTRUCTURED_CORR_t_matrix(results$TMB_params$cov_RE, ncol(TMB_data$Y))
# diag(fitted_cov_mat) = exp(results$TMB_params$logs_sd_RE)
# 
# (python_like_select_name(report$par.fixed, 'logs_sd_RE'))
# 
# 
# fitted_logR = simulate_from_correlated_binom(tmb_fit_object = report, full_RE = T, x_matrix=TMB_data$x, z_matrix=TMB_data$z, return_logratios=T)
# fitted_probs = simulate_from_correlated_binom(tmb_fit_object = report, full_RE = T, x_matrix=TMB_data$x, z_matrix=TMB_data$z, return_logratios=F)
# 
# dim(fitted_probs)[1] == dim(TMB_data$z)[1] ## it should be the same
# 
# colnames(fitted_probs) = colnames(TMB_data$Y)
# 
# ## the random effects are underspecified in the cases when it doesn't converge
# 
# ggplot(melt(fitted_probs))+
#   geom_raster( aes( x = Var2, y = Var1, fill = 1-value )  )+
#   scale_fill_gradient(low="blue", high="red")+
#   geom_point(data=melt((TMB_data$Y)) %>% filter(value == 0),
#              aes( x = Var2, y = Var1 ), col='white', shape=18, size=.3)+
#   labs(x='Signature', y='Sample')+
#   ggtitle('Probabilities of zero exposure and observed zero exposures')
# # ggsave("../results/zeros_modelling/correlated_binomial_2.pdf", width = 4, height = 4)
# 
# ##' it only converges if I use no correlations between signatures, i.e. a diagonal covariance matrix. A full
# ##' covariance matrix, which I am using here, doesn't work
# 
# results
# ggplot(cbind(fitted=melt(fitted_probs), true=melt((TMB_data$Y))), aes(x=fitted.value, y=true.value))+geom_point(aes(col=fitted.value>0.5))+
#   facet_wrap(.~fitted.Var2)
# 
# nsamples_sim = 200
# npatients_sim = 10
# x_sim = t(cbind(1, sample(c(0,1), size = nsamples_sim, replace = T)))
# z_sim = t(do.call('rbind', replicate(diag(10), n = nsamples_sim/npatients_sim, simplify = F)))
# A <- matrix(runif(d^2)*2-1, ncol=d) 
# cov_mat =  t(A) %*% A
# u_sim = MASS::mvrnorm(n = npatients_sim, mu = rep(0,ncol(TMB_data$Y)), Sigma = cov_mat)
# beta_sim = rbind(runif(n = d),runif(n = d))
# theta_sim = t(x_sim) %*% beta_sim + t(z_sim) %*% u_sim
# exposures_zeros_sim = t(apply(theta_sim, 1, function(i) sapply(i, function(j) sample(x = c(1,0), size = 1, prob = softmax(c(j, 0))))))
# TMB_data_sim = list(Y = exposures_zeros_sim,
#                     num_individuals = npatients_sim,
#                     num_sigs = ncol(exposures_zeros_sim),
#                     x = t(x_sim),
#                     z = t(z_sim))
# results_sim = give_res_report(TMB_data_sim, d = ncol(TMB_data_sim$Y))
# results_sim
# 
# 
# apply(z_sim, 1, sum)
# #### THIS IS NOT CORRECT!! THESE ARE THE INITIAL PARAMETERS, NOT THE FITTED ONES!
# fitted_cov_mat_sim = (give_UNSTRUCTURED_CORR_t_matrix(results_sim$TMB_params$cov_RE, ncol(TMB_data_sim$Y)))
# diag(fitted_cov_mat_sim) = exp(results_sim$TMB_params$logs_sd_RE)
# 
# par(mfrow=c(4,2))
# plot(0,0)
# image(fitted_cov_mat)
# image(cov_mat)
# image(fitted_cov_mat_sim)
# 
# #### WHY IS COV_RE ALWAYS POSITIVE? IT DOESN'T NEED TO BE!
# #### Are they around the initial value?
# hist(as.vector(fitted_cov_mat))
# hist(as.vector(fitted_cov_mat))
# 
# plot(as.vector(cov_mat), as.vector(fitted_cov_mat), col=factor(as.vector(diag(ncol(TMB_data$Y)))), pch=10)
# abline(coef = c(0,1))
# plot(as.vector(cov_mat), as.vector(fitted_cov_mat_sim),  col=factor(as.vector(diag(ncol(TMB_data$Y)))), pch=10)
# abline(coef = c(0,1))
# 
# ## Compare matrices of random effects
# par(mfrow=c(1,1))
# plot(as.vector(u_sim), as.vector(matrix(results_sim$report$par.random, ncol=7)))
# abline(coef = c(0,1))
# 
# # uuid = UUIDgenerate()
# # save.image(paste0("../out/inference/LN_fullRE_", uuid, '.Data'))
# # saveRDS(results, file = paste0("../out/inference/LN_fullRE_", uuid, '.RDS'))
