rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source("../helper/header.R")

exposures_transformed = t(sapply(as.vector(exposures_zeros), function(i) table(factor(i, levels=c(0,1)))))
d = ncol(exposures_zeros)

TMB_data = list(Y = exposures_zeros,
                num_individuals = num_indiv,
                num_sigs = ncol(exposures),
                x = cbind(1, as.numeric(as.factor(patient.meta$group))-1),
                z = give_z_matrix_from_labels(patient.meta$PATIENT_ID))

all(patient.meta$SAMPLE_ID == rownames(exposures))
all(patient.meta$SAMPLE_ID == rownames(exposures_zeros))

give_res_report = function(TMB_data, d){
  TMB_params = list(beta = (matrix(runif(d*2, min = -4, max = 4),
                                   nrow = 2, byrow=TRUE)),
                    u_large = matrix(rep(1, (d)*(TMB_data$num_individuals)), nrow=TMB_data$num_individuals),
                    logs_sd_RE=runif(n = d, min = 0, max = 2),
                    cov_RE_part = 2
  )
  
  
  obj <- MakeADFun(data = TMB_data, parameters = TMB_params, DLL="tmb_correlated_multinom_1", random = "u_large")
  obj$hessian <- TRUE
  opt <- do.call("optim", obj)
  opt
  opt$hessian ## <-- FD hessian from optim
  report = sdreport(obj)
  
  return(list(TMB_params=TMB_params, report=report))
  
}

read_cached = T
if(read_cached){
  results = readRDS("../out/inference/correlated_binom_1_7263302b-a1ea-4010-ba3d-415e655afc3e.RDS")
}else{
  ## This has worked in at least one occasion!! but it also doesn't converge in others
  results = give_res_report(TMB_data, d)
  results
}
report = results$report
TMB_params = results$TMB_params

##' replicate with different starting points. I had a case in which it converged when using independent
##' random effects but then all others didn't!
# multiple_results = replicate(n = 30, give_res_report())

# sapply(multiple_results[2,], function(j) j$pdHess)

image(give_z_matrix(n_times_2 = TMB_data$num_individuals))
length(report$par.random)
matrix(python_like_select_name(report$par.fixed, 'beta'), nrow=2) ## important fixed effects
re_vector_to_matrix(report$par.random, 7) ## negligible RE
(python_like_select_name(report$par.fixed, 'logs_sd_RE'))

fitted_logR = simulate_from_correlated_binom(tmb_fit_object = report, full_RE = T, x_matrix=TMB_data$x, z_matrix=TMB_data$z, return_logratios=T)
fitted_probs = simulate_from_correlated_binom(tmb_fit_object = report, full_RE = T, x_matrix=TMB_data$x, z_matrix=TMB_data$z, return_logratios=F)

dim(fitted_probs)[1] == dim(TMB_data$z)[1] ## it should be the same

colnames(fitted_probs) = colnames(TMB_data$Y)

## the random effects are underspecified in the cases when it doesn't converge

ggplot(melt(fitted_probs))+
  geom_raster( aes( x = Var2, y = Var1, fill = 1-value )  )+
  scale_fill_gradient(low="blue", high="red")+
  geom_point(data=melt((TMB_data$Y)) %>% filter(value == 0),
             aes( x = Var2, y = Var1 ), col='white', shape=18, size=.3)+
  labs(x='Signature', y='Sample')+
  ggtitle('Probabilities of zero exposure and observed zero exposures')
ggsave("../results/zeros_modelling/correlated_binomial_1_v2.pdf", width = 4, height = 4)
# ggsave("../results/zeros_modelling/correlated_binomial_1.pdf", width = 6, height = 6)

##' it only converges if I use no correlations between signatures, i.e. a diagonal covariance matrix. A full
##' covariance matrix doesn't work

results
ggplot(cbind(fitted=melt(fitted_probs), true=melt((TMB_data$Y))), aes(x=fitted.value, y=true.value))+geom_point(aes(col=fitted.value>0.5))+
  facet_wrap(.~fitted.Var2)
ggsave("../results/zeros_modelling/correlated_binomial_1_scatter.pdf", width = 6, height = 6)

## Note: the non-log-zero covariance I have chosen at random!
results$TMB_params$u_large ### THAT DOESN'T MAKE ANY SENSE

betas_inferred = matrix(python_like_select_name(results$report$par.fixed, 'beta'), nrow=2)
rownames(betas_inferred) = c('Intercept', 'Slope')
colnames(betas_inferred) = paste0('s', 1:7)
ggplot(melt(betas_inferred), aes(x=Var2, y=value))+geom_point()+facet_wrap(.~Var1)

## s1 has an extremely high coefficient despite being almost always non-zero.
xtable::xtable(sapply(split(exposures_zeros[,1], patient.meta$group), function(i) table(i)/length(i)*100))

lapply(1:7, function(j) fisher.test(sapply(split(exposures_zeros[,j], patient.meta$group), function(i) table(i))) )

uuid = UUIDgenerate()
save.image(paste0("../out/inference/correlated_binom_1_", uuid, '.Data'))
saveRDS(results, file = paste0("../out/inference/correlated_binom_1_", uuid, '.RDS'))

### Without signature s5
TMB_data_nos5 = TMB_data
TMB_data_nos5$Y = TMB_data_nos5$Y[,-5]
TMB_data_nos5$num_sigs = TMB_data_nos5$num_sigs -1
results_nos5 = give_res_report(TMB_data = TMB_data_nos5, d = TMB_data_nos5$num_sigs)
results_nos5
wald_TMB_wrapper(results_nos5$report)

wald_TMB_wrapper(results$report)
