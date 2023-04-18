## originally in LN_modelling_zeros_partialILR_TMB.R

# Simulate data to check convergence

# num_individuals_sim = 100
# d_sim = 4
# n_sim = num_individuals_sim*2
# x_sim= cbind(1, c(rep(0, num_individuals_sim), rep(1, num_individuals_sim)))
# z_sim = give_z_matrix(n_times_2 = num_individuals_sim*2)
# 
# beta_sim = (matrix(runif(d_sim*2, min = -4, max = 4),
#                nrow = 2, byrow=TRUE))
# logs_sd_RE_sim=runif(n = d_sim, min = 0, max = 2)
# cov_RE_sim = runif(n = ((d_sim)*(d_sim)-(d_sim))/2, min = 0.1, max = 0.2)
# sigma_sim = fill_covariance_matrix(d_sim, exp(logs_sd_RE_sim), cov_RE_sim)
# 
# u_large_sim = mvtnorm::rmvnorm(num_individuals_sim, mean = rep(0,d_sim), sigma = sigma_sim)
# 
# ilr_sim = x_sim %*% beta_sim + z_sim %*% u_large_sim
# 
# ## Now convert to compositional data
# 
# clr_recovered = (ilr_sim %*% t(compositions::ilrBase(D = d_sim+1)) %*% diag(d_sim+1))
# 
# inv_clrs = function(clrs) exp((clrs) + log(mean(exp(clrs))))/sum(exp((clrs) + log(mean(exp(clrs)))))
# inv_clrs(clrs)
# 
# probs_recovered = t(apply(clr_recovered, 1, inv_clrs))
# 
# ## Now remove some instances at random
# 
# probs_recovered[sample(1:length(probs_recovered), 100, FALSE)] = 0
# 
# ## Now re-normalise
# 
# probs_recovered = sweep(probs_recovered, 1, rowSums(probs_recovered), '/')
# rowSums(probs_recovered)
# 
# ## Now convert back to the partial ilr
# 
# irl_with_zeros_sim = give_partial_irl(probs_recovered)
# 
# pheatmap(apply(irl_with_zeros_sim, 1, function(i) as.numeric(i==0)))
# 
# TMB_data_sim = list(Y = irl_with_zeros_sim,
#                 num_individuals = num_individuals_sim,
#                 d = d_sim,
#                 n = n_sim,
#                 x = x_sim,
#                 z = z_sim)
# TMB_params_sim = give_TMB_params(d_sim, num_individuals_sim)

simulation_function = function(fraction_zeros, dim=4){
  
  num_individuals_sim = 100
  d_sim = dim
  n_sim = num_individuals_sim*2
  x_sim= cbind(1, c(rep(0, num_individuals_sim), rep(1, num_individuals_sim)))
  z_sim = give_z_matrix(n_times_2 = num_individuals_sim*2)
  
  beta_sim = (matrix(runif(d_sim*2, min = -4, max = 4),
                     nrow = 2, byrow=TRUE))
  logs_sd_RE_sim=runif(n = d_sim, min = 0, max = 2)
  cov_RE_sim = runif(n = ((d_sim)*(d_sim)-(d_sim))/2, min = 0.1, max = 0.2)
  sigma_sim = fill_covariance_matrix(d_sim, exp(logs_sd_RE_sim), cov_RE_sim)
  
  u_large_sim = mvtnorm::rmvnorm(num_individuals_sim, mean = rep(0,d_sim), sigma = sigma_sim)
  
  ilr_sim = x_sim %*% beta_sim + z_sim %*% u_large_sim
  
  ## Now convert to compositional data
  
  clr_recovered = (ilr_sim %*% t(compositions::ilrBase(D = d_sim+1)) %*% diag(d_sim+1))
  
  inv_clrs = function(clrs) exp((clrs) + log(mean(exp(clrs))))/sum(exp((clrs) + log(mean(exp(clrs)))))
  
  probs_recovered = t(apply(clr_recovered, 1, inv_clrs))
  
  ## Now remove some instances at random
  
  probs_recovered[sample(1:length(probs_recovered), fraction_zeros*length(probs_recovered), FALSE)] = 0
  
  
  ## remove any observations in which fewer than two parts of the compositon are greater than 0
  keep_bool = (colSums(apply(probs_recovered, 1, function(i) i>0)) > 2)
  probs_recovered = probs_recovered[keep_bool,]
  x_sim = x_sim[keep_bool,]
  z_sim = z_sim[keep_bool,]
  
  ## remove individuals, if their observations have all been removed
  z_sim = z_sim[,colSums(z_sim) > 0]
  n_sim = sum(keep_bool)
  num_individuals_sim = ncol(z_sim)
  
  ## Now re-normalise
  
  probs_recovered = sweep(probs_recovered, 1, rowSums(probs_recovered), '/')
  rowSums(probs_recovered)
  
  ## Now convert back to the partial ilr
  
  irl_with_zeros_sim = give_partial_irl(probs_recovered)
  
  pheatmap(apply(irl_with_zeros_sim, 1, function(i) as.numeric(i==0)))
  
  
  TMB_data_sim = list(Y = irl_with_zeros_sim,
                      num_individuals = num_individuals_sim,
                      d = d_sim,
                      n = n_sim,
                      x = x_sim,
                      z = z_sim)
  cat('Number of individuals ', num_individuals_sim)
  TMB_params_sim = give_TMB_params(d_sim, num_individuals_sim)
  obj <- MakeADFun(data = TMB_data_sim, parameters = TMB_params_sim, DLL="tmb_MVN_partial_ILR", random = "u_large")
  opt <- do.call("optim", obj)
  opt
  opt$hessian ## <-- FD hessian from optim
  rep <- sdreport(obj)
  return(list(results_inference=rep, true_values=list(beta_sim=beta_sim, logs_sd_RE_sim=logs_sd_RE_sim, cov_RE_sim=cov_RE_sim)))
}
#-------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------#
## simulated data to check convergence
obj <- MakeADFun(data = TMB_data_sim, parameters = TMB_params_sim, DLL="tmb_MVN_partial_ILR", random = "u_large")
opt <- do.call("optim", obj)
opt
opt$hessian ## <-- FD hessian from optim
rep <- sdreport(obj)
rep
#-------------------------------------------------------------------------------------------#

res_sim_01=simulation_function(0.1)
ggplot(cbind.data.frame(true_sim=c(as.vector(res_sim_01$true_values$beta_sim),
                                   res_sim_01$true_values$logs_sd_RE_sim,
                                   res_sim_01$true_values$cov_RE_sim),
                        inferred=res_sim_01$results_inference$par.fixed, names=names(res_sim_01$results_inference$par.fixed)),
       aes(x=true_sim, y=inferred, col=names))+geom_point()+
  geom_abline(intercept = 0, slope = 1, lty='dashed')
ggsave(paste0(folder_images_out, "partialILR_recovery_simulation.pdf"), width = 6, height = 4)

re_run = F
if(re_run){
  # plot(as.vector(res_sim_01$true_values$beta_sim),
  # res_sim_multiple=lapply(rep(seq(0.05, 0.95, by=0.05), each=3), simulation_function)
  fracs_zeros_vector = rep(seq(0.05, 0.85, by=0.1), each=10)
  res_sim_multiple_4=lapply(fracs_zeros_vector, function(i) try(simulation_function(i, dim=4)))
  res_sim_multiple_5=lapply(fracs_zeros_vector, function(i) try(simulation_function(i, dim=5)))
  res_sim_multiple_7=lapply(fracs_zeros_vector, function(i) try(simulation_function(i, dim=7)))
  
  give_cor = function(res_it){
    if(typeof(res_it) != "character"){
      cor(as.vector(res_it$true_values$beta_sim),
          python_like_select_name(res_it$results_inference$par.fixed, 'beta'))
    }else{
      NA
    }
  }
  corrs_all_d4 = sapply(res_sim_multiple_4, give_cor)
  corrs_all_d5 = sapply(res_sim_multiple_5, give_cor)
  corrs_all_d7 = sapply(res_sim_multiple_7, give_cor)
  
  ggplot(cbind.data.frame(cor=corrs_all_d4, fraction_zeros=fracs_zeros_vector), aes(x=fraction_zeros, y=cor, group=fraction_zeros))+geom_boxplot()
  
  ggplot(melt(cbind.data.frame(d4=corrs_all_d4, d5=corrs_all_d5, d7=corrs_all_d7, fraction_zeros=fracs_zeros_vector), id.vars="fraction_zeros"),
         aes(x=fraction_zeros, y=value, col=variable, group=interaction(fraction_zeros, variable)))+geom_boxplot()+
    theme()+guides(col=guide_legend("Number of categories"))+labs(x='Fraction of zero entries', y='Pearson correlation of recovered betas')
  ggsave(paste0(folder_images_out, "recovery_partialILR_simulation.pdf"), width=7, height=3)
  
  saveRDS(list(res_sim_multiple_4, res_sim_multiple_5, res_sim_multiple_7),
          file = paste0(folder_out_RDS, "recovery_partialILR.RDS"))
}


normalise_rw(exposures[,4:7])[patient.meta$PATIENT_ID]

all(rownames(exposures) == patient.meta$SAMPLE_ID)

par(mfrow=c(1,2))
plot_ternary(exposures[patient.meta$group == "rlps",c(1,3,5)], main='arx')
plot_ternary(exposures[patient.meta$group == "arx",c(1,3,5)], main='rlps')

plot_ternary(exposures[patient.meta$group == "arx",c(1,3,4)], main='arx')
plot_ternary(exposures[patient.meta$group == "rlps",c(1,3,4)], main='rlps')

plot_ternary(exposures[patient.meta$group == "arx",c(1,3,6)], main='arx')
plot_ternary(exposures[patient.meta$group == "rlps",c(1,3,6)], main='rlps')

