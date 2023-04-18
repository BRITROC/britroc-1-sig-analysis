## Infence using stan_fit_LNM for data in the simplex
## Simple LN, no fixed or random effects

rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# set.seed(1234)

library(uuid)
library(ggplot2)
library(reshape2)
library(compositions)
library(ggrepel)
library(TMB)
library(gridExtra)
library(zCompositions) ## for model-based imputation
library(robCompositions) ## for model-based imputation

source("../../../../../GlobalDA/code/2_inference_TMB/helper_TMB.R")
source("../helper/functions.R")
source("../helper/header.R")
source("../../../../Vias_Brenton/copy_number_analysis_organoids/helper_functions.R")
#-------------------------------------------------------------------------------------------#

all(patient.meta$SAMPLE_ID == rownames(exposures))

folder_out_RDS <- "../../../out/inference/partialILR/"
folder_images_out <- "../../../results/partialILRmodelling_ME/"
system(paste0("mkdir -p ", folder_out_RDS))
system(paste0("mkdir -p ", folder_images_out))

#-------------------------------------------------------------------------------------------#

## Transform data with partial ILR
exposures_imput1em2 <- impute(exposures, 1e-2)
which_zero = t(apply(exposures_imput1em2, 1, function(i) as.numeric((i==0)) ))

all(rownames(exposures_imput1em2) == patient.meta$SAMPLE_ID)
give_pca(exposures_imput1em2,
         names=rownames(exposures_imput1em2), groups=patient.meta$group,
         additional_df = data.frame(patient=patient.meta$PATIENT_ID))+theme_bw()#+
  # geom_line(aes(group=patient))

exposures_imput1em2[1,]
which_zero[1,]
as.vector(compositions::ilr(exposures_imput1em2)[1,])
give_partial_ilr_basis(which(which_zero[1,] == 1), d=7)
give_partial_ilr_basis(which_zero_vector = 2, d = 5)

prepare_TMB_data_with_subset = function(subset_sigs){
  .exps = normalise_rw(exposures_imput1em2[,subset_sigs])
  .irl_with_zeros = give_partial_irl(.exps)
  .keep = (rowSums(.irl_with_zeros == 0) < (ncol(.exps) - 2) )
  .z_britroc = give_z_matrix_from_labels(patient.meta$PATIENT_ID[.keep])
  .irl_with_zeros = .irl_with_zeros[.keep,]  ## if a sample has too many zeros (i.e. all but 1 or 2) remove
  
  TMB_data_sim = list(Y = .irl_with_zeros,
                      num_individuals = ncol(.z_britroc),
                      d = length(subset_sigs)-1,
                      n = nrow(.exps),
                      x = cbind(1, as.numeric(as.factor(patient.meta$group[.keep]))-1),
                      z = .z_britroc)
  return(TMB_data_sim)
}

exposures_imput1em2_partial_irl = give_partial_irl(exposures_imput1em2)

dim(exposures_imput1em2)
dim(exposures_imput1em2_partial_irl)

## note that the partial ILR that I have  computed is different from the ILR setting zero entries as zero ILRs
## different
# compositions::clr(exposures_imput1em2[1,]) %*% give_partial_ilr_basis(which_zero[[1]])
# compositions::clr(exposures_imput1em2[1,]) %*% irl_base_complete

## i.e the zero entries are different

# To show that I am computing the ILR correctly given a basis
# plot(as.vector(compositions::ilr(exposures_imput1em2)[1,]),
#      compositions::clr(exposures_imput1em2[1,]) %*% irl_base_complete)
# abline(coef = c(0,1))

#-------------------------------------------------------------------------------------------#
# TMB::compile("../tmb_RE/tmb_MVN_partial_ILR.cpp", "-std=gnu++17")
# dyn.load(dynlib("../tmb_RE/tmb_MVN_partial_ILR"))
## with correlations
TMB::compile("../tmb_RE/tmb_MVN_ILR.cpp", "-std=gnu++17")
dyn.load(dynlib("../tmb_RE/tmb_MVN_ILR"))

## ME of ILR or ALR without zeros; no correlations
TMB::compile("../tmb_RE/tmb_RE_20220222.cpp", "-std=gnu++17")
dyn.load(dynlib("../tmb_RE/tmb_RE_20220222"))

#-------------------------------------------------------------------------------------------#


#-------------------------------------------------------------------------------------------#
## remove the samples in which there are too many zeros. There can be at most 4 zeros (i.e. there has to be at least one log-ratio)
keep_britroc_samples = (rowSums(exposures_imput1em2_partial_irl == 0) < (ncol(exposures_imput1em2) - 2) )
exposures_imput1em2_partial_irl = exposures_imput1em2_partial_irl[keep_britroc_samples,]
all(patient.meta$SAMPLE_ID[keep_britroc_samples] == rownames(exposures_imput1em2_partial_irl))
z_britroc = give_z_matrix_from_labels(patient.meta$PATIENT_ID[keep_britroc_samples])
TMB_data = list(Y = exposures_imput1em2_partial_irl,
                num_individuals = ncol(z_britroc),
                d = d,
                n = nrow(exposures_imput1em2_partial_irl),
                x = cbind(1, as.numeric(as.factor(patient.meta$group[keep_britroc_samples]))-1),
                z = z_britroc)

give_TMB_data <- function(imputation_val, transformation='ALR'){
  if(is.numeric(imputation_val)){
    exposures_impute <- impute(exposures, imputation_val)
  }else{
    if(imputation_val == 'multLN'){
      exposures_impute <- zCompositions::multLN(exposures, label = 0, dl = rep(0.5, 7))
    }else if(imputation_val == 'multRepl'){
      exposures_impute <- zCompositions::multRepl(exposures, label = 0, dl = rep(0.5, 7))
      exposures_impute[exposures_impute < 0] <- 1e-2
      exposures_impute <- normalise_rw(exposures_impute)
    }else if(imputation_val == 'imputeBDLs'){
      exposures_impute <- robCompositions::imputeBDLs(exposures,dl = rep(0.05, 7), eps=1, method="lmrob")
      
    }else{
      stop('Specify correct imputation_val')
    }
  }
  
  if(transformation == 'ILR'){
    exposures_impute_partial_irl <- as(compositions::ilr(exposures_impute), 'matrix')    
  }else if(transformation == 'ALR'){
    exposures_impute_partial_irl <- as(compositions::alr(exposures_impute), 'matrix')
  }else if(transformation == 'none'){
    exposures_impute_partial_irl <- exposures_impute
  }else{
    stop()
  }

  list(Y = exposures_impute_partial_irl,
                  num_individuals = ncol(z_britroc),
                  d = d,
                  n = nrow(exposures_impute_partial_irl),
                  x = cbind(1, as.numeric(as.factor(patient.meta$group[keep_britroc_samples]))-1),
                  z = z_britroc)
}

sapply(TMB_data, dim)
sapply(TMB_params, dim)
#-------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------#
## TRUE_data: full
TMB_params = give_TMB_params(d, TMB_data$num_individuals)
obj <- MakeADFun(data = TMB_data, parameters = TMB_params, DLL="tmb_MVN_ILR", random = "u_large")
opt <- do.call("optim", obj)
opt
opt$hessian ## <-- FD hessian from optim
rep <- sdreport(obj)
rep

give_pairs_with_mvn_wrapper(matrix(rep$par.random, ncol=TMB_data$d), common_lims = F)
rep_nlminb <- wrapper_run_TMB_use_nlminb(model = "tmb_MVN_ILR",
                                         object = TMB_data, use_nlminb = T, iter.max=10000)
TMB::compile("../tmb_RE/tmb_RE_20220222.cpp", "-std=gnu++17")
dyn.load(dynlib("../tmb_RE/tmb_RE_20220222"))

rep_nlminbnocor <- wrapper_run_TMB_use_nlminb(model = "tmb_MVN_ILRnocor",
                                         object = TMB_data, use_nlminb = T, iter.max=10000)
give_pairs_with_mvn_wrapper(matrix(rep_nlminb$par.random, ncol=TMB_data$d), common_lims = F)
give_pairs_with_mvn_wrapper(matrix(rep_nlminbnocor$par.random, ncol=TMB_data$d), common_lims = F)
compare_TMB_fit_to_data(rep_nlminb, TMB_data, remove_zeros=T, title='res_nlminb')
compare_TMB_fit_to_data(rep_nlminbnocor, TMB_data, remove_zeros=T, title='res_nlminb')
## there still quite a severe bias, e.g. in ILR1

# RE_res_nlminb <- matrix(rep_nlminb$par.random, ncol=TMB_data$d)
# fitted_res_nlminb <- TMB_data$x %*% matrix(python_like_select_name(rep_nlminb$par.fixed, 'beta'), nrow=2)+
#   TMB_data$z %*% RE_res_nlminb
# compare_matrices(mat1 = TMB_data$Y,
#                  mat2 = add_rownames(add_colnames(fitted_res_nlminb, colnames(TMB_data$Y)),
#                                      rownames(TMB_data$Y)), remove_zeros=T) ## removing zero entries
## much better fit than when using partial ILR

compare_TMB_fit_to_data(rep_nlminb, TMB_data)

repeat_full_with_different_initial = function(){
  TMB_params = give_TMB_params(d, TMB_data$num_individuals)
  obj <- MakeADFun(data = TMB_data, parameters = TMB_params, DLL="tmb_MVN_partial_ILR", random = "u_large")
  opt <- do.call("optim", obj)
  opt
  opt$hessian ## <-- FD hessian from optim
  rep <- sdreport(obj)
  rep
}

run_several_imputation <- function(imputation_val){
  exposures_imput1em2 <- impute(exposures_imput1em2, imputation_val)
  exposures_imput1em2_partial_irl = give_partial_irl(exposures_imput1em2)
  
  TMB_data = list(Y = exposures_imput1em2_partial_irl,
                  num_individuals = ncol(z_britroc),
                  d = d,
                  n = nrow(exposures_imput1em2_partial_irl),
                  x = cbind(1, as.numeric(as.factor(patient.meta$group))-1),
                  z = z_britroc)

  TMB_params = give_TMB_params(d, TMB_data$num_individuals)
  rep <- list()
  rep$pdHess <- FALSE
  ntries <- 0
  while((!rep$pdHess) & ntries < 10){
    try({
    obj <- MakeADFun(data = TMB_data, parameters = TMB_params, DLL="tmb_MVN_ILR", random = "u_large")
    opt <- do.call("optim", obj)
    opt
    opt$hessian ## <-- FD hessian from optim
    rep <- try(sdreport(obj))
    if(is.character(rep)){
      rep <- list()
      rep$pdHess <- FALSE
    }
    })
    ntries <- ntries+1
  }
  return(rep)
}

imput_vec <- seq(0.001, 0.2, length.out = 20)
run_several_imputation_res <- lapply(imput_vec, run_several_imputation)

plot(sapply(run_several_imputation_res, function(i) i$par.fixed[1]), type='l')
plot(sapply(run_several_imputation_res, function(i) i$par.fixed[2]), type='l')

par(mfrow=c(5,7), mar=c(0,0,0,0))
sapply(1:length(i$par.fixed), function(j){
  plot(sapply(run_several_imputation_res, function(i) i$par.fixed[j]), type='l')
})

pairs(sapply(run_several_imputation_res, function(i) i$par.fixed))

sapply(run_several_imputation_res, wald_TMB_wrapper)
sapply(run_several_imputation_res, wald_TMB_wrapper, fail_non_converged = F)
plot(1:length(run_several_imputation_res), sapply(run_several_imputation_res, wald_TMB_wrapper, fail_non_converged = F), type='l')

repeat_full_with_different_initial_obs = replicate(n = 10, repeat_full_with_different_initial())

saveRDS(repeat_full_with_different_initial_obs, paste0(folder_out_RDS, "repeat_full_with_different_initial_obs.RDS"))
## we have one pdHess true!!
unlist(apply(repeat_full_with_different_initial_obs, 2, `[`, 'pdHess'))

repeat_full_with_different_initial_obs[,9]$pdHess

## Beta values are the same regardless of whether the runs have converged or not
pdf(paste0(folder_images_out, "full_partialILR_several_runs.pdf"))
pairs(repeat_full_with_different_initial_obs['par.fixed',])
dev.off()
pdf(paste0(folder_images_out, "full_partialILR_several_runs_only_beta.pdf"))
pairs(sapply(repeat_full_with_different_initial_obs['par.fixed',], python_like_select_name, grep_substring="beta"),
      main='Only beta (intercept+slope)')
dev.off()

## Analyse the results of the converged run
results_full0 = repeat_full_with_different_initial_obs[,9]
results_full = TMB::summary.sdreport(results_full0)
intercept_est = select_intercept(python_like_select_rownames(results_full, 'beta')[,1]) ## estimate
intercept_err = select_intercept(python_like_select_rownames(results_full, 'beta')[,2]) ## std err
ggplot(cbind.data.frame(name=paste0('ILR', 1:d), intercept_est, intercept_err),
       aes(x=name, y=intercept_est))+
  geom_errorbar(aes(ymin=intercept_est-intercept_err, ymax=intercept_est+intercept_err), width=.1)+
  geom_line() + geom_point()+theme_linedraw()+labs(x='Name of beta intercept', y='Estimate for beta value')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(paste0(folder_images_out, "full_partialILR_intercept_errorbar.png"), width = 3, height = 2.5)

slope_est = select_slope_2(python_like_select_rownames(results_full, 'beta')[,1],v=F) ## estimate
slope_err = select_slope_2(python_like_select_rownames(results_full, 'beta')[,2],v=F) ## std err
ggplot(cbind.data.frame(name=paste0('ILR', 1:d), slope_est, slope_err),
       aes(x=name, y=slope_est))+
  geom_errorbar(aes(ymin=slope_est-slope_err, ymax=slope_est+slope_err), width=.1)+
  geom_line() + geom_point()+theme_linedraw()+labs(x='Name of beta slope', y='Estimate for beta value')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(paste0(folder_images_out, "full_partialILR_slope_errorbar.png"), width = 3, height = 2.5)

wald_TMB_wrapper(results_full0, verbatim = FALSE)


## TRUE_data: subset of signatures. Removing s5
TMB_data_with_subset = prepare_TMB_data_with_subset(c(1:4,6:7))
rownames(TMB_data_with_subset$Y)
TMB_data_with_subset$x
pheatmap(TMB_data_with_subset$Y)
TMB_params_with_subset= give_TMB_params(arg_d = TMB_data_with_subset$d, arg_num_individuals = TMB_data_with_subset$num_individuals)
obj <- MakeADFun(data = TMB_data_with_subset, parameters = TMB_params_with_subset, DLL="tmb_MVN_partial_ILR", random = "u_large")
opt <- do.call("optim", obj)
opt
opt$hessian ## <-- FD hessian from optim
rep <- sdreport(obj)
rep
saveRDS(rep, paste0(folder_out_RDS, "repeat_full_nos5.RDS"))
wald_TMB_wrapper(rep)

ggplot(cbind.data.frame(name=paste0('ILR', c(1:3, 5:6)),
                        intercept_est=select_intercept(python_like_select_rownames(TMB::summary.sdreport(rep), 'beta')[,1]),
                        intercept_err=select_intercept(python_like_select_rownames(TMB::summary.sdreport(rep), 'beta')[,2])),
       aes(x=name, y=intercept_est))+
  geom_errorbar(aes(ymin=intercept_est-intercept_err, ymax=intercept_est+intercept_err), width=.1)+
  geom_line() + geom_point()+theme_linedraw()+labs(x='Name of beta intercept', y='Estimate for beta value')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(paste0(folder_images_out, "full_partialILR_nos5_intercept_errorbar.png"), width = 3, height = 2.5)
ggplot(cbind.data.frame(name=paste0('ILR', c(1:3, 5:6)),
                        slope_est=select_slope_2(python_like_select_rownames(TMB::summary.sdreport(rep), 'beta')[,1],v=F),
                        slope_err=select_slope_2(python_like_select_rownames(TMB::summary.sdreport(rep), 'beta')[,2],v=F)),
       aes(x=name, y=slope_est))+
  geom_errorbar(aes(ymin=slope_est-slope_err, ymax=slope_est+slope_err), width=.1)+
  geom_line() + geom_point()+theme_linedraw()+labs(x='Name of beta slope', y='Estimate for beta value')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(paste0(folder_images_out, "full_partialILR_nos5_slope_errorbar.png"), width = 3, height = 2.5)
#-------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------#
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

## Without random effects; just fixed, and without correlations between categories

TMB::compile("../tmb_RE/mvn_beta_no_cor.cpp", "-std=gnu++17")
dyn.load(dynlib("../tmb_RE/mvn_beta_no_cor"))

TMB_params = give_TMB_params(d, TMB_data$num_individuals)
TMB_params$logs_sd <- TMB_params$logs_sd_RE
obj <- MakeADFun(data = TMB_data, parameters = TMB_params, DLL="mvn_beta_no_cor")
opt <- do.call("optim", obj)
opt
opt$hessian ## <-- FD hessian from optim
rep <- sdreport(obj)
rep
plot_betas(rep)
obj$fn()


TMB_res2 <- wrapper_run_TMB_use_nlminb(give_TMB_data(imputation_val = 1e-3, transformation = 'ALR'), model = "mvn_beta_no_cor")
TMB_res3 <- wrapper_run_TMB_use_nlminb(give_TMB_data(imputation_val = 1e-4, transformation = 'ALR'), model = "mvn_beta_no_cor")
TMB_res4 <- wrapper_run_TMB_use_nlminb(give_TMB_data(imputation_val = 1e-3, transformation = 'ILR'), model = "mvn_beta_no_cor")
TMB_res5 <- wrapper_run_TMB_use_nlminb(give_TMB_data(imputation_val = 'multLN', transformation = 'ALR'), model = "mvn_beta_no_cor")
TMB_res6 <- wrapper_run_TMB_use_nlminb(give_TMB_data(imputation_val = 'multRepl', transformation = 'ALR'), model = "mvn_beta_no_cor")
TMB_res6 <- wrapper_run_TMB_use_nlminb(give_TMB_data(imputation_val = 'imputeBDLs', transformation = 'ALR'), model = "mvn_beta_no_cor")

plot(as.vector(unlist(give_TMB_data(imputation_val = 1e-3,
                                    transformation = 'ALR')$Y)),
     as.vector(unlist(give_TMB_data(imputation_val = 'multRepl',
                                    transformation = 'ALR')$Y)),
     col=rep(c(1:6), each=nrow(exposures)))
abline(coef = c(0,1), lty='dashed')

grid.arrange(plot_betas(TMB_res2, title='Imputation 1e-3, ALR'),
             plot_betas(TMB_res3, title='Imputation 1e-4, ALR'),
             plot_betas(TMB_res4, title='Imputation 1e-3, ILR'),
             plot_betas(TMB_res5, title='multLN imputation, ALR'),
             plot_betas(TMB_res6, title='multRepl imputation, ALR'))


wald_TMB_wrapper(TMB_res2)
wald_TMB_wrapper(TMB_res3)

#-------------------------------------------------------------------------------------------------------------#
## Differential overdispersion


add_accessory_lambda <- function(i){
  i$lambda_accessory_mat = t(sapply(cbind(1, as.numeric(as.factor(patient.meta$group[keep_britroc_samples]))-1)[,2],
                                    function(i) table(factor(i, levels=c(0,1)))))
  i
}

keep_britroc_samples = rep(T, nrow(exposures))
z_britroc = give_z_matrix_from_labels(patient.meta$PATIENT_ID[keep_britroc_samples])
wrapper_run_TMB_use_nlminb(add_accessory_lambda(give_TMB_data(imputation_val = 1e-3, transformation = 'ALR')), model = "tmb_MVN_with_mean_2_scaledsdpergroup")

dataaa <- add_accessory_lambda(give_TMB_data(imputation_val = 1e-3, transformation = 'ALR'))
paraaams <- give_TMB_params(d, dataaa$num_individuals)

paraaams$log_lambda = 0

dataaa$num_individuals
dataaa$d
dataaa$n
sapply(dataaa, dim)
sapply(paraaams, dim)
length(paraaams$logs_sd_RE)

TMB::compile("../tmb_RE/tmb_MVN_with_mean_2_scaledsdpergroup.cpp", "-std=gnu++17")
dyn.load(dynlib("../tmb_RE/tmb_MVN_with_mean_2_scaledsdpergroup"))

obj <- MakeADFun(data = dataaa, parameters = paraaams, DLL="tmb_MVN_with_mean_2_scaledsdpergroup")
opt <- do.call("optim", obj)
opt
opt$hessian ## <-- FD hessian from optim
rep <- sdreport(obj)
rep

rep$par.fixed

#-------------------------------------------------------------------------------------------------------------#

exposures
imputation_methods <- list(list(Y=exposures),
                           give_TMB_data(imputation_val = 1e-2, transformation = 'none'),
                           give_TMB_data(imputation_val = 1e-3, transformation = 'none'),
                           give_TMB_data(imputation_val = 1e-4, transformation = 'none'),
                           give_TMB_data(imputation_val = 'multLN', transformation = 'none'),
                           give_TMB_data(imputation_val = 'multRepl', transformation = 'none'))
                           #give_TMB_data(imputation_val = 'imputeBDLs', transformation = 'none'))

lapply(imputation_methods, '[', 'Y')
lapply(imputation_methods, function(i) min(i$Y))

## find distance between these datasets
## using aitchison's distance for each pair

dists_imputation_methods <- outer(1:length(imputation_methods), 1:length(imputation_methods), Vectorize(function(i,j){
  if(i>j){
    .clri <- as(compositions::clr(imputation_methods[[i]]$Y), 'matrix')
    .clrj <- as(compositions::clr(imputation_methods[[j]]$Y), 'matrix')
    sum(sapply(1:nrow(.clri), function(rw) dist(rbind(.clri[rw,], .clrj[rw,]))))
  }else{
    NA
  }
}))

diag(dists_imputation_methods) = 0
for(i in 1:length(imputation_methods)){
  for(j in 1:length(imputation_methods)){
    if(j>i){
      dists_imputation_methods[i,j] = dists_imputation_methods[j,i]
    }
  }
}
rownames(dists_imputation_methods) <- colnames(dists_imputation_methods) <- c('original', 'imp 1e-2', 'imp 1e-3', 'imp 1e-4',
                                                                              'multLN', 'multRepl')

pheatmap(dists_imputation_methods)

imputation_methods_unlist <- data.frame(sapply(lapply(imputation_methods, '[', 'Y'), unlist))
colnames(imputation_methods_unlist) <- rownames(dists_imputation_methods)
imputation_methods_unlist$colour <- rep(c(1:7), each=nrow(imputation_methods[[1]]$Y))

#pairs(imputation_methods_unlist)

library(GGally)

mycorrelations <- function(data,mapping,...){
  data2 = data
  data2$x = as.numeric(data[,as_label(mapping$x)])
  data2$y = as.numeric(data[,as_label(mapping$y)])
  corstat <-   list(estimate = round(as.numeric(cor.test(data2$x,data2$y,
                                                         method="spearman")$estimate),2),
              pvalue = cor.test(data2$x,data2$y,method="spearman")$p.value)
  corstat$pvalue_star = as.character(symnum(corstat$pvalue, corr = FALSE, na = FALSE,
                                                cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                                                symbols = c("***", "**", "*", "'", " ")))
  
  correlation_df <- data.frame(group=1, x=1, y=1)
  ggplot(data=correlation_df, aes(x=1,y=group,color=group))+
    geom_text(aes(label=paste0('Corr:',": ",corstat$estimate,corstat$pvalue_star)), size=2)
}

imputation_methods_unlist$colour <- as.factor(imputation_methods_unlist$colour)
ggpairs(imputation_methods_unlist,
        mapping = ggplot2::aes(color=(colour), alpha=0.2), 
        upper = list(continuous = mycorrelations), alpha=0.2,
        columns = colnames(imputation_methods_unlist)[-length(colnames(imputation_methods_unlist))])
ggsave(paste0(folder_images_out, "ggpairs_imputation_vals.pdf"), height=5, width = 5)

#-------------------------------------------------------------------------------------------------------------#
## check that ILR is multivariate-normally distributed
## MVN distrib
exposures

## check function <give_pairs_with_mvn_wrapper> which generalises the plots below

exposures_imput1em2_ilr <- as(compositions::ilr(exposures_imput1em2), 'matrix')
exposures_imput1em2_alr <- as(compositions::alr(exposures_imput1em2), 'matrix')
exposures_partialILR <- give_partial_irl(exposures)
exposures_partialILR[exposures_partialILR == 0] <- NA

pdf(paste0(folder_images_out, 'normality_assumption_pairs_ILR_imput1em2.pdf'))
par(mfrow=c(ncol(exposures_imput1em2_ilr), ncol(exposures_imput1em2_ilr)), mar=c(0,0,0,0))
for(ii in 1:ncol(exposures_imput1em2_ilr)){
  for(jj in 1:ncol(exposures_imput1em2_ilr)){
    if(ii != jj){
      give_pairs_with_mvn(exposures_imput1em2_ilr[,c(ii, jj)])
    }else{plot.new()}
  }
}
dev.off()

pdf(paste0(folder_images_out, 'normality_assumption_pairs_ALR_imput1em2.pdf'))
par(mfrow=c(ncol(exposures_imput1em2_alr), ncol(exposures_imput1em2_alr)), mar=c(0,0,0,0))
for(ii in 1:ncol(exposures_imput1em2_alr)){
  for(jj in 1:ncol(exposures_imput1em2_alr)){
    if(ii != jj){
      give_pairs_with_mvn(exposures_imput1em2_alr[,c(ii, jj)])
    }else{plot.new()}
  }
}
dev.off()

pdf(paste0(folder_images_out, 'normality_assumption_pairs_partialILR.pdf'))
par(mfrow=c(ncol(exposures_partialILR), ncol(exposures_partialILR)), mar=c(0,0,0,0))
for(ii in 1:ncol(exposures_partialILR)){
  for(jj in 1:ncol(exposures_partialILR)){
    if(ii != jj){
      give_pairs_with_mvn(remove_all_NA(exposures_partialILR[,c(ii, jj)]))
    }else{plot.new()}
  }
}
dev.off()

