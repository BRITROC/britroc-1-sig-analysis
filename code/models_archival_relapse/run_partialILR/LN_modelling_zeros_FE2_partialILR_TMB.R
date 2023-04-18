## FE using mixed effects code

rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# set.seed(1234)

library(uuid)
library(ggplot2)
library(reshape2)
library(compositions)
library(TMB)
library(ComplexHeatmap)
library(gridExtra)

# source("../../../GlobalDA/code/2_inference/helper/helper_DA_stan.R")
# source("../../../GlobalDA/code/3_analysis/helper/helper_analyse_posteriors.R")
source("../../../../../GlobalDA/code/2_inference_TMB/helper_TMB.R")
source("../helper/functions.R")
source("../helper/header.R")

#-------------------------------------------------------------------------------------------#

all(patient.meta$SAMPLE_ID == rownames(exposures))

folder_out_RDS <- "../../../out/inference/partialILR_FE2/"
folder_images_out <- "../../../results/partialILRmodelling_FE2/"

system(paste0("mkdir -p ", folder_out_RDS))
system(paste0("mkdir -p ", folder_images_out))

#-------------------------------------------------------------------------------------------#

## Transform data with partial ILR
which_zero = t(apply(exposures, 1, function(i) as.numeric((i==0)) ))


exposures[1,]
which_zero[1,]
as.vector(compositions::ilr(exposures)[1,])
give_partial_ilr_basis(which(which_zero[1,] == 1), d=7)
give_partial_ilr_basis(which_zero_vector = 2, d = 5)


give_partial_ilr_basis(which_zero_vector = 2, d = 4)

prepare_TMB_data_with_subset = function(subset_sigs, additional_keep_bool=NULL){
  .exps = normalise_rw(exposures[,subset_sigs])
  .irl_with_zeros = give_partial_irl(.exps)
  if(is.null(additional_keep_bool)){
    additional_keep_bool = rep(T, nrow(.exps))
  }
  .keep = (rowSums(.irl_with_zeros == 0) < (ncol(.exps) - 2) ) & additional_keep_bool
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

exposures_partial_irl = give_partial_irl(exposures)

dim(exposures)
dim(exposures_partial_irl)

exposures_partial_irl_with_NAzeros <- exposures_partial_irl
exposures_partial_irl_with_NAzeros[exposures_partial_irl_with_NAzeros == 0] <- NA
pairs(exposures_partial_irl_with_NAzeros)

ComplexHeatmap::Heatmap(exposures_partial_irl_with_NAzeros, cluster_rows = F, cluster_columns = F)

## note that the partial ILR that I have  computed is different from the ILR setting zero entries as zero ILRs
## different
# compositions::clr(exposures[1,]) %*% give_partial_ilr_basis(which_zero[[1]])
# compositions::clr(exposures[1,]) %*% irl_base_complete

## i.e the zero entries are different

# To show that I am computing the ILR correctly given a basis
# plot(as.vector(compositions::ilr(exposures)[1,]),
#      compositions::clr(exposures[1,]) %*% irl_base_complete)
# abline(coef = c(0,1))
# 
# #-------------------------------------------------------------------------------------------#
# TMB::compile("../tmb_RE/tmb_MVN_partial_ILR.cpp", "-std=gnu++17")
# dyn.load(dynlib("../tmb_RE/tmb_MVN_partial_ILR"))
# # TMB::compile("../tmb_RE/tmb_MVN_ILR.cpp", "-std=gnu++17")
# # dyn.load(dynlib("../tmb_RE/tmb_MVN_ILR"))
# TMB::compile("../tmb_RE/tmb_MVN_partial_ILR_noularge.cpp", "-std=gnu++17")
# dyn.load(dynlib("../tmb_RE/tmb_MVN_partial_ILR_noularge"))
# TMB::compile("../tmb_RE/tmb_MVN_partial_ILR_FEc.cpp", "-std=gnu++17")
# dyn.load(dynlib("../tmb_RE/tmb_MVN_partial_ILR_FEc"))
# #-------------------------------------------------------------------------------------------#
# 
# 
# #-------------------------------------------------------------------------------------------#
# ## remove the samples in which there are too many zeros. There can be at most 4 zeros (i.e. there has to be at least one log-ratio)
# keep_britroc_samples = (rowSums(exposures_partial_irl == 0) < (ncol(exposures) - 2) )
# table(keep_britroc_samples)
# exposures_partial_irl = exposures_partial_irl[keep_britroc_samples,]
# all(patient.meta$SAMPLE_ID[keep_britroc_samples] == rownames(exposures_partial_irl))
# z_britroc = give_z_matrix_from_labels(1:sum(keep_britroc_samples))
# dim(z_britroc)
# TMB_data = list(Y = exposures_partial_irl,
#                 d = d,
#                 n = nrow(exposures_partial_irl),
#                 x = cbind(1, as.numeric(as.factor(patient.meta$group[keep_britroc_samples]))-1))
# 
# sapply(TMB_data, dim)
# sapply(TMB_params, dim)
# #-------------------------------------------------------------------------------------------#
# 
# #-------------------------------------------------------------------------------------------#
# ## TRUE_data: full
# pheatmap(TMB_data$Y, cluster_cols = F)
# TMB_params = give_TMB_params(d, arg_num_individuals = 1) ## dummy arg_num_individuals
# TMB_params$u_large <- NULL
# obj <- MakeADFun(data = TMB_data, parameters = TMB_params, DLL="tmb_MVN_partial_ILR_noularge")
# opt <- do.call("optim", obj)
# opt
# opt$hessian ## <-- FD hessian from optim
# rep <- sdreport(obj)
# rep
# 
# repeat_full_with_different_initial = function(){
#   obj <- MakeADFun(data = TMB_data, parameters = TMB_params, DLL="tmb_MVN_partial_ILR_noularge")
#   opt <- do.call("optim", obj)
#   opt
#   opt$hessian ## <-- FD hessian from optim
#   rep <- sdreport(obj)
#   return(rep)
# }
# 
# repeat_full_with_different_initial_obs = replicate(n = 20, tryCatch(repeat_full_with_different_initial()))
# 
# # saveRDS(repeat_full_with_different_initial_obs, paste0(folder_out_RDS, "repeat_full_with_different_initial_obs.RDS"))
# ## we have one pdHess true!!
# 
# unlist(sapply(repeat_full_with_different_initial_obs, function(i) try(i$pdHess) ))
# unlist(apply(repeat_full_with_different_initial_obs, 2, function(i) try(i$pdHess) ))
# 
# repeat_full_wit_different_initial_obs[,1]
# 
# repeat_full_with_different_initial_obs <- repeat_full_with_different_initial_obs[sapply(repeat_full_with_different_initial_obs, typeof) == "list"]
# repeat_full_with_different_initial_obs[,9]$pdHess
# sapply(repeat_full_with_different_initial_obs, wald_TMB_wrapper, verbatim = FALSE)
# 
# plot_betas(TMB::summary.sdreport(repeat_full_with_different_initial_obs[,1]))
# plot_betas(repeat_full_with_different_initial_obs[,2])
# 
# ## Beta values are the same regardless of whether the runs have converged or not
# # pdf(paste0(folder_images_out, "full_partialILR_several_runs.pdf"))
# # pairs(repeat_full_with_different_initial_obs['par.fixed',])
# # dev.off()
# # pdf(paste0(folder_images_out, "full_partialILR_several_runs_only_beta.pdf"))
# # pairs(sapply(repeat_full_with_different_initial_obs['par.fixed',], python_like_select_name, grep_substring="beta"),
# #       main='Only beta (intercept+slope)')
# # dev.off()
# 
# which(unlist(sapply(repeat_full_with_different_initial_obs, function(i) try(i$pdHess) )))
# 
# 
# ## Analyse the results of the converged run
# results_full0 = repeat_full_with_different_initial_obs[[5]]
# results_full = TMB::summary.sdreport(results_full0)
# intercept_est = select_intercept(python_like_select_rownames(results_full, 'beta')[,1]) ## estimate
# intercept_err = select_intercept(python_like_select_rownames(results_full, 'beta')[,2]) ## std err
# ggplot(cbind.data.frame(name=paste0('ILR', 1:d), intercept_est, intercept_err),
#        aes(x=name, y=intercept_est))+
#   geom_errorbar(aes(ymin=intercept_est-intercept_err, ymax=intercept_est+intercept_err), width=.1)+
#   geom_line() + geom_point()+theme_linedraw()+labs(x='Name of beta intercept', y='Estimate for beta value')+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
# # ggsave(paste0(folder_images_out, "full_partialILR_intercept_errorbar.png"), width = 3, height = 2.5)
# 
# slope_est = select_slope_2(python_like_select_rownames(results_full, 'beta')[,1],v=F) ## estimate
# slope_err = select_slope_2(python_like_select_rownames(results_full, 'beta')[,2],v=F) ## std err
# ggplot(cbind.data.frame(name=paste0('ILR', 1:d), slope_est, slope_err),
#        aes(x=name, y=slope_est))+
#   geom_errorbar(aes(ymin=slope_est-slope_err, ymax=slope_est+slope_err), width=.1)+
#   geom_line() + geom_point()+theme_linedraw()+labs(x='Name of beta slope', y='Estimate for beta value')+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
# # ggsave(paste0(folder_images_out, "full_partialILR_slope_errorbar.png"), width = 3, height = 2.5)
# 
# wald_TMB_wrapper(results_full0, verbatim = FALSE)
# 
# ##--------------------------------------------------------------------------------------------##
# ## TRUE_data: subset of signatures. Removing s5
# # TMB_data_with_subset = prepare_TMB_data_with_subset(c(1:4,6:7), additional_keep_bool = (1:nrow(exposures) %in% 1:200)) #rep(T, nrow(exposures)))
# TMB_data_with_subset = prepare_TMB_data_with_subset(c( 3, 6, 7), additional_keep_bool = rep(T, nrow(exposures)))
# pairs(TMB_data_with_subset$Y)
# TMB_data_with_subset$z <- NULL
# TMB_data_with_subset$num_individuals <- NULL
# rownames(TMB_data_with_subset$Y)
# TMB_data_with_subset$x
# pheatmap(TMB_data_with_subset$Y)
# # TMB_params_with_subset= give_TMB_params(arg_d = TMB_data_with_subset$d, arg_num_individuals = 1)
# # TMB_params_with_subset$u_large <- NULL
# # obj <- MakeADFun(data = TMB_data_with_subset, parameters = TMB_params_with_subset, DLL="tmb_MVN_partial_ILR_noularge")
# # obj$maxit <- 10000
# # opt <- do.call("optim",obj)
# # opt
# # opt$hessian ## <-- FD hessian from optim
# # rep <- sdreport(obj)
# # rep
# # # saveRDS(rep, paste0(folder_out_RDS, "repeat_full_nos5.RDS"))
# # wald_TMB_wrapper(rep)
# # wald_TMB_wrapper(rep, fail_non_converged=F)
# # rep$cov.fixed
# # i=rep
# # idx_beta = select_slope_2(which(names(i$par.fixed) == "beta"), verbatim=verbatim)
# # wald_generalised(v = i$par.fixed[idx_beta], sigma = i$cov.fixed[idx_beta,idx_beta])
# # i$par.fixed[idx_beta]
# # plot_betas(rep)
# pheatmap(i$cov.fixed[idx_beta,idx_beta], cluster_rows = F, cluster_cols = F)
# 
# # TMB_data_with_subset_for_nlminb <- prepare_TMB_data_with_subset(c(1, 3, 6, 7), additional_keep_bool = rep(T, nrow(exposures)))
# # TMB_data_with_subset_for_nlminb$z <- diag(nrow(TMB_data_with_subset_for_nlminb$z))
# # TMB_data_with_subset_for_nlminb$d = TMB_data_with_subset_for_nlminb$d+1
# # # res_nlminb <- wrapper_run_TMB_use_nlminb(model = "tmb_MVN_partial_ILR_noularge", object = TMB_data_with_subset_for_nlminb,
# # #                                          use_nlminb = T)
# # # res_nlminb$pdHess
# # ## this model doesn't make any sense
# # # pheatmap(res_nlminb$cov.fixed, cluster_rows = F, cluster_cols = F)
# # # pheatmap(res_nlminb$cov.fixed[idx_beta,idx_beta], cluster_rows = F, cluster_cols = F)
# # 
# # TMB::compile("../tmb_RE/tmb_MVN_partial_ILR_FEc.cpp", "-std=gnu++17")
# # dyn.load(dynlib("../tmb_RE/tmb_MVN_partial_ILR_FEc"))
# # 
# # res_nlminb <- wrapper_run_TMB_use_nlminb(model = "tmb_MVN_partial_ILR_FEc", object = TMB_data_with_subset_for_nlminb,
# #                                          use_nlminb = T)
# # res_nlminb$pdHess
# # pheatmap(res_nlminb$cov.fixed, cluster_rows = F, cluster_cols = F)
# 
# ggplot(cbind.data.frame(name=paste0('ILR', c(1:3, 5:6)),
#                         intercept_est=select_intercept(python_like_select_rownames(TMB::summary.sdreport(rep), 'beta')[,1]),
#                         intercept_err=select_intercept(python_like_select_rownames(TMB::summary.sdreport(rep), 'beta')[,2])),
#        aes(x=name, y=intercept_est))+
#   geom_errorbar(aes(ymin=intercept_est-intercept_err, ymax=intercept_est+intercept_err), width=.1)+
#   geom_line() + geom_point()+theme_linedraw()+labs(x='Name of beta intercept', y='Estimate for beta value')+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
# # ggsave(paste0(folder_images_out, "full_partialILR_nos5_intercept_errorbar.png"), width = 3, height = 2.5)
# ggplot(cbind.data.frame(name=paste0('ILR', c(1:3, 5:6)),
#                         slope_est=select_slope_2(python_like_select_rownames(TMB::summary.sdreport(rep), 'beta')[,1],v=F),
#                         slope_err=select_slope_2(python_like_select_rownames(TMB::summary.sdreport(rep), 'beta')[,2],v=F)),
#        aes(x=name, y=slope_est))+
#   geom_errorbar(aes(ymin=slope_est-slope_err, ymax=slope_est+slope_err), width=.1)+
#   geom_line() + geom_point()+theme_linedraw()+labs(x='Name of beta slope', y='Estimate for beta value')+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
# # ggsave(paste0(folder_images_out, "full_partialILR_nos5_slope_errorbar.png"), width = 3, height = 2.5)
# #-------------------------------------------------------------------------------------------#
# 
# #-------------------------------------------------------------------------------------------#
# 
# normalise_rw(exposures[,4:7])[patient.meta$PATIENT_ID]
# 
# all(rownames(exposures) == patient.meta$SAMPLE_ID)
# 
# plot_ternary(exposures[patient.meta$group == "rlps",c(1,3,5)])
# plot_ternary(exposures[patient.meta$group == "arx",c(1,3,5)])
# 
# par(mfrow=c(1,2))
# plot_ternary(remove_only_zero_rows(exposures[patient.meta$group == "arx",c(2,3,7)]))
# plot_ternary(remove_only_zero_rows(exposures[patient.meta$group == "rlps",c(2,3,7)]))
# 
# ## plot by patients which have no relapse sample
# 
# dev.off()
# pdf(paste0(folder_images_out, "ternary_plot_samples_archival_relapse.pdf"), width=8, height=3)
# par(mfrow=c(1,3))
# ## patients with at least one relapse sample
# ## patients with no relapse sample
# plot_ternary(remove_only_zero_rows(exposures[!(patient.meta$PATIENT_ID %in% patient.meta$PATIENT_ID[patient.meta$group == "rlps"]),
#                                              c(2,3,7)]), main = "Patients without relapse sample.\n Archival samples")
# plot_ternary(remove_only_zero_rows(exposures[(patient.meta$PATIENT_ID %in% patient.meta$PATIENT_ID[patient.meta$group == "rlps"] &
#                                                 patient.meta$group == "arx"),
#                                              c(2,3,7)]), main = "Patients with relapse sample.\n Archival samples.")
# plot_ternary(remove_only_zero_rows(exposures[(patient.meta$PATIENT_ID %in% patient.meta$PATIENT_ID[patient.meta$group == "rlps"] &
#                                                 patient.meta$group == "rlps"),
#                                              c(2,3,7)]), main = "Patients with relapse sample.\n Relapse samples.")
# dev.off()
# 
# pdf(paste0(folder_images_out, "ternary_plot_samples_archival_relapse_2.pdf"), width=8, height=3)
# par(mfrow=c(1,3))
# ## patients with at least one relapse sample
# ## patients with no relapse sample
# plot_ternary(remove_only_zero_rows(exposures[!(patient.meta$PATIENT_ID %in% patient.meta$PATIENT_ID[patient.meta$group == "rlps"]),
#                                              c(4,3,7)]), main = "Patients without relapse sample.\n Archival samples")
# plot_ternary(remove_only_zero_rows(exposures[(patient.meta$PATIENT_ID %in% patient.meta$PATIENT_ID[patient.meta$group == "rlps"] &
#                                                 patient.meta$group == "arx"),
#                                              c(4,3,7)]), main = "Patients with relapse sample.\n Archival samples.")
# plot_ternary(remove_only_zero_rows(exposures[(patient.meta$PATIENT_ID %in% patient.meta$PATIENT_ID[patient.meta$group == "rlps"] &
#                                                 patient.meta$group == "rlps"),
#                                              c(4,3,7)]), main = "Patients with relapse sample.\n Relapse samples.")
# dev.off()
# 
# pdf(paste0(folder_images_out, "ternary_plot_samples_archival_relapse_PhilAmplificationSamples.pdf"), width=6, height=5)
# par(mfrow=c(2,2), mar=c(0.2,0.2,0.8,0.2))
# subset_sigs <- c(4,3,7)
# plot_ternary(remove_only_zero_rows(exposures[patient.meta$group == "arx",subset_sigs]), main = "All patients, arx")
# plot_ternary(remove_only_zero_rows(exposures[patient.meta$group == "rlps", subset_sigs]), main = "All patients, rel")
# plot_ternary(remove_only_zero_rows(exposures[(patient.meta$PATIENT_ID %in% paste0('BRITROC-', c(65, 241, 23, 267, 216, 74, 274, 209)))&
#                                                (patient.meta$group == "arx"),
#                                              subset_sigs]), main = "Interesting patients, arx")
# plot_ternary(remove_only_zero_rows(exposures[(patient.meta$PATIENT_ID %in% paste0('BRITROC-', c(65, 241, 23, 267, 216, 74, 274, 209)))&
#                                                (patient.meta$group == "rlps"),
#                                              subset_sigs]), main = "Interesting patients, rel")
# dev.off()
# 
# ###--- with imputation
# impute_value <- 0.01
# exposures_partial_irl_imput = give_partial_irl(exposures+impute_value)
# 
# keep_britroc_samples_imput = (rowSums(exposures_partial_irl == 0) < (ncol(exposures+impute_value) - 2) )
# exposures_partial_irl_imput = exposures_partial_irl_imput[keep_britroc_samples_imput,]
# patient.meta_imput <- patient.meta
# patient.meta_imput <- patient.meta_imput[match(rownames(exposures_partial_irl_imput), patient.meta_imput$SAMPLE_ID),]
# 
# all(patient.meta_imput$SAMPLE_ID[keep_britroc_samples_imput] == rownames(exposures_partial_irl_imput))
# z_britroc_imput = give_z_matrix_from_labels(1:sum(keep_britroc_samples_imput))
# TMB_data_imput = list(Y = exposures_partial_irl_imput,
#                 num_individuals = ncol(z_britroc_imput),
#                 d = d,
#                 n = nrow(exposures_partial_irl_imput),
#                 x = cbind(1, as.numeric(as.factor(patient.meta_imput$group))-1),
#                 z = z_britroc_imput)
# 
# ## TRUE_data: full
# pheatmap(TMB_data_imput$Y, cluster_cols = F)
# TMB_params_imput = give_TMB_params(d, TMB_data_imput$num_individuals)
# obj_imput <- MakeADFun(data = TMB_data_imput, parameters = TMB_params_imput, DLL="tmb_MVN_partial_ILR")
# obj_imput <- MakeADFun(data = TMB_data_imput, parameters = TMB_params_imput, DLL="tmb_MVN_partial_ILR")
# opt_imput <- do.call("optim", obj_imput)
# opt_imput
# opt_imput$hessian ## <-- FD hessian from optim
# rep_imput <- sdreport(obj_imput)
# rep_imput
# 
# pheatmap(TMB_data$Y, cluster_cols = F)
# pheatmap(TMB_data$Y, cluster_cols = F)
# 
# pairs(TMB_data$Y) 
# pairs(TMB_data_imput$Y) ## interesting patters of what might be samples from the same patient, or due to zero artefacts, or imputation
# 
# GGally::ggpairs(data.frame(TMB_data_imput$Y, patient=patient.meta_imput$PATIENT_ID[keep_britroc_samples_imput]),
#                 aes(col=patient))
# 
# obj_imput <- MakeADFun(data = TMB_data_imput, parameters = TMB_params_imput, DLL="tmb_MVN_ILR")
# obj_imput <- MakeADFun(data = TMB_data_imput, parameters = TMB_params_imput, DLL="tmb_MVN_ILR")
# opt_imput <- do.call("optim", obj_imput)
# opt_imput
# opt_imput$hessian ## <-- FD hessian from optim
# rep_imput <- sdreport(obj_imput)
# rep_imput
# 
# plot_betas(rep_imput) 
# 
# colSums(exposures_zeros[x[2,] == 0,])/nrow(exposures_zeros)
# colSums(exposures_zeros[x[2,] == 1,])/nrow(exposures_zeros)
# 
# par(mfrow=c(1,1))
# hclust_britroc <- hclust(dist(exposures))
# plot(hclust_britroc)
# 
# source("../../../../Vias_Brenton/copy_number_analysis_organoids/helper_functions.R")
# library(ggdendro)
# library(ggrepel)
# library(gridExtra)
# extra_expand <- 0.2
# extra_expand_v2 <- 0.2
# give_dendrogram_from_imputation(impute_VALUE=0.01, plot=T, exposures=exposures)
# give_dendrogram_from_imputation(impute_VALUE=0.03, plot=T, exposures=exposures)
# 

##--------
TMB_data_with_subset_for_nlminb <- prepare_TMB_data_with_subset(c(1:4, 6, 7), additional_keep_bool = rep(T, nrow(exposures)))
TMB_data_with_subset_for_nlminb$z <- diag(nrow(TMB_data_with_subset_for_nlminb$z)) ## need to keep this for size of ularge
# TMB_data_with_subset_for_nlminb$z <- NULL
TMB_data_with_subset_for_nlminb$num_individuals <- NULL
# TMB_data_with_subset_for_nlminb$d = TMB_data_with_subset_for_nlminb$d+1

TMB::compile("../tmb_RE/tmb_MVN_partial_ILR_FEd.cpp", "-std=gnu++17")
dyn.load(dynlib("../tmb_RE/tmb_MVN_partial_ILR_FEd"))

# res_nlminb <- wrapper_run_TMB_use_nlminb(model = "tmb_MVN_partial_ILR_FEd", object = TMB_data_with_subset_for_nlminb,
#                                          use_nlminb = T)
res_nlminb_Fed_integrateularge <- wrapper_run_TMB_use_nlminb(model = "tmb_MVN_partial_ILR_FEd_integrateularge",
                                                             object = TMB_data_with_subset_for_nlminb,
                                                             use_nlminb = T)

res_nlminb_Fed_integrateularge$pdHess
# res_nlminb$pdHess
TMB::sdreport(res_nlminb)
wald_TMB_wrapper(res_nlminb)
wald_TMB_wrapper(res_nlminb, fail_non_converged = F)
wald_TMB_wrapper(res_nlminb_Fed_integrateularge, fail_non_converged = F)
wald_TMB_wrapper(res_nlminb_Fed_integrateularge, fail_non_converged = T)
saveRDS(res_nlminb_Fed_integrateularge, paste0(folder_out_RDS, "partialILR_FEd_nos5_nlminb.RDS"))

pheatmap(res_nlminb$cov.fixed, cluster_rows = F, cluster_cols = F)

grid.arrange(plot_betas(res_nlminb),
plot_betas(res_nlminb_Fed_integrateularge)) ## why so different?

ComplexHeatmap::Heatmap(res_nlminb_Fed_integrateularge$cov.fixed)
ComplexHeatmap::Heatmap(res_nlminb$cov.fixed[!grepl('u_large', rownames(res_nlminb$cov.fixed)),
                                             !grepl('u_large', colnames(res_nlminb$cov.fixed))])
res_nlminb_Fed_integrateularge_summary = TMB::summary.sdreport(res_nlminb_Fed_integrateularge)

intercept_est_partialILR_FEdint = select_intercept(python_like_select_rownames(res_nlminb_Fed_integrateularge_summary, 'beta')[,1]) ## estimate
intercept_err_partialILR_FEdint = select_intercept(python_like_select_rownames(res_nlminb_Fed_integrateularge_summary, 'beta')[,2]) ## std err
ggplot(cbind.data.frame(name=paste0('ILR', c(1:d)[-5]), intercept_est_partialILR_FEdint, intercept_err_partialILR_FEdint),
       aes(x=name, y=intercept_est_partialILR_FEdint))+
  geom_errorbar(aes(ymin=intercept_est_partialILR_FEdint-intercept_err_partialILR_FEdint,
                    ymax=intercept_est_partialILR_FEdint+intercept_err_partialILR_FEdint), width=.1)+
  geom_line() + geom_point()+theme_linedraw()+labs(x='Name of beta intercept', y='Estimate for beta value')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(paste0(folder_images_out, "full_partialILR_FEdint_intercept_errorbar.png"), width = 3, height = 2.5)

slope_est_partialILR_FEdint = select_slope_2(python_like_select_rownames(res_nlminb_Fed_integrateularge_summary, 'beta')[,1],v=F) ## estimate
slope_err_partialILR_FEdint = select_slope_2(python_like_select_rownames(res_nlminb_Fed_integrateularge_summary, 'beta')[,2],v=F) ## std err
ggplot(cbind.data.frame(name=paste0('ILR', c(1:d)[-5]), slope_est_partialILR_FEdint, slope_err_partialILR_FEdint),
       aes(x=name, y=slope_est_partialILR_FEdint))+
  geom_errorbar(aes(ymin=slope_est_partialILR_FEdint-slope_err_partialILR_FEdint,
                    ymax=slope_est_partialILR_FEdint+slope_err_partialILR_FEdint), width=.1)+
  geom_line() + geom_point()+theme_linedraw()+labs(x='Name of beta slope', y='Estimate for beta value')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(paste0(folder_images_out, "full_partialILR_FEdint_slope_errorbar.png"), width = 3, height = 2.5)

