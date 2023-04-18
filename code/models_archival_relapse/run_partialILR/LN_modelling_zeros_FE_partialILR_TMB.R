# 
# rm(list = ls())
# setwd(dirname(rstudioapi::getSourceEditorContext()$path))
# 
# # set.seed(1234)
# 
# library(uuid)
# library(ggplot2)
# library(reshape2)
# library(compositions)
# library(TMB)
# # source("../../../GlobalDA/code/2_inference/helper/helper_DA_stan.R")
# # source("../../../GlobalDA/code/3_analysis/helper/helper_analyse_posteriors.R")
# source("../../../../../GlobalDA/code/2_inference_TMB/helper_TMB.R")
# source("../helper/functions.R")
# source("../helper/header.R")
# 
# folder_out_RDS <- "../../../out/inference/partialILR_FE/"
# folder_images_out <- "../../../results/partialILRmodelling_FE/"
# system(paste0("mkdir -p ", folder_out_RDS))
# system(paste0("mkdir -p ", folder_images_out))
# 
# #-------------------------------------------------------------------------------------------#
# 
# ## Transform data with partial ILR
# which_zero = t(apply(exposures, 1, function(i) as.numeric((i==0)) ))
# 
# 
# exposures[1,]
# which_zero[1,]
# as.vector(compositions::ilr(exposures)[1,])
# give_partial_ilr_basis(which(which_zero[1,] == 1), d=7)
# give_partial_ilr_basis(which_zero_vector = 2, d = 5)
# 
# 
# give_partial_ilr_basis(which_zero_vector = 2, d = 4)
# 
# tmb_MVN_partial_ILR_FEb = function(subset_sigs, .keep_additional){
#   .exps = normalise_rw(exposures[,subset_sigs])
#   .irl_with_zeros = give_partial_irl(.exps)
#   .keep = ((rowSums(.irl_with_zeros == 0) < (ncol(.exps) - 2) )) & .keep_additional
#   .z_britroc = give_z_matrix_from_labels(patient.meta$PATIENT_ID[.keep])
#   .irl_with_zeros = .irl_with_zeros[.keep,]  ## if a sample has too many zeros (i.e. all but 1 or 2) remove
# 
#   TMB_data_sim = list(Y = .irl_with_zeros,
#                       num_individuals = ncol(.z_britroc),
#                       d = length(subset_sigs)-1,
#                       n = nrow(.exps),
#                       x = cbind(1, as.numeric(as.factor(patient.meta$group[.keep]))-1),
#                       z = .z_britroc)
#   return(TMB_data_sim)
# }
# 
# exposures_partial_irl = give_partial_irl(exposures)
# 
# dim(exposures)
# dim(exposures_partial_irl)
# 
# #-------------------------------------------------------------------------------------------#
# 
# #-------------------------------------------------------------------------------------------#
# ## We are only using non-correlated MVN (tmb_MVN_partial_ILR_FEb.cpp). There is no correlated
# ## version as it would be quite tricky (we would have to subset the covariance matrix in each
# ## observation, and I haven't implemented it. The sparsecov versions of the DM models are also
# ## simpler in that there is one single subset)
# TMB::compile("../tmb_RE/tmb_MVN_partial_ILR_FEb.cpp", "-std=gnu++17")
# dyn.load(dynlib("../tmb_RE/tmb_MVN_partial_ILR_FEb"))
# # TMB::compile("tmb_RE/tmb_MVN_partial_ILR.cpp", "-std=gnu++17")
# # dyn.load(dynlib("tmb_RE/tmb_MVN_partial_ILR"))
# #-------------------------------------------------------------------------------------------#
# 
# 
# #-------------------------------------------------------------------------------------------#
# ## remove the samples in which there are too many zeros. There can be at most 4 zeros (i.e. there has to be at least one log-ratio)
# keep_britroc_samples = (rowSums(exposures_partial_irl == 0) < (ncol(exposures) - 2) )
# exposures_partial_irl = exposures_partial_irl[keep_britroc_samples,]
# z_britroc = give_z_matrix_from_labels(1:sum(keep_britroc_samples))
# TMB_data = list(Y = exposures_partial_irl,
#                 num_individuals = ncol(z_britroc),
#                 d = d,
#                 n = nrow(exposures_partial_irl),
#                 x = cbind(1, as.numeric(as.factor(patient.meta$group[keep_britroc_samples]))-1),
#                 z = z_britroc)
# 
# sapply(TMB_data, dim)
# #-------------------------------------------------------------------------------------------#
# 
# #-------------------------------------------------------------------------------------------#
# ## TRUE_data: full
# TMB_params = give_TMB_params(d, TMB_data$num_individuals)
# TMB_params$logsd = rep(1, d)
# TMB_params$logs_sd_RE <- NULL
# TMB_params$cov_RE <- NULL
# sapply(TMB_params, dim)
# obj <- MakeADFun(data = TMB_data, parameters = TMB_params, DLL="tmb_MVN_partial_ILR_FEb")
# opt <- do.call("optim", obj)
# opt
# opt$hessian ## <-- FD hessian from optim
# rep <- sdreport(obj)
# rep
# 
# ## Analyse the results of the converged run
# results_full0 = rep
# results_full = TMB::summary.sdreport(results_full0)
# intercept_est = select_intercept(python_like_select_rownames(results_full, 'beta')[,1]) ## estimate
# intercept_err = select_intercept(python_like_select_rownames(results_full, 'beta')[,2]) ## std err
# 
# plot_betas(results_full0)
# 
# ggplot(cbind.data.frame(name=paste0('ILR', 1:d), intercept_est, intercept_err),
#        aes(x=name, y=intercept_est))+
#   geom_errorbar(aes(ymin=intercept_est-intercept_err, ymax=intercept_est+intercept_err), width=.1)+
#   geom_line() + geom_point()+theme_linedraw()+labs(x='Name of beta intercept', y='Estimate for beta value')+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
# ggsave(paste0(folder_images_out, "full_partialILR_intercept_errorbar.png"), width = 3, height = 2.5)
# 
# slope_est = select_slope_2(python_like_select_rownames(results_full, 'beta')[,1],v=F) ## estimate
# slope_err = select_slope_2(python_like_select_rownames(results_full, 'beta')[,2],v=F) ## std err
# ggplot(cbind.data.frame(name=paste0('ILR', 1:d), slope_est, slope_err),
#        aes(x=name, y=slope_est))+
#   geom_errorbar(aes(ymin=slope_est-slope_err, ymax=slope_est+slope_err), width=.1)+
#   geom_line() + geom_point()+theme_linedraw()+labs(x='Name of beta slope', y='Estimate for beta value')+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
# ggsave(paste0(folder_images_out, "full_partialILR_slope_errorbar.png"), width = 3, height = 2.5)
# 
# wald_TMB_wrapper(results_full0, verbatim = FALSE)
# 
# 
# ## TRUE_data: subset of signatures. Removing s5
# TMB_data_with_subset = tmb_MVN_partial_ILR_FEb(c(1:4,6:7))
# pheatmap(TMB_data_with_subset$Y)
# TMB_params_with_subset= give_TMB_params(arg_d = TMB_data_with_subset$d, arg_num_individuals = TMB_data_with_subset$num_individuals)
# TMB_params_with_subset$logsd = rep(1, d-1)
# obj <- MakeADFun(data = TMB_data_with_subset, parameters = TMB_params_with_subset, DLL="tmb_MVN_partial_ILR_FEb")
# opt <- do.call("optim", obj)
# opt
# opt$hessian ## <-- FD hessian from optim
# rep <- sdreport(obj)
# rep
# # saveRDS(rep, "../out/inference/partialILR/repeat_full_nos5.RDS")
# wald_TMB_wrapper(rep)
# 
# ggplot(cbind.data.frame(name=paste0('ILR', c(1:3, 5:6)),
#                         intercept_est=select_intercept(python_like_select_rownames(TMB::summary.sdreport(rep), 'beta')[,1]),
#                         intercept_err=select_intercept(python_like_select_rownames(TMB::summary.sdreport(rep), 'beta')[,2])),
#        aes(x=name, y=intercept_est))+
#   geom_errorbar(aes(ymin=intercept_est-intercept_err, ymax=intercept_est+intercept_err), width=.1)+
#   geom_line() + geom_point()+theme_linedraw()+labs(x='Name of beta intercept', y='Estimate for beta value')+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
# ggsave(paste0(folder_images_out, "full_partialILR_nos5_intercept_errorbar.png"), width = 3, height = 2.5)
# ggplot(cbind.data.frame(name=paste0('ILR', c(1:3, 5:6)),
#                         slope_est=select_slope_2(python_like_select_rownames(TMB::summary.sdreport(rep), 'beta')[,1],v=F),
#                         slope_err=select_slope_2(python_like_select_rownames(TMB::summary.sdreport(rep), 'beta')[,2],v=F)),
#        aes(x=name, y=slope_est))+
#   geom_errorbar(aes(ymin=slope_est-slope_err, ymax=slope_est+slope_err), width=.1)+
#   geom_line() + geom_point()+theme_linedraw()+labs(x='Name of beta slope', y='Estimate for beta value')+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
# ggsave(paste0(folder_images_out, "full_partialILR_nos5_slope_errorbar.png"), width = 3, height = 2.5)
# 
# ## only matched samples
# ## TRUE_data: subset of signatures. Removing s5
# both_arx_rlps <- sapply(unique(patient.meta$PATIENT_ID), function(i) all(c('rlps', 'arx') %in% 
#                                                                            patient.meta$group[patient.meta$PATIENT_ID == i]))
# both_arx_rlps <- rownames(exposures) %in% patient.meta$SAMPLE_ID[patient.meta$PATIENT_ID %in% patient.meta$PATIENT_ID[both_arx_rlps]]
# dim(exposures)
# length(both_arx_rlps)
# TMB_data_with_subset_only_matched = tmb_MVN_partial_ILR_FEb(c(1:4,6:7),
#                                                             .keep_additional = both_arx_rlps)
# rownames(TMB_data_with_subset_only_matched$Y)
# TMB_data_with_subset_only_matched$x
# pheatmap(TMB_data_with_subset_only_matched$Y)
# TMB_params_with_subset_only_matched= give_TMB_params(arg_d = TMB_data_with_subset_only_matched$d,
#                                                      arg_num_individuals = TMB_data_with_subset_only_matched$num_individuals)
# TMB_params_with_subset_only_matched$logsd = TMB_params_with_subset_only_matched$logs_sd_RE
# TMB_params_with_subset_only_matched$cov_RE <- NULL
# TMB_params_with_subset_only_matched$logs_sd_RE <- NULL
# obj_only_matched <- MakeADFun(data = TMB_data_with_subset_only_matched,
#                               parameters = TMB_params_with_subset_only_matched, DLL="tmb_MVN_partial_ILR_FEb")
# opt_only_matched <- do.call("optim", obj_only_matched)
# opt_only_matched
# opt_only_matched$hessian ## <-- FD hessian from optim
# rep_only_matched <- sdreport(obj_only_matched)
# rep_only_matched
# wald_TMB_wrapper(rep_only_matched)
# saveRDS(rep_only_matched, paste0(folder_out_RDS, "repeat_full_nos5_onlymatched_FE.RDS"))
# 
# dim(TMB_data_with_subset$Y)
# dim(TMB_data_with_subset_only_matched$Y)
# 
# #-------------------------------------------------------------------------------------------#
# 
# pdf(paste0(folder_images_out, "exposures_ternary.pdf"), width = 6, height = 2.5)
# par(mfrow=c(1,2), mar=c(0,0,1,0))
# plot_ternary(exposures[patient.meta$group == "arx",c(1,3,5)], legend_on = F, main='Archival')
# plot_ternary(exposures[patient.meta$group == "rlps",c(1,3,5)], legend_on = F, main='Relapse')
# dev.off()
# 
# remove_NA <- function(i) i[!apply(apply(i, 1, is.na),2, any),]
# pdf(paste0(folder_images_out, "exposures_ternary_2.pdf"), width = 6, height = 2.5)
# par(mfrow=c(1,2), mar=c(0,0,1,0))
# plot_ternary(remove_NA(normalise_rw(exposures[patient.meta$group == "arx",c(3,5,7)])),
#              legend_on = F, main='Archival')
# plot_ternary(remove_NA(normalise_rw(exposures[patient.meta$group == "rlps",c(3,5,7)])),
#              legend_on = F, main='Relapse')
# dev.off()
# 
# pdf(paste0(folder_images_out, "exposures_ternary_3.pdf"), width = 6, height = 2.5)
# par(mfrow=c(1,2), mar=c(0,0,1,0))
# plot_ternary(remove_NA(normalise_rw(exposures[patient.meta$group == "arx",c(1,3,7)])),
#              legend_on = F, main='Archival')
# plot_ternary(remove_NA(normalise_rw(exposures[patient.meta$group == "rlps",c(1,3,7)])),
#              legend_on = F, main='Relapse')
# dev.off()
# 
# exposures
# 
# match(rownames(exposures),
# patient.meta$SAMPLE_ID)
# 
# rownames(exposures_partial_irl)
# patient.meta$SAMPLE_ID
# 
# all(rownames(exposures) == patient.meta$SAMPLE_ID)
# exposures_df <- (melt(exposures))
# exposures_df$group <- patient.meta$group[match(exposures_df$Var1, patient.meta$SAMPLE_ID)]
# exposures_df$patient <- patient.meta$PATIENT_ID[match(exposures_df$Var1, patient.meta$SAMPLE_ID)]
# 
# exposures_df$group[exposures_df$group == 'arx'] <- 'Archival'
# exposures_df$group[exposures_df$group == 'rlps'] <- 'Relapse'
# ggplot(exposures_df, aes(x=Var2, y=value, col=group, group=interaction(Var2,group)))+
#   geom_boxplot()+ geom_point(position=position_jitterdodge(),
#                              alpha=0.2, size=0.9, width = 0.2)+
#   theme_bw()+labs(x='Signature', y='Exposure', col='Group')
