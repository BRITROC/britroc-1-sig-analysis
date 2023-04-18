## Using partialILR to deal with zero exposures, applied to the BriTROC-1 dataset

rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# set.seed(1234)

library(uuid)
library(ggplot2)
library(reshape2)
library(compositions)
library(TMB)
library(gridExtra)
library(Ternary)

# source("../../../GlobalDA/code/2_inference/helper/helper_DA_stan.R")
# source("../../../GlobalDA/code/3_analysis/helper/helper_analyse_posteriors.R")
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
which_zero = t(apply(exposures, 1, function(i) as.numeric((i==0)) ))

exposures[1,]
which_zero[1,]
as.vector(compositions::ilr(exposures)[1,])
give_partial_ilr_basis(which(which_zero[1,] == 1), d=7)
give_partial_ilr_basis(which_zero_vector = 2, d = 5)

give_partial_ilr_basis(which_zero_vector = 2, d = 4)

.xx <- give_partial_ilr_basis(which_zero_vector = 2, d = 5)
sum(.xx[,1]*.xx[,2])
sum(.xx[,1]*.xx[,3])
sum(.xx[,1]*.xx[,4])
for(i in 1:4) print(sum(.xx[,i]*.xx[,i]))

exposures_partial_irl0 = give_partial_irl(exposures)

dim(exposures)
dim(exposures_partial_irl0)

## note that the partial ILR that I have  computed is different from the ILR setting zero entries as zero ILRs
## different
# compositions::clr(exposures[1,]) %*% give_partial_ilr_basis(which_zero[[1]])
# compositions::clr(exposures[1,]) %*% irl_base_complete

## i.e the zero entries are different

# To show that I am computing the ILR correctly given a basis
# plot(as.vector(compositions::ilr(exposures)[1,]),
#      compositions::clr(exposures[1,]) %*% irl_base_complete)
# abline(coef = c(0,1))

#-------------------------------------------------------------------------------------------#
## random effects partial ILR
TMB::compile("../tmb_RE/tmb_MVN_partial_ILR.cpp", "-std=gnu++17")
dyn.load(dynlib("../tmb_RE/tmb_MVN_partial_ILR"))
## fixed effects partial ILR
TMB::compile("../tmb_RE/tmb_MVN_partial_ILR_FEb.cpp", "-std=gnu++17")
dyn.load(dynlib("../tmb_RE/tmb_MVN_partial_ILR_FEb"))
## fitting a MVN with FE on the mean
TMB::compile("../tmb_RE/tmb_MVN_with_mean_2.cpp", "-std=gnu++17")
dyn.load(dynlib("../tmb_RE/tmb_MVN_with_mean_2"))
## MVN for ILR or ALR without zeros (ME, correlated)
TMB::compile("../tmb_RE/tmb_MVN_ILR.cpp", "-std=gnu++17")
dyn.load(dynlib("../tmb_RE/tmb_MVN_ILR"))

## random effects partial ILR, no correlations
TMB::compile("../tmb_RE/tmb_MVN_partial_ILR_notcor.cpp", "-std=gnu++17")
dyn.load(dynlib("../tmb_RE/tmb_MVN_partial_ILR_notcor"))
## random effects partial ILR, no correlations, overdispersion
TMB::compile("../tmb_RE/tmb_MVN_partial_ILR_notcor_overdisp.cpp", "-std=gnu++17")
dyn.load(dynlib("../tmb_RE/tmb_MVN_partial_ILR_notcor_overdisp"))
## fixed-effects partial ILR, with correlations
TMB::compile("../tmb_RE/tmb_MVN_partial_ILR_FEe.cpp", "-std=gnu++17")
dyn.load(dynlib("../tmb_RE/tmb_MVN_partial_ILR_FEe"))

TMB::compile("../tmb_RE/tmb_MVN_partial_ILR_notcor_outsidesd.cpp", "-std=gnu++17")
dyn.load(dynlib("../tmb_RE/tmb_MVN_partial_ILR_notcor_outsidesd"))
TMB::compile("../tmb_RE/tmb_MVN_partial_ILR_outsidesd.cpp", "-std=gnu++17")
dyn.load(dynlib("../tmb_RE/tmb_MVN_partial_ILR_outsidesd"))

#-------------------------------------------------------------------------------------------#


#-------------------------------------------------------------------------------------------#
## remove the samples in which there are too many zeros. There can be at most 4 zeros (i.e. there has to be at least one log-ratio)
keep_britroc_samples = (rowSums(exposures_partial_irl0 == 0) < (ncol(exposures) - 2) )

exposures_partial_irl0[!keep_britroc_samples,] ## sample with only s1 and s3
exposures[!keep_britroc_samples,] ## sample with only s1

exposures_partial_irl = exposures_partial_irl0[keep_britroc_samples,]
colnames(exposures_partial_irl) <- paste0('ILR', 1:ncol(exposures_partial_irl))
all(patient.meta$SAMPLE_ID[keep_britroc_samples] == rownames(exposures_partial_irl))
z_britroc = give_z_matrix_from_labels(patient.meta$PATIENT_ID[keep_britroc_samples])

## the number of samples
nrow(exposures_partial_irl)

## the number of samples per patient
table(colSums(z_britroc))

TMB_data = list(Y = exposures_partial_irl,
                num_individuals = ncol(z_britroc),
                d = d,
                n = nrow(exposures_partial_irl),
                x = cbind(1, as.numeric(as.factor(patient.meta$group[keep_britroc_samples]))-1),
                z = z_britroc)

TMB_data_all_samples = list(Y = exposures_partial_irl0,
                num_individuals = ncol(give_z_matrix_from_labels(patient.meta$PATIENT_ID)),
                d = d,
                n = nrow(exposures_partial_irl0),
                x = cbind(1, as.numeric(as.factor(patient.meta$group))-1),
                z = give_z_matrix_from_labels(patient.meta$PATIENT_ID))
give_summary_dimensions_TMBobj(TMB_data)

#-------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------#
## TRUE_data: full
TMB_params = give_TMB_params(d, TMB_data$num_individuals)
sapply(TMB_params, dim)

## below: to optimisise with optim. I insetad use a wrapper that uses nlminb 
# obj <- MakeADFun(data = TMB_data, parameters = TMB_params, DLL="tmb_MVN_partial_ILR", random = "u_large")
# opt <- do.call("optim", obj)
# opt
# opt$hessian ## <-- FD hessian from optim
# rep <- sdreport(obj)
# rep

## <partialILR> uses tmb_MVN_partial_ILR.cpp
res_nlminb <- wrapper_run_TMB_use_nlminb(model = "partialILR", object = TMB_data,
                                         use_nlminb = T, iter.max=10000)
res_nlminb
saveRDS(res_nlminb, paste0(folder_out_RDS, "partialILRcor_nlminb_allsigs.RDS"))

## it has not converged, seemingly because of the covariances in the random effects
## note that the values of the matrix L (to be transformed into covariance) are of very different
## magnitudes, and have high errors

## apparent problem with the last category: it is very highly correlated with the first ILR
give_pairs_with_mvn_wrapper(matrix(res_nlminb$par.random, ncol=TMB_data$d), common_lims = T)
give_pairs_with_mvn_wrapper(matrix(res_nlminb$par.random, ncol=TMB_data$d), common_lims = F) ## no problem with differences in magnitudes

par(mfrow=c(1,1), mar=c(1,1,1,1))
plot(matrix(res_nlminb$par.random, ncol=TMB_data$d)[,c(1,6)])

## it does not seem to be a problem particular to any sample, i.e. there is no clear outlier

## there are no obvious correlations in the ILR values used as input
give_pairs_with_mvn_wrapper(TMB_data$Y)
give_pairs_with_mvn_wrapper(TMB_data$Y, zero_to_NA = T, common_lims = F)
give_pairs_with_mvn_wrapper(TMB_data$Y, zero_to_NA = T, common_lims = T)

RE_res_nlminb <- matrix(res_nlminb$par.random, ncol=TMB_data$d)

fitted_res_nlminb <- TMB_data$x %*% matrix(python_like_select_name(res_nlminb$par.fixed, 'beta'), nrow=2)+
  TMB_data$z %*% RE_res_nlminb

## the RE of the first ILR correlate highly with many, e.g. with the last one
# give_pairs_with_mvn_wrapper(matrix(res_nlminb$par.random, ncol=TMB_data$d), common_lims = F)
give_pairs_with_mvn_wrapper(RE_res_nlminb, common_lims = F)

compare_matrices(mat1 = TMB_data$Y,
                 mat2 = add_rownames(add_colnames(fitted_res_nlminb, colnames(TMB_data$Y)),
                                     rownames(TMB_data$Y)), remove_zeros=T, facets = T,
                 groups=TMB_data$x[,2]) ## removing zero entries

## running the partialILR model without correlations, the model converges well
res_nlminb_nocoroutsidesd <- wrapper_run_TMB_use_nlminb(model = "partialILRnocoroutsidesd", object = TMB_data,
                                                        use_nlminb = T, iter.max=10000)
res_nlminb_outsidesd <- wrapper_run_TMB_use_nlminb(model = "partialILRoutsidesd", object = TMB_data,
                                                        use_nlminb = T, iter.max=10000)
res_nlminb_nocor <- wrapper_run_TMB_use_nlminb(model = "partialILRnocor", object = TMB_data,
                                               use_nlminb = T, iter.max=10000)
res_nlminb_nocor
saveRDS(res_nlminb_nocor, paste0(folder_out_RDS, "partialILRnocor_nlminb_allsigs.RDS"))
saveRDS(TMB_data, paste0(folder_out_RDS, "partialILRnocor_nlminb_allsigs_TMBdata.RDS"))

## very clearly some log-ratios have larger intercepts than others. ILR4 especially important
give_pairs_with_mvn_wrapper(matrix(res_nlminb_nocor$par.random, ncol=TMB_data$d), common_lims = T)

grid.arrange(plot_betas(res_nlminb), plot_betas(res_nlminb_nocor))

res_nlminb_nocoroverdisp <- wrapper_run_TMB_use_nlminb(model = "partialILRnocoroverdisp", object = TMB_data,
                                               use_nlminb = T, iter.max=10000)
res_nlminb_nocoroverdisp
res_nlminb_nocoroverdisp$par.fixed

grid.arrange(plot_betas(res_nlminb), plot_betas(res_nlminb_nocor), plot_betas(res_nlminb_nocoroverdisp),
             plot_betas(res_nlminb_outsidesd), ncol=1)
## betas are nearly identical
wald_TMB_wrapper(res_nlminb_nocor)
wald_TMB_wrapper(res_nlminb_nocoroverdisp)
wald_TMB_wrapper(res_nlminb_outsidesd)

res_nlminb_FEcor <- wrapper_run_TMB_use_nlminb(model = "tmb_MVN_partial_ILR_FEe", object = TMB_data,
                                                       use_nlminb = T, iter.max=10000)
res_nlminb_FEcor

## the random effect model was not statistically significant whereas the FE model with no correlation was
## this is contrary to the results in SNV, where adding RE increases power. In CN in the model assessment the FE
## consistently have higher sensitivity than the ME models, when simulating using britroc data
wald_TMB_wrapper(res_nlminb_FEcor)
wald_TMB_wrapper(res_nlminb_nocor)
wald_TMB_wrapper(res_nlminb_nocoroutsidesd)

grid.arrange(compare_TMB_fit_to_data(res_nlminb, TMB_data, remove_zeros=T, title='res_nlminb'),
             compare_TMB_fit_to_data(res_nlminb_nocor, TMB_data, remove_zeros=T, title='res_nlminb_nocor'),
             compare_TMB_fit_to_data(res_nlminb_nocoroverdisp, TMB_data, remove_zeros=T, title='res_nlminb_nocoroverdisp'))
## arguably the best fit is that of the overdispersed, uncorrelated, version

pairs(matrix(res_nlminb_nocor$par.random, ncol=TMB_data$d))
pairs(matrix(res_nlminb_nocoroverdisp$par.random, ncol=TMB_data$d))
compare_TMB_fits(res_nlminb_nocor, res_nlminb_nocoroverdisp, TMB_data)

##----------------------------------------------------------------------------------------------##


##----------------------------------------------------------------------------------------------##
## Keeping the ILR as they are, i.e. without re-running them, fit model using a subset of ILRs  
TMB_data_subset <- TMBobj_partialILR_remove_single_obs(give_subset_sigs_TMBobj(give_subset_samples_TMBobj_2(TMB_data, 1:nrow(TMB_data$Y)),
                                                            c('ILR2', 'ILR3', 'ILR5'), remove_zero_rows = F))
give_summary_dimensions_TMBobj(TMB_data_subset)
give_pairs_with_mvn_wrapper(TMB_data_subset$Y)
give_pairs_with_mvn_wrapper(TMB_data_subset$Y, zero_to_NA = T)
res_nlminb_subset <- wrapper_run_TMB_use_nlminb(model = "partialILR", object = TMB_data_subset,
                                                           use_nlminb = T, iter.max=1000)

## tmb_MVN_ILR: MVN for ILR or ALR without zeros (ME, correlated).
sum(TMB_data_subset$Y == 0) ## there are some zeros
give_summary_dimensions_TMBobj(TMB_data_subset)

## There are zeros, so I should not be using this model
# res_nlminb_subset_b <- wrapper_run_TMB_use_nlminb(model = "tmb_MVN_ILR", object = TMB_data_subset,
#                                                 use_nlminb = T, iter.max=10000)

## extracting the random effect coefficients from the partialILR model
RE_res_nlminb_subset <- matrix(res_nlminb_subset$par.random, ncol=TMB_data_subset$d)
give_pairs_with_mvn_wrapper(RE_res_nlminb_subset, common_lims = T) ## extremely correlated
pheatmap(RE_res_nlminb_subset)
apply(RE_res_nlminb_subset, 2, max) ## there are no random effects for the third category, i.e. no variance between patients

plot(TMB_data_subset$Y[,1:2], col=factor(TMB_data_subset$x[,2]))
plot(TMB_data_subset$Y[,1:3], col=factor(TMB_data_subset$x[,2]))

## random effect coefficients from the tmb_MVN_ILR model (only run if there are no zeros)
# RE_res_nlminb_subset_b <- matrix(res_nlminb_subset_b$par.random, ncol=TMB_data_subset$d)
## RE are all the same
# pheatmap(RE_res_nlminb_subset_b, cluster_cols = F) ## only zeros in the second RE
# give_pairs_with_mvn_wrapper(RE_res_nlminb_subset_b, common_lims = F)
# dev.off()
# plot(RE_res_nlminb_subset_b[,1:2])
# plot(RE_res_nlminb_subset_b[,1:3])
## an option is that the random effects are negligible compared to the fixed effects

## getting fitted values
fitted_res_nlminb_subset <- TMB_data_subset$x %*% matrix(python_like_select_name(res_nlminb_subset$par.fixed, 'beta'), nrow=2)+
  TMB_data_subset$z %*% RE_res_nlminb_subset
# fitted_res_nlminb_subset_b <- TMB_data_subset$x %*% matrix(python_like_select_name(res_nlminb_subset_b$par.fixed, 'beta'), nrow=2)+
#   TMB_data_subset$z %*% RE_res_nlminb_subset_b

## some correlations but none too extreme
sort(RE_res_nlminb_subset[,1], decreasing = T)[1:10]
par(mfrow=c(1,3))
hist(RE_res_nlminb_subset[,1], breaks = 20)
hist(RE_res_nlminb_subset[,2], breaks = 20)
hist(RE_res_nlminb_subset[,3], breaks = 20) ## they all look good and following a normal distrib
# dev.off()
# par(mfrow=c(1,3))
# plot(density(RE_res_nlminb_subset_b[,1]))
# plot(density(RE_res_nlminb_subset_b[,2]))
# plot(density(RE_res_nlminb_subset_b[,3]))

## comparing fitted values to observed values
##' a large fraction of ILR4 are zero (or very close to zero) in the observed category. These are the
##' observed zeros which have not been taken into account in the partialILR analysis
compare_matrices(mat1 = TMB_data_subset$Y,
                 mat2 = add_rownames(add_colnames(fitted_res_nlminb_subset, colnames(TMB_data_subset$Y)),
                                     rownames(TMB_data_subset$Y)), remove_zeros=F)
compare_matrices(mat1 = TMB_data_subset$Y,
                 mat2 = add_rownames(add_colnames(fitted_res_nlminb_subset, colnames(TMB_data_subset$Y)),
                                     rownames(TMB_data_subset$Y)), remove_zeros=T) ## removing zero entries
# compare_matrices(mat1 = TMB_data_subset$Y,
#                  mat2 = add_rownames(add_colnames(fitted_res_nlminb_subset_b, colnames(TMB_data_subset$Y)),
#                                      rownames(TMB_data_subset$Y)))
## there are a few samples with very extreme random intercepts for the first category


## no clear correlations in the observed exposures
give_pairs_with_mvn_wrapper(TMB_data_subset$Y, zero_to_NA=T)
give_pairs_with_mvn_wrapper(RE_res_nlminb_subset, zero_to_NA=T) ## clear correlations in RE intercepts
# give_pairs_with_mvn_wrapper(RE_res_nlminb_subset_b, zero_to_NA=T) ## clear correlations, and all zero in the second category

## Comparison of RE intercepts in the two models
# compare_matrices(RE_res_nlminb_subset, RE_res_nlminb_subset_b) ## difference in the second cat, which is always zero in RE_res_nlminb_subset_b

## see how many combinations of two ILRs we have in each group
outer(1:ncol(TMB_data_subset$Y), 1:ncol(TMB_data_subset$Y), Vectorize(function(i,j){
  .x <- TMB_data_subset$Y[TMB_data_subset$x[,2] == 0,]
  sum ((.x[,i]  != 0) &  (.x[,j] != 0) )
})) ## for the first group
outer(1:ncol(TMB_data_subset$Y), 1:ncol(TMB_data_subset$Y), Vectorize(function(i,j){
  .x <- TMB_data_subset$Y[TMB_data_subset$x[,2] == 1,]
  sum ((.x[,i]  != 0) &  (.x[,j] != 0) )
})) ## for the second group

pheatmap(TMB_data_subset$Y)


.cov_est <- L_to_cov(python_like_select_name(res_nlminb$par.fixed, 'cov_RE'),
                     d =  ncol(TMB_data_subset$Y)) ## estimates
.cov_est
.cov_est <- L_to_cov(python_like_select_rownames(summary(res_nlminb), 'cov_RE')[,2],
                     d =  ncol(TMB_data_subset$Y)) ## se
.cov_est
##----------------------------------------------------------------------------------------------##


##----------------------------------------------------------------------------------------------##
### run with all sigs
# res_nlminb0 <- wrapper_run_TMB_use_nlminb(model = "partialILR",
#                                          object = TMBobj_partialILR_remove_single_obs(TMB_data), use_nlminb = T,
#                                          iter.max=1000)
# res_nlminb <- wrapper_run_TMB_use_nlminb(model = "partialILR",
#                                          object = TMBobj_partialILR_remove_single_obs(TMB_data), use_nlminb = T,
#                                          iter.max=1000, initial_params = list(beta = matrix(python_like_select_name(res_nlminb0$par.fixed, 'beta'), nrow=2),
#                                                                               u_large = matrix(res_nlminb0$par.random, ncol=ncol(TMB_data$Y)),
#                                                                               logs_sd_RE=python_like_select_name(res_nlminb0$par.fixed, 'logs_sd_RE'),
#                                                                               cov_RE = python_like_select_name(res_nlminb0$par.fixed, 'cov_RE')))
# res_nlminb
# saveRDS(res_nlminb, paste0(folder_out_RDS, "res_nlminb.RDS"))
# res_nlminb_from_saved <- readRDS(paste0(folder_out_RDS, "res_nlminb.RDS"))
# give_pairs_with_mvn_wrapper(matrix(res_nlminb_from_saved$par.random, ncol=TMB_data$d), common_lims = F)

# repeat_full_with_different_initial = function(){
#   TMB_params = give_TMB_params(d, TMB_data$num_individuals)
#   obj <- MakeADFun(data = TMB_data, parameters = TMB_params, DLL="tmb_MVN_partial_ILR", random = "u_large")
#   opt <- do.call("optim", obj)
#   opt
#   opt$hessian ## <-- FD hessian from optim
#   rep <- sdreport(obj)
#   rep
# }
# 
# repeat_full_with_different_initial_obs = replicate(n = 10, repeat_full_with_different_initial())
# 
# saveRDS(repeat_full_with_different_initial_obs, paste0(folder_out_RDS, "repeat_full_with_different_initial_obs.RDS"))
## we have one pdHess true!!
# unlist(apply(repeat_full_with_different_initial_obs, 2, `[`, 'pdHess'))

# repeat_full_with_different_initial_obs[,9]$pdHess

## Beta values are the same regardless of whether the runs have converged or not
# pdf(paste0(folder_images_out, "full_partialILR_several_runs.pdf"))
# pairs(repeat_full_with_different_initial_obs['par.fixed',])
# dev.off()
# pdf(paste0(folder_images_out, "full_partialILR_several_runs_only_beta.pdf"))
# pairs(sapply(repeat_full_with_different_initial_obs['par.fixed',], python_like_select_name, grep_substring="beta"),
#       main='Only beta (intercept+slope)')
# dev.off()


## Analyse the results of the converged run
results_full0 = res_nlminb#repeat_full_with_different_initial_obs[,9]
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

##-------------------------------------------------------------------------------------------##
## TRUE_data: subset of signatures. Removing s5
both_arx_rlps <- sapply(unique(patient.meta$PATIENT_ID), function(i) all(c('rlps', 'arx') %in% 
                                                                           patient.meta$group[patient.meta$PATIENT_ID == i]))
both_arx_rlps <- rownames(exposures) %in% patient.meta$SAMPLE_ID[patient.meta$PATIENT_ID %in% patient.meta$PATIENT_ID[both_arx_rlps]]

## subset of samples for non-duplicated patients (i.e. for fixed effects, so that we don't include replicates). Selecting the first
single_group_sample_patient <- sapply(unique(patient.meta$PATIENT_ID), function(i) sapply(c('rlps', 'arx'), function(j) which( 
                                                                           (patient.meta$PATIENT_ID == i) & (patient.meta$group == j))[1]))
single_group_sample_patient_bool <- 1:nrow(patient.meta) %in% single_group_sample_patient

all(as.character(patient.meta$SAMPLE_ID) == rownames(exposures))
# both_arx_rlps <- rep(T, length(both_arx_rlps))

table(both_arx_rlps)
dim(patient.meta)
# both_arx_rlps[which(both_arx_rlps)[1]] <- F
##-------------------------------------------------------------------------------------------##

##-------------------------------------------------------------------------------------------##
## Only samples with matched archival and relapse
# keep_only_both <- patient.meta$PATIENT_ID %in% 
TMB_data_only_matched = prepare_TMB_data_with_subset(c(1:7), exposures=exposures,
                                                                 .keep_additional = both_arx_rlps)
saveRDS(TMB_data_only_matched, paste0("../../../out/inference/partialILR/TMB_data_only_matched_res_nlminb_nocoroutsidesd_only_matched_allsigs_RData.RDS"))
TMB_data_with_subset_only_matched = prepare_TMB_data_with_subset(c(1:4,6:7), exposures=exposures,
                                                                 .keep_additional = both_arx_rlps)
TMB_data_with_subset_nonrepeated = prepare_TMB_data_with_subset(c(1:4,6:7), exposures=exposures, ## without repeated patients in the same group
                                                                 .keep_additional = single_group_sample_patient_bool)
nrow(TMB_data_with_subset_only_matched$Y)
nrow(TMB_data_with_subset_nonrepeated$Y)
rownames(TMB_data_with_subset_only_matched$Y)
TMB_data_with_subset_only_matched$x
pheatmap(TMB_data_with_subset_only_matched$Y)
TMB_params_with_subset_only_matched= give_TMB_params(arg_d = TMB_data_with_subset_only_matched$d,
                                        arg_num_individuals = TMB_data_with_subset_only_matched$num_individuals)
obj_only_matched <- MakeADFun(data = TMB_data_with_subset_only_matched,
                              parameters = TMB_params_with_subset_only_matched, DLL="tmb_MVN_partial_ILR", random = "u_large")
opt_only_matched  <- do.call("optim", obj_only_matched)
opt_only_matched 
opt_only_matched$hessian ## <-- FD hessian from optim
rep_only_matched  <- sdreport(obj_only_matched )
rep_only_matched

saveRDS(rep, paste0(folder_out_RDS, "repeat_full_nos5_onlymatched.RDS"))
wald_TMB_wrapper(rep_only_matched)

res_nlminb_nocoroutsidesd_only_matched_allsigs <- wrapper_run_TMB_use_nlminb(model = "partialILRnocoroutsidesd",
                                                                     object = TMB_data_only_matched,
                                                                     use_nlminb = T, iter.max=10000)
res_nlminb_outsidesd_only_matched_allsigs <- wrapper_run_TMB_use_nlminb(model = "partialILRoutsidesd",
                                                                object = TMB_data_only_matched,
                                                                use_nlminb = T, iter.max=10000)
saveRDS(res_nlminb_nocoroutsidesd_only_matched_allsigs, paste0(folder_out_RDS, "res_nlminb_nocoroutsidesd_only_matched_allsigs.RDS"))
saveRDS(res_nlminb_outsidesd_only_matched_allsigs, paste0(folder_out_RDS, "res_nlminb_outsidesd_only_matched_allsigs.RDS"))
wald_TMB_wrapper(res_nlminb_nocoroutsidesd_only_matched_allsigs)
# [,1]
# [1,] 0.007668402
## only matched , without s5
res_nlminb_nocoroutsidesd_only_matched <- wrapper_run_TMB_use_nlminb(model = "partialILRnocoroutsidesd",
                                                                     object = TMB_data_with_subset_only_matched,
                                                        use_nlminb = T, iter.max=10000)
res_nlminb_outsidesd_only_matched <- wrapper_run_TMB_use_nlminb(model = "partialILRoutsidesd",
                                                                     object = TMB_data_with_subset_only_matched,
                                                                     use_nlminb = T, iter.max=10000)
saveRDS(res_nlminb_nocoroutsidesd_only_matched, paste0(folder_out_RDS, "res_nlminb_nocoroutsidesd_only_matched.RDS"))
saveRDS(res_nlminb_outsidesd_only_matched, paste0(folder_out_RDS, "res_nlminb_outsidesd_only_matched.RDS"))

res_nlminb_nocoroutsidesd_only_matched
wald_TMB_wrapper(res_nlminb_nocoroutsidesd_only_matched)
wald_TMB_wrapper(res_nlminb_outsidesd_only_matched) # 0.05207677
res_nlminb_outsidesd_only_matched_covariance <- add_diag(add_rownames(add_colnames(L_to_cov(python_like_select_name(res_nlminb_outsidesd_only_matched$par.fixed, 'cov_RE'),
                                            d =  ncol(TMB_data_with_subset_only_matched$Y)), paste0('ILR', c(1:3, 5:6))), paste0('ILR', c(1:3, 5:6))),
         exp(python_like_select_name(res_nlminb_outsidesd_only_matched$par.fixed, 'logs_sd_RE')))
pdf(paste0(folder_images_out, "res_nlminb_outsidesd_only_matched_covariance.pdf"), width = 3.5, height = 3)
print(pheatmap(res_nlminb_outsidesd_only_matched_covariance))
dev.off()
png(paste0(folder_images_out, "res_nlminb_outsidesd_only_matched_covariance.png"), width = 3.5, height = 3, units = "in", res=300)
print(pheatmap(res_nlminb_outsidesd_only_matched_covariance))
dev.off()

give_pairs_with_mvn_wrapper(matrix(rep_only_matched$par.random, ncol=TMB_data_with_subset_only_matched$d), common_lims = F)
## very high correlation of RE between ILR3 and ILR4, and several with ILR1

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
##-------------------------------------------------------------------------------------------##


##----------------------------------------------------------------------------------------------##

##----------------------------------------------------------------------------------------------##

##-------------------------------------------------------------------------------------------##
## Using only samples for which we have both archival and relapse
## excluding s5
TMB_data_with_subset = prepare_TMB_data_with_subset(c(1:4,6:7), exposures=exposures)
dim(TMB_data_with_subset$Y)

give_summary_dimensions_TMBobj(TMB_data_with_subset)

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

give_pairs_with_mvn_wrapper(matrix(rep$par.random, ncol=TMB_data_with_subset$d), common_lims = F)

## same as above, but running nlminb
res_nlminb_nos5 <- wrapper_run_TMB_use_nlminb(model = "partialILR", object = TMB_data_with_subset, use_nlminb = T,
                                         iter.max=1000)
res_nlminb_nos5nocor <- wrapper_run_TMB_use_nlminb(model = "partialILRnocor", object = TMB_data_with_subset, use_nlminb = T,
                                              iter.max=1000)
res_nlminb_nos5
res_nlminb_nos5nocoroutsidesd <- wrapper_run_TMB_use_nlminb(model = "partialILRnocoroutsidesd", object = TMB_data_with_subset, use_nlminb = T,
                                                   iter.max=1000)
res_nlminb_nos5outsidesd <- wrapper_run_TMB_use_nlminb(model = "partialILRoutsidesd", object = TMB_data_with_subset, use_nlminb = T,
                                                            iter.max=1000)

res_nlminb_nos5_nocoroverdisp <- wrapper_run_TMB_use_nlminb(model = "partialILRnocoroverdisp", object = TMB_data_with_subset,
                                                       use_nlminb = T, iter.max=10000)

give_pairs_with_mvn_wrapper(TMB_data_with_subset$Y, zero_to_NA = T)

L_to_cov_wrapper_TMB(res_nlminb_nos5)
summary(res_nlminb_nos5)

random_intercepts <- matrix(res_nlminb_nos5$par.random, ncol=length(python_like_select_name(res_nlminb_nos5$par.fixed, 'beta'))/2)
give_pairs_with_mvn_wrapper(random_intercepts) ## the estimates are very highly correlated

give_pairs_with_mvn_wrapper(RE_res_nlminb) ## the estimates are very highly correlated

wald_TMB_wrapper(rep, fail_non_converged = F)
wald_TMB_wrapper(rep_only_matched, fail_non_converged = F)

res_nlminb_FEcor_nos5 <- wrapper_run_TMB_use_nlminb(model = "tmb_MVN_partial_ILR_FEe", object = TMB_data_with_subset,
                                                    use_nlminb = T, iter.max=10000)
res_nlminb_FEcor_nos5_nonrepeated <- wrapper_run_TMB_use_nlminb(model = "tmb_MVN_partial_ILR_FEe", object = TMB_data_with_subset_nonrepeated,
                                                    use_nlminb = T, iter.max=10000)


res_nlminb_FEcor_nos5
wald_TMB_wrapper(res_nlminb_FEcor_nos5) ## FE cor: 0.007; still statistically significant even after removing s5
wald_TMB_wrapper(res_nlminb_nos5nocoroutsidesd) ## ME no cor (better): 0.006; still statistically significant even after removing s5
wald_TMB_wrapper(res_nlminb_nos5outsidesd) ## ME cor:  0.00299376
wald_TMB_wrapper(res_nlminb_nos5) ## ME cor: does no converge
wald_TMB_wrapper(res_nlminb_nos5nocor) ## ME no cor: 0.09
wald_TMB_wrapper(res_nlminb_nos5_nocoroverdisp) # 0.07879901
wald_TMB_wrapper(res_nlminb_FEcor_nos5_nonrepeated) ## 0.0683059

saveRDS(res_nlminb_FEcor_nos5_nonrepeated, paste0(folder_out_RDS, "res_nlminb_FEcor_nos5_nonrepeated.RDS"))

dev.off()
plot(python_like_select_name(res_nlminb_FEcor_nos5$par.fixed, 'beta'),
     python_like_select_name(res_nlminb_nos5nocor$par.fixed, 'beta'), col=c(1,2)); abline(coef = c(0,1), col='blue', lty='dashed')
plot(python_like_select_name(res_nlminb_FEcor_nos5$par.fixed, 'logsd'),
     (python_like_select_name(res_nlminb_nos5nocor$par.fixed, 'logs_sd_RE'))); abline(coef = c(0,1), col='blue', lty='dashed')

pdf((paste0(folder_images_out, "partialILR_coefficients_comparison_models.pdf")), width = 2.5, height = 7, onefile = F)
grid.arrange(plot_betas(res_nlminb_nos5, title='ME cor (prev)'),
             plot_betas(res_nlminb_nos5nocor, title='ME no cor (prev)'),
             plot_betas(res_nlminb_FEcor_nos5, title = "FE cor"),
             plot_betas(res_nlminb_nos5outsidesd, title = "ME cor (current)"),
             plot_betas(res_nlminb_nos5nocoroutsidesd, title = "ME (current)"),
             ncol=1)
dev.off()

## why do we have a lower power when we are including random intercepts? in Chapter 2 it is the opposite
image(TMB_data_with_subset$z)


random_intercepts_nos5_nocor <- matrix(res_nlminb_nos5outsidesd$par.random, ncol=ncol(TMB_data_with_subset$Y))
colnames(random_intercepts_nos5_nocor) <- paste0('ILR', 1:ncol(random_intercepts_nos5_nocor))
give_pairs_with_mvn_wrapper(random_intercepts_nos5_nocor)
a <- pheatmap(give_rownames(random_intercepts_nos5_nocor, (colnames(TMB_data_with_subset$z))))
dev.off(); a
## very high stratification across patients
dim(TMB_data_with_subset$z)

all(rownames(exposures) == patient.meta$SAMPLE_ID)

## signatures sorted by the dendrogram of random effects
patient.meta$SAMPLE_ID <- as.character(patient.meta$SAMPLE_ID)
patient.meta$PATIENT_ID <- as.character(patient.meta$PATIENT_ID)
# RE_sample_ids <- patient.meta$SAMPLE_ID[unlist(sapply(labels(a$tree_row)[a$tree_row$order], function(i) which(patient.meta$PATIENT_ID == i)))]

first_sample_patients <- apply(TMB_data_with_subset$z, 2, function(i) rownames(TMB_data_with_subset$z)[which(i == 1)[1]]) ## get the name of the first sample in each patient
first_sample_patients <- rbind(first_sample_patients,
                               apply(TMB_data_with_subset$z, 2, function(i) paste0(names(table(patient.meta$group[patient.meta$SAMPLE_ID %in% rownames(TMB_data_with_subset$z)[which(i == 1)]])), collapse = '+')))
first_sample_patients
createBarplot(exposures[first_sample_patients[1,],])
rownames(random_intercepts_nos5_nocor) <- first_sample_patients[1,]

## doesn't say much; all the signatures look similar
ggplot(melt(cbind.data.frame(as(exposures[first_sample_patients[1,],], 'matrix'), sample=first_sample_patients[1,],
                           sampletype=first_sample_patients[2,]), id.vars = c("sample", "sampletype")),
       aes(x=sample, y=value, fill=variable))+geom_bar(stat = "identity")+facet_wrap(.~sampletype, scales = "free_x")+theme_bw()

ggplot(melt(cbind.data.frame(as(random_intercepts_nos5_nocor[first_sample_patients[1,],], 'matrix'), sample=first_sample_patients[1,],
                             sampletype=first_sample_patients[2,]), id.vars = c("sample", "sampletype")),
       aes(x=sample, y=variable, fill=value))+geom_tile(stat = "identity")+facet_wrap(.~sampletype, scales = "free_x")+theme_bw()

ComplexHeatmap::Heatmap(matrix = random_intercepts_nos5_nocor[first_sample_patients[1,],],
                        left_annotation = ComplexHeatmap::HeatmapAnnotation(df= data.frame(first_sample_patients[2,]),
                                                                            name = 'Sample availability', which='row',
                                                                            col = list("arx+rlps"='blue', "arx"='red', "rlps"='orange')))
##' clustering from pheatmap, which is nicer. there is no important trend of how many samples
##' are included in each strata in the RE intercepts
# pdf((paste0(folder_images_out, "partialILRME_no_cor_randomintercepts_heatmap.pdf")), width = 5.5, height = 5)
# print(ComplexHeatmap::Heatmap(matrix = random_intercepts_nos5_nocor[first_sample_patients[1,],][a$tree_row$order,], cluster_rows = F,
#                         left_annotation = ComplexHeatmap::HeatmapAnnotation(df= data.frame(`Sample availability`=first_sample_patients[2,][a$tree_row$order]),
#                                                                             name = 'Sample availability', which='row'), show_row_names = F))
# dev.off()

# theme(axis.title.x=element_blank(),
#       axis.text.x=element_blank(),
#       axis.ticks.x=element_blank())

## there are some patients for whom we only have relapse, but they also have an intercept
## we have 158 random intercepts

pdf("../../../results/partialILRmodelling_ME/pheatmap_random_intercepts_res_nlminb_nos5outsidesd.pdf", width = 4, height = 6)
print(pheatmap(add_rownames(add_colnames(matrix(res_nlminb_nos5outsidesd$par.random, ncol=ncol(TMB_data_with_subset$Y)),
                      paste0('ILR', c(1:3, 5:6))), colnames(TMB_data_with_subset$z))))
dev.off()

tp53 <- read.csv("../../../data/restricted/tp53_opt_panel.csv")
meta_archival <- read.csv("../../../data/restricted/Archival_samples.csv")
clinbrit <- read.csv("../../../data/restricted/clinical/britroc_cohort_patient_data.tsv", sep = "\t")
# clinbrit2 <- read.csv("../../../data/restricted/clinical/britroc_patient_treatment_data.tsv", sep = "\t")

clinbrit <- clinbrit[match(gsub("BRITROC-", "", colnames(TMB_data_with_subset$z)), clinbrit$britroc_number),]
meta_archival <- meta_archival[match(gsub("BRITROC-", "", colnames(TMB_data_with_subset$z)), meta_archival$BritRoc.No),]
patient.meta_summary_at_archival <- patient.meta %>% group_by(PATIENT_ID) %>% filter(group == "arx") %>%
  summarize(meanTP53cn = mean(TP53cn), meanpurity=mean(purity), meanploidy=mean(ploidy))
patient.meta_summary_at_archival <- patient.meta_summary_at_archival[match(colnames(TMB_data_with_subset$z),
                                                                           as.character(patient.meta_summary_at_archival$PATIENT_ID)),]

pdf("../../../results/partialILRmodelling_ME/pheatmap_random_intercepts_res_nlminb_nos5outsidesd_anno.pdf", width = 7, height = 6)
print(ComplexHeatmap::Heatmap(add_rownames(add_colnames(matrix(res_nlminb_nos5outsidesd$par.random, ncol=ncol(TMB_data_with_subset$Y)),
                          paste0('ILR', c(1:3, 5:6))), colnames(TMB_data_with_subset$z)),
                        right_annotation = ComplexHeatmap::HeatmapAnnotation(df = data.frame(age=clinbrit$age,
                                   `Sample availab`=first_sample_patients[2,][colnames(TMB_data_with_subset$z)],
                                   num_samples=colSums(TMB_data_with_subset$z),
                                   patient.meta_summary_at_archival[,-1]),
                                   ploidyHigher2p3=patient.meta_summary_at_archival$meanploidy > 2.3,
                                   which = "row"), show_row_names = F, col = ))
dev.off()
##-------------------------------------------------------------------------------------------##
## Amalgamations for convergence
TMB_obj_exposures <- prepare_TMB_data_with_subset(subset_sigs = colnames(exposures), ilr_trans = F, exposures = exposures)
amalgamated_vec1 <- c( 's5')
amalgamated_vec2 <- c( 's7')
amalgamated_vec3 <- c( 's1', 's6')

amalgamation_1 <- give_amalgamated_exposures_TMBobj(sig_obj = TMB_obj_exposures,
        list_groupings = (c(list(amalgamated_vec1), list(amalgamated_vec2), list(amalgamated_vec3),
                            as.list(colnames(exposures)[!(colnames(exposures) %in% c(amalgamated_vec1, amalgamated_vec2, amalgamated_vec3))]))))
head(amalgamation_1$Y)
TMB_obj_amalgamation <- prepare_TMB_data_with_subset(subset_sigs = colnames(amalgamation_1$Y),
                                                     exposures = amalgamation_1$Y, ilr_trans = T)
colnames(amalgamation_1$Y)

give_pairs_with_mvn_wrapper(TMB_obj_amalgamation$Y, zero_to_NA = T, )
TMB_obj_amalgamation
res_nlminb_amalgamation1 <- wrapper_run_TMB_use_nlminb(model = "partialILR", object = TMB_obj_amalgamation, use_nlminb = T,
                                              iter.max=1000)
res_nlminb_amalgamation1
random_intercepts_amalgamation1 <- matrix(res_nlminb_amalgamation1$par.random, ncol=length(python_like_select_name(res_nlminb_amalgamation1$par.fixed, 'beta'))/2)
give_pairs_with_mvn_wrapper(random_intercepts_amalgamation1) ## the estimates are very highly correlated
random_intercepts_amalgamation1
res_nlminb_amalgamation1

plot_betas(res_nlminb_amalgamation1)
wald_TMB_wrapper(res_nlminb_amalgamation1)

##-------------------------------------------------------------------------------------------##
## Amalgamations for convergence, without s5
TMB_obj_exposures_nos5 <- prepare_TMB_data_with_subset(subset_sigs = colnames(exposures)[-5], ilr_trans = F,
                                                       exposures = exposures)
amalgamated_vec1 <- c( 's3', 's7')
amalgamated_vec2 <- c( 's2', 's4')
amalgamated_vec3 <- c( 's1', 's6')
amalgamation_1_no_s5 <- give_amalgamated_exposures_TMBobj(sig_obj = TMB_obj_exposures_nos5,
                                                    list_groupings = (c(list(amalgamated_vec1), list(amalgamated_vec2), list(amalgamated_vec3),
                                                                        as.list(colnames(TMB_obj_exposures_nos5$Y)[!(colnames(TMB_obj_exposures_nos5$Y) %in% c(amalgamated_vec1, amalgamated_vec2, amalgamated_vec3))]))))
head(amalgamation_1_no_s5$Y)
TMB_obj_amalgamation_no_s5 <- prepare_TMB_data_with_subset(subset_sigs = colnames(amalgamation_1_no_s5$Y),
                                                     exposures = amalgamation_1_no_s5$Y, ilr_trans = T)
colnames(amalgamation_1_no_s5$Y)

give_pairs_with_mvn_wrapper(TMB_obj_amalgamation_no_s5$Y, zero_to_NA = T)
TMB_obj_amalgamation
res_nlminb_amalgamation1_nos5 <- wrapper_run_TMB_use_nlminb(model = "partialILR", object = TMB_obj_amalgamation_no_s5,
                                                            use_nlminb = T, iter.max=1000)
res_nlminb_amalgamation1_nos5
random_intercepts_amalgamation1_nos5 <- matrix(res_nlminb_amalgamation1_nos5$par.random, ncol=length(python_like_select_name(res_nlminb_amalgamation1_nos5$par.fixed, 'beta'))/2)
give_pairs_with_mvn_wrapper(random_intercepts_amalgamation1_nos5) ## the estimates are very highly correlated
random_intercepts_amalgamation1_nos5

plot_betas(res_nlminb_amalgamation1_nos5)
wald_TMB_wrapper(res_nlminb_amalgamation1_nos5)

grid.arrange(plot_betas(res_nlminb), plot_betas(res_nlminbnozeros), plot_betas(res_nlminb2), plot_betas(res_nlminb3))


## quite good
# amalgamated_vec1 <- c( 's3', 's7')
# amalgamated_vec2 <- c( 's1')
amalgamated_vec1 <- c( 's1')
amalgamated_vec2 <- c( 's2')
amalgamated_vec3 <- c( 's3')
no_amalgamation_1_no_s5 <- give_amalgamated_exposures_TMBobj(sig_obj = TMB_obj_exposures_nos5,
                                                          list_groupings = (c(list(amalgamated_vec1), list(amalgamated_vec2), list(amalgamated_vec3),
                                                                              as.list(colnames(TMB_obj_exposures_nos5$Y)[!(colnames(TMB_obj_exposures_nos5$Y) %in% c(amalgamated_vec1, amalgamated_vec2, amalgamated_vec3))]))))
TMB_obj_no_amalgamation_no_s5 <- prepare_TMB_data_with_subset(subset_sigs = colnames(no_amalgamation_1_no_s5$Y),
                                                           exposures = no_amalgamation_1_no_s5$Y, ilr_trans = T)
res_nlminb_no_amalgamation1_nos5 <- wrapper_run_TMB_use_nlminb(model = "partialILR", object = TMB_obj_no_amalgamation_no_s5,
                                                            use_nlminb = T, iter.max=1000)

res_nlminb_no_amalgamation1_nos5
random_intercepts_no_amalgamation1_nos5 <- matrix(res_nlminb_no_amalgamation1_nos5$par.random,
                                                  ncol=length(python_like_select_name(res_nlminb_no_amalgamation1_nos5$par.fixed, 'beta'))/2)
give_pairs_with_mvn_wrapper(random_intercepts_no_amalgamation1_nos5) ## the estimates are very highly correlated

## plot with imputation
pairs(impute(TMB_obj_exposures$Y, 1e-2))
pairs(as(compositions::alr(impute(TMB_obj_exposures$Y, 1e-2)), 'matrix'), col=scales::alpha('black', 0.2), pch=19)
pairs(as(compositions::alr(impute(TMB_obj_exposures$Y, 1e-2), ivar = 1), 'matrix'), col=scales::alpha('black', 0.2), pch=19)
pairs(as(compositions::alr(impute(TMB_obj_exposures$Y, 1e-2), ivar = 2), 'matrix'), col=scales::alpha('black', 0.2), pch=19)

give_pairs_with_mvn_wrapper(as(compositions::alr(impute(TMB_obj_exposures$Y, 1e-2), ivar = 1), 'matrix'))
give_pairs_with_mvn_wrapper(as(compositions::alr(impute(TMB_obj_exposures$Y, 1e-2), ivar = 2), 'matrix'))
give_pairs_with_mvn_wrapper(as(compositions::alr(impute(TMB_obj_exposures$Y, 1e-2), ivar = 4), 'matrix'))
give_pairs_with_mvn_wrapper(as(compositions::alr(impute(TMB_obj_exposures$Y, 1e-2), ivar = 7), 'matrix'))
#-------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------#
## imputation
all(match(patient.meta[keep_britroc_samples,'SAMPLE_ID'], rownames(exposures_partial_irl))  == 1:nrow(exposures_partial_irl))
patient.meta_keep <- patient.meta[keep_britroc_samples,]
dim(exposures_partial_irl)
dim(patient.meta_keep)

exposures_inputation1em2_ILR <- as(compositions::ilr(impute(exposures, inputation_value = 1e-2)), 'matrix')
exposures_inputation1em2_ALR <- as(compositions::alr(impute(exposures, inputation_value = 1e-2)), 'matrix')
exposures_partial_irl_NA <- exposures_partial_irl
exposures_partial_irl_NA[exposures_partial_irl_NA == 0] <- NA

patients_both_samples <- unique(patient.meta_keep$PATIENT_ID)[which(sapply(unique(patient.meta_keep$PATIENT_ID), function(i) all(c('arx', 'rlps') %in% patient.meta_keep$group[patient.meta_keep$PATIENT_ID == i])))]

differences_ILR_fun <- function(exposures_mat){
  lapply(patients_both_samples, function(patient_it){
  ct = 1
  difference_ILR= list()
  for(arx_it in which(patient.meta_keep[patient.meta_keep$PATIENT_ID == patient_it,'group'] == "arx")){
    for(rlps_it in which(patient.meta_keep[patient.meta_keep$PATIENT_ID == patient_it,'group'] == "rlps")){
      sam_arx <- patient.meta_keep[patient.meta_keep$PATIENT_ID == patient_it,][arx_it,'SAMPLE_ID']
      sam_rlps <- patient.meta_keep[patient.meta_keep$PATIENT_ID == patient_it,][rlps_it,'SAMPLE_ID']
      difference_ILR[[paste0(sam_arx, '/', sam_rlps)]]= exposures_mat[as.character(sam_rlps),]- exposures_mat[as.character(sam_arx),]
      ct = ct+1
    }
  }
  return(difference_ILR)
})
}
differences_ILR <- differences_ILR_fun(exposures_partial_irl_NA)
differences_ILR_imput <- differences_ILR_fun(exposures_inputation1em2_ILR)

hist(exposures[patient.meta_keep$group == 'arx','s3'] )
dev.off(); par(mfrow=c(1,3))
plot(density(exposures[patient.meta_keep$group == 'arx','s3'] ))
plot(density(exposures[patient.meta_keep$group == 'arx','s4'] ))
plot(density(exposures[patient.meta_keep$group == 'arx','s7'] ))
abline(v = 0.1)
arx_without_s3 <- names(which(exposures[patient.meta_keep$group == 'arx','s3'] == 0))
arx_wihout_s7 <- names(which(exposures[patient.meta_keep$group == 'arx','s7'] <= 0.1))
patients_without_arx_hrd <- patient.meta_keep$PATIENT_ID[match(arx_without_s3, patient.meta_keep$SAMPLE_ID)]
patients_without_arx_hrd

differences_ILR_rbind <- do.call('rbind', do.call('c', differences_ILR))
differences_ILR_melt <- melt(differences_ILR_rbind)
differences_ILR_melt$Var2 = paste0('ILR', differences_ILR_melt$Var2)
differences_ILR_melt$arx_without_s3 = sapply(differences_ILR_melt$Var1, function(j) any(sapply(arx_without_s3, function(i) grepl(i, j))) )
differences_ILR_melt$arx_wihout_s7 = sapply(differences_ILR_melt$Var1, function(j) any(sapply(arx_wihout_s7, function(i) grepl(i, j))) )

ggplot(melt(differences_ILR_melt), aes(x=value, group=interaction(Var2,arx_without_s3),
                                       lty=arx_without_s3))+geom_density()+facet_wrap(.~Var2)+theme_bw()+
  geom_vline(xintercept = 0, lty='dashed')+theme(legend.position = "bottom")
ggsave("../../../results/partialILRmodelling_ME/differences_in_partialILR_HRDs3stratification.pdf", width = 6, height = 4)

ggplot(melt(differences_ILR_melt), aes(x=value, group=interaction(Var2,arx_wihout_s7),
                                       lty=arx_wihout_s7))+geom_density()+facet_wrap(.~Var2)+theme_bw()+
  geom_vline(xintercept = 0, lty='dashed')
ggsave("../../../results/partialILRmodelling_ME/differences_in_partialILR_HRDs7stratification.pdf", width = 6, height = 4)

table(differences_ILR_melt$arx_without_s3)
table(patient.meta_keep$PATIENT_ID %in% patients_without_arx_hrd)

## in samples where the archival is diploid (i.e. zero s7), the relapse can either gain ILR3 (undergo WGD?) or stay diploid
## samples where the archival is already WGD (high s7 in archival) there can be a loss of s4 due to (possibly) increase in some other signature
## ILR1 (s1/s2) is more extreme in changes in the group where there is WGD in the archival
## in samples where there is no WGD (no s7) in the archival, there is an increase of s7 in the relapse (i.e. undergo WGD?)

pairs(differences_ILR_rbind)

plot(differences_ILR_rbind[,c(2,6)])


#--- With imputated data
differences_ILR_rbind_imput <- do.call('rbind', do.call('c', differences_ILR_imput))
differences_ILR_melt_imput <- melt(differences_ILR_rbind_imput)
differences_ILR_melt_imput$Var2 = paste0('ILR', differences_ILR_melt_imput$Var2)
differences_ILR_melt_imput$arx_without_s3 = sapply(differences_ILR_melt_imput$Var1, function(j) any(sapply(arx_without_s3, function(i) grepl(i, j))) )
differences_ILR_melt_imput$arx_wihout_s7 = sapply(differences_ILR_melt_imput$Var1, function(j) any(sapply(arx_wihout_s7, function(i) grepl(i, j))) )

ggplot(melt(differences_ILR_melt_imput), aes(x=value, col=Var2, group=interaction(Var2,arx_without_s3),
                                       lty=arx_without_s3))+geom_density()+facet_wrap(.~Var2)+theme_bw()+
  geom_vline(xintercept = 0, lty='dashed')
ggsave("../../../results/partialILRmodelling_ME/differences_in_partialILR_HRDs3stratification_imput.pdf", width = 6, height = 4)

ggplot(melt(differences_ILR_melt_imput), aes(x=value, group=interaction(Var2,arx_wihout_s7),
                                       lty=arx_wihout_s7))+geom_density()+facet_wrap(.~Var2)+theme_bw()+
  geom_vline(xintercept = 0, lty='dashed')
ggsave("../../../results/partialILRmodelling_ME/differences_in_partialILR_HRDs7stratification_imput.pdf", width = 6, height = 4)

pairs(exposures_inputation1em2_ILR, col=scales::alpha('black', 0.2))


## Running imputation
## ILR
TMB_params_imputation1em2_ILR = TMB_data
TMB_params_imputation1em2_ILR$Y <- exposures_inputation1em2_ILR[keep_britroc_samples,]
give_summary_dimensions_TMBobj(TMB_params_imputation1em2_ILR)
res_nlminb_imputation1em2_ILR <- wrapper_run_TMB_use_nlminb(model = "partialILR",
    object = TMB_params_imputation1em2_ILR, use_nlminb = T, iter.max=10000)
res_nlminb_imputation1em2_ILR
give_pairs_with_mvn_wrapper(matrix(res_nlminb_imputation1em2_ILR$par.random, ncol=TMB_params_imputation1em2_ILR$d), common_lims = F) ## no problem with differences in magnitudes

## ALR
TMB_params_imputation1em2_ALR = TMB_data
TMB_params_imputation1em2_ALR$Y <- exposures_inputation1em2_ALR[keep_britroc_samples,]
give_summary_dimensions_TMBobj(TMB_params_imputation1em2_ALR)
res_nlminb_imputation1em2_ALR <- wrapper_run_TMB_use_nlminb(model = "partialILRnocoroutsidesd",
                                                            object = TMB_params_imputation1em2_ALR, use_nlminb = T, iter.max=10000)
res_nlminb_imputation1em2_ALR
give_pairs_with_mvn_wrapper(matrix(res_nlminb_imputation1em2_ALR$par.random, ncol=TMB_params_imputation1em2_ALR$d), common_lims = F) ## no problem with differences in magnitudes

object_to_cov_heatmap(res_nlminb_imputation1em2_ALR, TMB_params_imputation1em2_ALR$d) ## extremely correlated
object_to_cov_heatmap(res_nlminb_imputation1em2_ILR, TMB_params_imputation1em2_ILR$d,
                      names_cov = paste0('ILR', 1:6)) ## correlated to a lesser extent

pdf("../../../results/partialILRmodelling_ME/res_nlminb_partialILRnocoroutsidesd_imputation1em2_ALR.pdf", width = 6, height = 4)
plot_betas(res_nlminb_imputation1em2_ALR, names_cats = paste0(paste0('ALR(', 1:6), '/7)'))
dev.off()

##------------------------------------------------------------------------------------------------##
## Using the previous partial ILR definition
TMB_obj_exposures_previouspartialILR <- prepare_TMB_data_with_subset(subset_sigs = colnames(exposures),
                                                                     ilr_trans = T, exposures = exposures, previous_implementation = T)
TMB_obj_exposures_currentpartialILR <- prepare_TMB_data_with_subset(subset_sigs = colnames(exposures),
                                                                    ilr_trans = T, exposures = exposures, previous_implementation = F)


res_nlminb_previouspartialILR <- wrapper_run_TMB_use_nlminb(model = "partialILR", object = TMB_obj_exposures_previouspartialILR, use_nlminb = T,
                                                           iter.max=1000)

res_nlminb_currentpartialILR <- wrapper_run_TMB_use_nlminb(model = "partialILR", object = TMB_obj_exposures_currentpartialILR, use_nlminb = T,
                                                       iter.max=1000)

res_nlminb_previouspartialILR$pdHess
res_nlminb_currentpartialILR$pdHess

## comparison of partial ILR
plot(as.vector(TMB_obj_exposures_previouspartialILR$Y),
     as.vector(TMB_obj_exposures_currentpartialILR$Y)) ## extremely similar, for the most part

## comparing to estimates
ggplot(data.frame(prev=res_nlminb_previouspartialILR$par.fixed, cur=res_nlminb_currentpartialILR$par.fixed,
                  names=names(res_nlminb_currentpartialILR$par.fixed)), aes(x=prev, y=cur))+
  geom_point()+facet_wrap(.~names, scales = "free")+theme_bw()+geom_abline(intercept = 0, slope = 1)
## very similar




