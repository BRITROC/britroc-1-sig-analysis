## Infence using stan_fit_LNM for data in the simplex
## Simple LN, no fixed or random effects

rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# set.seed(1234)

library(uuid)
library(ggplot2)
library(reshape2)
library(compositions)
# source("../../../GlobalDA/code/2_inference/helper/helper_DA_stan.R")
# source("../../../GlobalDA/code/3_analysis/helper/helper_analyse_posteriors.R")
source("../../../GlobalDA/code/2_inference_TMB/helper_TMB.R")
library(TMB)

source("../../../GlobalDA/code/2_inference_TMB/helper_TMB.R")
source("functions.R")
source("header.R")

#-------------------------------------------------------------------------------------------#
give_short_names = function(i){
  if(grepl('TCGA', i)){
    paste0(strsplit(i, split = '-')[[1]][1:3], collapse='-')
  }else{
    i
  }
}

give_TMB_params = function(arg_d, arg_num_individuals){
  list(beta = (matrix(runif(arg_d*2, min = -4, max = 4),
                      nrow = 2, byrow=TRUE)),
       u_large = matrix(rep(1, (arg_d)*(arg_num_individuals)), nrow=arg_num_individuals),
       logs_sd_RE=runif(n = arg_d, min = 0, max = 2),
       cov_RE = runif(n = ((arg_d)*(arg_d)-(arg_d))/2, min = 0.1, max = 0.2))
}

load("../data/britroc_30kb_signature_data.rds")
# patient.meta
# sample_by_component
# sig_quants

exposures = t(sig_quants)

## random order
patient.meta <- patient.meta[sample(1:nrow(patient.meta), size = nrow(patient.meta), replace = F),]
# patient.meta <- patient.meta[match(rownames(exposures),
#                                    patient.meta$SAMPLE_ID),]

all(patient.meta$SAMPLE_ID == rownames(exposures))

#-------------------------------------------------------------------------------------------#

## Transform data with partial ILR
which_zero = t(apply(exposures, 1, function(i) as.numeric((i==0)) ))


exposures[1,]
which_zero[1,]
as.vector(compositions::ilr(exposures)[1,])
give_partial_ilr_basis(which(which_zero[1,] == 1), d=7)
give_partial_ilr_basis(which_zero_vector = 2, d = 5)


give_partial_ilr_basis(which_zero_vector = 2, d = 4)

prepare_TMB_data_with_subset = function(subset_sigs){
  .exps = normalise_rw(exposures[,subset_sigs])
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

exposures_partial_irl = give_partial_irl(exposures)

dim(exposures)
dim(exposures_partial_irl)

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
TMB::compile("tmb_RE/tmb_MVN_partial_ILR.cpp", "-std=gnu++17")
dyn.load(dynlib("tmb_RE/tmb_MVN_partial_ILR"))
#-------------------------------------------------------------------------------------------#


#-------------------------------------------------------------------------------------------#
## remove the samples in which there are too many zeros. There can be at most 4 zeros (i.e. there has to be at least one log-ratio)
keep_britroc_samples = (rowSums(exposures_partial_irl == 0) < (ncol(exposures) - 2) )
exposures_partial_irl = exposures_partial_irl[keep_britroc_samples,]
all(patient.meta$SAMPLE_ID[keep_britroc_samples] == rownames(exposures_partial_irl))
z_britroc = give_z_matrix_from_labels(patient.meta$PATIENT_ID[keep_britroc_samples])
TMB_data = list(Y = exposures_partial_irl,
                num_individuals = ncol(z_britroc),
                d = d,
                n = nrow(exposures_partial_irl),
                x = cbind(1, as.numeric(as.factor(patient.meta$group[keep_britroc_samples]))-1),
                z = z_britroc)

sapply(TMB_data, dim)
sapply(TMB_params, dim)
#-------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------#
## TRUE_data: full
TMB_params = give_TMB_params(d, TMB_data$num_individuals)
obj <- MakeADFun(data = TMB_data, parameters = TMB_params, DLL="tmb_MVN_partial_ILR", random = "u_large")
opt <- do.call("optim", obj)
opt
opt$hessian ## <-- FD hessian from optim
rep <- sdreport(obj)
rep

repeat_full_with_different_initial = function(){
  TMB_params = give_TMB_params(d, TMB_data$num_individuals)
  obj <- MakeADFun(data = TMB_data, parameters = TMB_params, DLL="tmb_MVN_partial_ILR", random = "u_large")
  opt <- do.call("optim", obj)
  opt
  opt$hessian ## <-- FD hessian from optim
  rep <- sdreport(obj)
  rep
}

repeat_full_with_different_initial_obs = replicate(n = 10, repeat_full_with_different_initial())
# saveRDS(repeat_full_with_different_initial_obs, "../out/inference/partialILR/repeat_full_with_different_initial_obs.RDS")
system("mkdir ../results/UNSORTED_DEBUG_partialILRmodelling/")

# repeat_full_with_different_initial_obs = readRDS("../out/inference/partialILR/repeat_full_with_different_initial_obs.RDS")
## we have one pdHess true!!
unlist(apply(repeat_full_with_different_initial_obs, 2, `[`, 'pdHess'))

repeat_full_with_different_initial_obs[,5]$pdHess

## Beta values are the same regardless of whether the runs have converged or not
pdf("../results/UNSORTED_DEBUG_partialILRmodelling/full_partialILR_several_runs.pdf")
pairs(repeat_full_with_different_initial_obs['par.fixed',])
dev.off()
pdf("../results/UNSORTED_DEBUG_partialILRmodelling/full_partialILR_several_runs_only_beta.pdf")
pairs(sapply(repeat_full_with_different_initial_obs['par.fixed',], python_like_select_name, grep_substring="beta"),
      main='Only beta (intercept+slope)')
dev.off()

## Analyse the results of the converged run
results_full0 = repeat_full_with_different_initial_obs[,3]
results_full = TMB::summary.sdreport(results_full0)
intercept_est = select_intercept(python_like_select_rownames(results_full, 'beta')[,1]) ## estimate
intercept_err = select_intercept(python_like_select_rownames(results_full, 'beta')[,2]) ## std err
ggplot(cbind.data.frame(name=paste0('ILR', 1:d), intercept_est, intercept_err),
       aes(x=name, y=intercept_est))+
  geom_errorbar(aes(ymin=intercept_est-intercept_err, ymax=intercept_est+intercept_err), width=.1)+
  geom_line() + geom_point()+theme_linedraw()+labs(x='Name of beta intercept', y='Estimate for beta value')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave("../results/UNSORTED_DEBUG_partialILRmodelling/full_partialILR_intercept_errorbar.png", width = 3, height = 2.5)

slope_est = select_slope_2(python_like_select_rownames(results_full, 'beta')[,1],v=F) ## estimate
slope_err = select_slope_2(python_like_select_rownames(results_full, 'beta')[,2],v=F) ## std err
ggplot(cbind.data.frame(name=paste0('ILR', 1:d), slope_est, slope_err),
       aes(x=name, y=slope_est))+
  geom_errorbar(aes(ymin=slope_est-slope_err, ymax=slope_est+slope_err), width=.1)+
  geom_line() + geom_point()+theme_linedraw()+labs(x='Name of beta slope', y='Estimate for beta value')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave("../results/UNSORTED_DEBUG_partialILRmodelling/full_partialILR_slope_errorbar.png", width = 3, height = 2.5)

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
# saveRDS(rep, "../out/UNSORTED_DEBUG_inference/partialILR/repeat_full_nos5.RDS")
wald_TMB_wrapper(rep)

ggplot(cbind.data.frame(name=paste0('ILR', c(1:3, 5:6)),
                        intercept_est=select_intercept(python_like_select_rownames(TMB::summary.sdreport(rep), 'beta')[,1]),
                        intercept_err=select_intercept(python_like_select_rownames(TMB::summary.sdreport(rep), 'beta')[,2])),
       aes(x=name, y=intercept_est))+
  geom_errorbar(aes(ymin=intercept_est-intercept_err, ymax=intercept_est+intercept_err), width=.1)+
  geom_line() + geom_point()+theme_linedraw()+labs(x='Name of beta intercept', y='Estimate for beta value')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave("../results/UNSORTED_DEBUG_partialILRmodelling/full_partialILR_nos5_intercept_errorbar.png", width = 3, height = 2.5)
ggplot(cbind.data.frame(name=paste0('ILR', c(1:3, 5:6)),
                        slope_est=select_slope_2(python_like_select_rownames(TMB::summary.sdreport(rep), 'beta')[,1],v=F),
                        slope_err=select_slope_2(python_like_select_rownames(TMB::summary.sdreport(rep), 'beta')[,2],v=F)),
       aes(x=name, y=slope_est))+
  geom_errorbar(aes(ymin=slope_est-slope_err, ymax=slope_est+slope_err), width=.1)+
  geom_line() + geom_point()+theme_linedraw()+labs(x='Name of beta slope', y='Estimate for beta value')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave("../results/UNSORTED_DEBUG_partialILRmodelling/full_partialILR_nos5_slope_errorbar.png", width = 3, height = 2.5)
#-------------------------------------------------------------------------------------------#
