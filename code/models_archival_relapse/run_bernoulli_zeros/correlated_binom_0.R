rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source("../helper/header.R")

exposures_transformed = t(sapply(as.vector(exposures_zeros), function(i) table(factor(i, levels=c(0,1)))))
d = ncol(exposures_zeros)

## 20220314: 1 is now zero, and 0 is non-zero (to unify)
exposures_zeros <- 1-exposures_zeros
colSums(exposures_zeros) ## e.g. very few zeros at s1

## How many zeros are there in each group?
groups_numeric = as.numeric(as.factor(patient.meta$group))
apply(exposures_zeros[patient.meta$group == 'arx',], 2, table)[1,]/sum(patient.meta$group == 'arx')
apply(exposures_zeros[patient.meta$group == 'rlps',], 2, table)[1,]/sum(patient.meta$group == 'rlps')

stopifnot(all(rownames(exposures)== patient.meta$SAMPLE_ID))

TMB_data = list(Y = exposures_zeros,
                num_individuals = num_indiv,
                num_sigs = ncol(exposures),
                x = cbind(1, as.numeric(as.factor(patient.meta$group))-1),
                z = give_z_matrix_from_labels(patient.meta$PATIENT_ID))

# give_res_report = function(){
#   TMB_params = list(beta = (matrix(runif(d*2, min = -4, max = 4),
#                                    nrow = 2, byrow=TRUE)),
#                     u_large = matrix(rep(1, (d)*(TMB_data$num_individuals)), nrow=TMB_data$num_individuals),
#                     logs_sd_RE=runif(n = d, min = 0, max = 2)
#   )
#   
#   
#   obj <- MakeADFun(data = TMB_data, parameters = TMB_params, DLL="tmb_correlated_multinom", random = "u_large")
#   obj$hessian <- TRUE
#   opt <- do.call("optim", obj)
#   opt
#   opt$hessian ## <-- FD hessian from optim
#   report = sdreport(obj)
#   
#   return(list(TMB_params=TMB_params, report=report))
#   
# }

results = wrapper_run_TMB_use_nlminb(model = 'bernoulliMEnocor', object = TMB_data)
results <- list(report=results, TMB_params=NULL)
plot_betas(results$report)

# results = give_res_report()
results
report = results$report
TMB_params = results$TMB_params

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
  geom_raster( aes( x = Var2, y = Var1, fill = value )  )+
  scale_fill_gradient(low="blue", high="red")+
  geom_point(data=melt((TMB_data$Y)) %>% filter(value == 0),
             aes( x = Var2, y = Var1 ), col='white', shape=18, size=.3)+
  labs(x='Signature', y='Sample')+
  ggtitle('Probabilities of zero exposure and observed zero exposures')
# ggsave("../results/zeros_modelling/correlated_binomial_0.pdf", width = 4, height = 4)

rownames(TMB_data$Y) <- rownames(fitted_probs) <- rownames(exposures)

fitted_zeros <- ggplot(melt(fitted_probs[order(patient.meta$group),]))+
  geom_tile( aes( x = Var2, y = Var1, fill = value )  )+
  scale_fill_gradient(low="blue", high="red")+
  geom_tile(data = reshape2::melt(TMB_data$Y[order(patient.meta$group),])  %>% filter(value == 1), aes(x = Var2, y = Var1),
            fill='white', width=.3)+
  labs(x='Signatures' , y='Samples', fill='Probability')+
  theme(legend.position = "bottom", axis.title.y=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank())+labs(col='Probability of zero')
fitted_zeros
ggsave("../../../results/zeros_modelling/correlated_binomial_0_fitted_zeros.pdf", width = 2.5, height = 3)


##' it only converges if I use no correlations between signatures, i.e. a diagonal covariance matrix. A full
##' covariance matrix doesn't work

results
ggplot(cbind(fitted=melt(fitted_probs), true=melt((TMB_data$Y))), aes(x=fitted.value, y=true.value))+geom_point(aes(col=fitted.value>0.5))+
  facet_wrap(.~fitted.Var2)


uuid = UUIDgenerate()
save.image(paste0("../../../out/inference/correlated_binom_0_", uuid, '.Data'))
saveRDS(results, file = paste0("../../../out/inference/correlated_binom_0_", uuid, '.RDS'))

##--------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------##
## Only including samples for patients with both groups
both_arx_rlps <- sapply(unique(patient.meta$PATIENT_ID), function(i) all(c('rlps', 'arx') %in% 
                                                                           patient.meta$group[patient.meta$PATIENT_ID == i]))
both_arx_rlps <- rownames(exposures) %in% patient.meta$SAMPLE_ID[patient.meta$PATIENT_ID %in% patient.meta$PATIENT_ID[both_arx_rlps]]

table(TMB_data$Y)
table(TMB_data_only_matched$Y)

# TMB_data_only_matched = give_subset_samples_TMBobj(TMB_data,samples_to_remove = which(!both_arx_rlps))
TMB_data_only_matched = give_subset_samples_TMBobj(TMB_data,samples_to_remove = as.character(patient.meta$SAMPLE_ID[!both_arx_rlps]))
dim(TMB_data$Y)
dim(TMB_data_only_matched$Y)
TMB_data$Y['JBLAB-4979',]
TMB_data_only_matched$Y['JBLAB-4979',]
rownames(TMB_data_only_matched$Y) <- rownames(exposures)[both_arx_rlps]

TMB_data_only_matched$Y
results_only_matched = wrapper_run_TMB_use_nlminb(model = 'bernoulliMEnocor', object = TMB_data_only_matched)
results_only_matched <- list(report=results_only_matched, TMB_params=TMB_data_only_matched)
plot_betas(results_only_matched$report)
saveRDS(results_only_matched, file = paste0("../../../out/inference/correlated_binom_0_", uuid, '_only_matched.RDS'))

fitted_logR_only_matched = simulate_from_correlated_binom(tmb_fit_object = results_only_matched$report, full_RE = T,
                                                          x_matrix=results_only_matched$TMB_params$x,
                                                          z_matrix=results_only_matched$TMB_params$z, return_logratios=T)
fitted_probs_only_matched = simulate_from_correlated_binom(tmb_fit_object = results_only_matched$report, full_RE = T,
                                                           x_matrix=results_only_matched$TMB_params$x, z_matrix=results_only_matched$TMB_params$z, return_logratios=F)
# _only_matched
colnames(fitted_probs_only_matched) <- colnames(results_only_matched$TMB_params$Y)

rownames(results_only_matched$TMB_params$Y) <- rownames(fitted_probs_only_matched) <- rownames(exposures)[both_arx_rlps]
fitted_zeros_only_matched <- ggplot(melt(fitted_probs_only_matched[order(patient.meta$group[both_arx_rlps]),]))+
  geom_tile( aes( x = Var2, y = Var1, fill = value )  )+
  scale_fill_gradient(low="blue", high="red")+
  geom_tile(data = reshape2::melt(results_only_matched$TMB_params$Y[order(patient.meta$group[both_arx_rlps]),])  %>% 
              filter(value == 1), aes(x = Var2, y = Var1),
            fill='white', width=.3)+
  labs(x='Signatures' , y='Samples', fill='Probability')+
  theme(legend.position = "bottom", axis.title.y=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank())+labs(col='Probability of zero')
fitted_zeros_only_matched
ggsave("../../../results/zeros_modelling/correlated_binomial_0_fitted_zeros_only_matched.pdf", width = 2.5, height = 3)


ggplot()+
  geom_tile(data = reshape2::melt(results_only_matched$TMB_params$Y[order(patient.meta$group[both_arx_rlps]),])  %>% filter(value == 1), aes(x = Var2, y = Var1),
            fill='white', width=.3)
##--------------------------------------------------------------------------------##


##--------------------------------------------------------------------------------##
## read partialILR
results <- readRDS(file = paste0("../../../out/inference/correlated_binom_0_63ddb8a4-d766-4b28-bb2d-3d4b0c2ebf62.RDS"))

# partialILRnocor_nlminb_allsigs <- readRDS(paste0("../../../out/inference/partialILR/partialILRnocor_nlminb_allsigs.RDS"))
# partialILRnocor_nlminb_allsigs_TMBdata <- readRDS(paste0("../../../out/inference/partialILR/partialILRnocor_nlminb_allsigs_TMBdata.RDS"))
partialILRnocor_nlminb_allsigs_only_matched <- readRDS(paste0(folder_out_RDS, "res_nlminb_nocoroutsidesd_only_matched_allsigs.RDS"))
partialILRnocor_nlminb_allsigs_TMBdata_only_matched <- readRDS(paste0("../../../out/inference/partialILR/TMB_data_only_matched_res_nlminb_nocoroutsidesd_only_matched_allsigs_RData.RDS"))


ggplot(cbind.data.frame(label=colnames(exposures_zeros),
                        col=(colnames(exposures_zeros) == 's1'), ## s1 is used as baseline, somewhat
                        partialILR=rbind(0, python_like_select_rownames(summary(partialILRnocor_nlminb_allsigs), 'beta')[c(F,T),])[match(colnames(exposures_zeros), colnames(exposures)),],
                        zeros=python_like_select_rownames(summary(results$report), 'beta')[c(F,T),]),
       aes(x=partialILR.Estimate, y=-zeros.Estimate, label=label, col=col))+geom_point()+
  geom_hline(yintercept = 0, col='blue', lty='dashed', size=0.3)+
  geom_vline(xintercept = 0, col='blue', lty='dashed', size=0.3)+
  geom_errorbar(aes(xmin=`partialILR.Estimate`-`partialILR.Std. Error`,xmax=`partialILR.Estimate`+`partialILR.Std. Error`), width=.1)+
  geom_errorbar(aes(ymin=-`zeros.Estimate`-`zeros.Std. Error`, ymax=-`zeros.Estimate`+`zeros.Std. Error`), width=.1)+
  # geom_label_repel(aes(col=col))+
  facet_wrap(.~label, nrow=1, scales = "free")+guides(col='none')+
  labs(x='Partial ILR estimate', y='Bernoulli estimate')+
  # ggtitle('Summary of both models')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("../../../results/zeros_modelling/britroc_summary_both_sets_coefs_nocor.pdf", width = 8, height = 1.8)

ggplot(cbind.data.frame(label=colnames(exposures_zeros[both_arx_rlps,]),
                        col=(colnames(exposures_zeros[both_arx_rlps,]) == 's1'), ## s1 is used as baseline, somewhat
                        partialILR=rbind(0, python_like_select_rownames(summary(partialILRnocor_nlminb_allsigs_only_matched), 'beta')[c(F,T),])[match(colnames(exposures_zeros[both_arx_rlps,]), colnames(exposures[both_arx_rlps,])),],
                        zeros=python_like_select_rownames(summary(results_only_matched$report), 'beta')[c(F,T),]),
       aes(x=partialILR.Estimate, y=-zeros.Estimate, label=label, col=col))+geom_point()+
  geom_hline(yintercept = 0, col='blue', lty='dashed', size=0.3)+
  geom_vline(xintercept = 0, col='blue', lty='dashed', size=0.3)+
  geom_errorbar(aes(xmin=`partialILR.Estimate`-`partialILR.Std. Error`,xmax=`partialILR.Estimate`+`partialILR.Std. Error`), width=.1)+
  geom_errorbar(aes(ymin=-`zeros.Estimate`-`zeros.Std. Error`, ymax=-`zeros.Estimate`+`zeros.Std. Error`), width=.1)+
  facet_wrap(.~label, nrow=1, scales = "free")+guides(col='none')+
  labs(x='Partial ILR estimate', y='Bernoulli estimate')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("../../../results/zeros_modelling/britroc_summary_both_sets_coefs_nocor_only_matched.pdf", width = 8, height = 1.8)

ggplot(cbind.data.frame(label=colnames(exposures_zeros[both_arx_rlps,]),
                        col=(colnames(exposures_zeros[both_arx_rlps,]) == 's1'), ## s1 is used as baseline, somewhat
                        partialILR=rbind(0, python_like_select_rownames(summary(partialILRnocor_nlminb_allsigs_only_matched), 'beta')[c(F,T),])[match(colnames(exposures_zeros[both_arx_rlps,]), colnames(exposures[both_arx_rlps,])),],
                        zeros=python_like_select_rownames(summary(results_only_matched$report), 'beta')[c(F,T),]),
       aes(x=partialILR.Estimate, y=-zeros.Estimate, label=label, col=col))+geom_point()+
  geom_hline(yintercept = 0, col='blue', lty='dashed', size=0.3)+
  geom_vline(xintercept = 0, col='blue', lty='dashed', size=0.3)+
  geom_errorbar(aes(xmin=`partialILR.Estimate`-`partialILR.Std. Error`,xmax=`partialILR.Estimate`+`partialILR.Std. Error`), width=.1)+
  geom_errorbar(aes(ymin=-`zeros.Estimate`-`zeros.Std. Error`, ymax=-`zeros.Estimate`+`zeros.Std. Error`), width=.1)+
  facet_wrap(.~label, nrow=1)+guides(col='none')+
  labs(x='Partial ILR estimate', y='Bernoulli estimate')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("../../../results/zeros_modelling/britroc_summary_both_sets_coefs_nocor_only_matched_sharedaxes.pdf", width = 8, height = 1.8)


##--------------------------------------------------------------------------------##

