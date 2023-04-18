rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source("header.R")
library(gridExtra)
library(ggrepel)

correlated_binom_0_fit = readRDS(file = paste0("../../../out/inference/correlated_binom_0_251afb3a-a218-475b-950c-7f1503acc5e1.RDS"))
correlated_binom_1_fit = readRDS(file = paste0("../../../out/inference/correlated_binom_1_7263302b-a1ea-4010-ba3d-415e655afc3e.RDS"))
correlated_binom_2_fit= readRDS(file = paste0("../../../out/inference/correlated_binom_2_a8acbdff-a73e-4bdf-8525-f72803583e2d.RDS"))


give_fitted_probs = function(results_TMB, TMB_data){
  fitted_logR = simulate_from_correlated_binom(tmb_fit_object = results_TMB, full_RE = T, x_matrix=TMB_data$x, z_matrix=TMB_data$z, return_logratios=T)
  fitted_probs = simulate_from_correlated_binom(tmb_fit_object = results_TMB, full_RE = T, x_matrix=TMB_data$x, z_matrix=TMB_data$z, return_logratios=F)
  colnames(fitted_probs) = colnames(TMB_data$Y)
  
}

## betas
plot(python_like_select_name(correlated_binom_1_fit$report$par.fixed, 'beta'),
     python_like_select_name(correlated_binom_2_fit$report$par.fixed, 'beta'))
abline(coef = c(0,1))

pairs(cbind(indepRE=python_like_select_name(correlated_binom_0_fit$report$par.fixed, 'beta'),
      partialRE=python_like_select_name(correlated_binom_1_fit$report$par.fixed, 'beta'),
      fullRE=python_like_select_name(correlated_binom_2_fit$report$par.fixed, 'beta')))

df_slopes_1_0 = cbind.data.frame(partialRE=python_like_select_name(correlated_binom_1_fit$report$par.fixed, 'beta'),
                 indepRE=python_like_select_name(correlated_binom_0_fit$report$par.fixed, 'beta'),
                 label=rep(paste0('s', 1:7), each=2), type_beta=rep(c('intercept', 'slope')))

df_slopes_2_1 = cbind.data.frame(partialRE=python_like_select_name(correlated_binom_1_fit$report$par.fixed, 'beta'),
                 fullRE=python_like_select_name(correlated_binom_2_fit$report$par.fixed, 'beta'),
                 label=rep(paste0('s', 1:7), each=2),type_beta=rep(c('intercept', 'slope')))

pdf("../results/zeros_modelling/correlated_binomials_betas.pdf")
do.call('grid.arrange', list(
  ggplot(df_slopes_1_0,
         aes(x=indepRE, y=partialRE, label=label, col=type_beta))+geom_point()+facet_wrap(.~type_beta)+
    geom_abline(intercept = 0, slope = 1)+geom_label_repel(),
  ggplot(df_slopes_2_1,
         aes(x=partialRE, y=fullRE, label=label, col=type_beta))+ geom_abline(intercept = 0, slope = 1)+
    facet_wrap(.~type_beta)+geom_point()+geom_label_repel()))
dev.off()

pdf("../results/zeros_modelling/correlated_binomials_betas.pdf")
do.call('grid.arrange', list(
  ggplot(df_slopes,
         aes(x=indepRE, y=partialRE, label=label, col=type_beta))+geom_point()+facet_wrap(.~type_beta)+
    geom_abline(intercept = 0, slope = 1)+geom_label_repel(),
  ggplot(cbind.data.frame(partialRE=python_like_select_name(correlated_binom_1_fit$report$par.fixed, 'beta'),
                          fullRE=python_like_select_name(correlated_binom_2_fit$report$par.fixed, 'beta'),
                          label=paste0('s', 1:7),type_beta=rep(c('intercept', 'slope'))),
         aes(x=partialRE, y=fullRE, label=label, col=type_beta))+ geom_abline(intercept = 0, slope = 1)+
    facet_wrap(.~type_beta)+geom_point()+geom_label_repel()))
dev.off()

## random effects standard deviations
plot(python_like_select_name(correlated_binom_1_fit$report$par.fixed, 'logs_sd_RE'),
     python_like_select_name(correlated_binom_2_fit$report$par.fixed, 'logs_sd_RE'))
abline(c(0,1))

pairs(cbind(indepRE=python_like_select_name(correlated_binom_0_fit$report$par.fixed, 'logs_sd_RE'),
            partialRE=python_like_select_name(correlated_binom_1_fit$report$par.fixed, 'logs_sd_RE'),
            fullRE=python_like_select_name(correlated_binom_2_fit$report$par.fixed, 'logs_sd_RE')))

ggplot(cbind.data.frame(partialRE=python_like_select_name(correlated_binom_1_fit$report$par.fixed, 'logs_sd_RE'),
                        indepRE=python_like_select_name(correlated_binom_0_fit$report$par.fixed, 'logs_sd_RE'),
                        label=paste0('s', 1:7)),
       aes(x=indepRE, y=partialRE, label=label))+geom_point()+geom_label()+
  geom_abline(intercept = 0, slope = 1)

ggplot(cbind.data.frame(partialRE=python_like_select_name(correlated_binom_1_fit$report$par.fixed, 'logs_sd_RE'),
       fullRE=python_like_select_name(correlated_binom_2_fit$report$par.fixed, 'logs_sd_RE'),
       label=paste0('s', 1:7)),
       aes(x=partialRE, y=fullRE, label=label))+geom_point()+geom_label()+
  geom_abline(intercept = 0, slope = 1)

## What is the covariance matrix like?
fitted_cov_mat_sim_2 = give_UNSTRUCTURED_CORR_t_matrix(python_like_select_name(correlated_binom_2_fit$report$par.fixed, 'cov_RE'),
                                                     7)
diag(fitted_cov_mat_sim_2) = exp(python_like_select_name(correlated_binom_2_fit$report$par.fixed, 'logs_sd_RE'))

fitted_cov_mat_sim_0 = matrix(0, 7, 7)
diag(fitted_cov_mat_sim_0) = exp(python_like_select_name(correlated_binom_0_fit$report$par.fixed, 'logs_sd_RE'))
colnames(fitted_cov_mat_sim_2) = colnames(fitted_cov_mat_sim_0) = paste0('s', 1:7)
rownames(fitted_cov_mat_sim_2) = rownames(fitted_cov_mat_sim_0) = paste0('s', 1:7)

pdf("../results/zeros_modelling/correlated_binomial_2_covmat.pdf")
ComplexHeatmap::Heatmap(fitted_cov_mat_sim_2, cluster_rows = FALSE, cluster_columns = FALSE)
dev.off()
pdf("../results/zeros_modelling/correlated_binomial_0_covmat.pdf")
ComplexHeatmap::Heatmap(fitted_cov_mat_sim_0, cluster_rows = FALSE, cluster_columns = FALSE)
dev.off()

## Is there differential abundance?
wald_TMB_wrapper(correlated_binom_0_fit$report)
wald_TMB_wrapper(correlated_binom_2_fit$report)
i=correlated_binom_2_fit$report; verbatim=F
idx_beta = select_slope_2(which(names(i$par.fixed) == "beta"), verbatim=verbatim)
wald_generalised(v = i$par.fixed[idx_beta], sigma = i$cov.fixed[idx_beta,idx_beta])

correlated_binom_0_fit$report$cov.fixed[idx_beta,idx_beta]
correlated_binom_1_fit$report$cov.fixed[idx_beta,idx_beta]
correlated_binom_2_fit$report$cov.fixed[idx_beta,idx_beta]
fitted_cov_mat_sim_2

ComplexHeatmap::Heatmap(correlated_binom_2_fit$report$cov.fixed, cluster_rows = FALSE, cluster_columns = FALSE)
ComplexHeatmap::Heatmap(correlated_binom_2_fit$report$cov.fixed[idx_beta,idx_beta], cluster_rows = FALSE, cluster_columns = FALSE)

pdf("../results/zeros_modelling/correlated_binomial_2_covmat_2.pdf")
ComplexHeatmap::Heatmap(correlated_binom_2_fit$report$cov.fixed[idx_beta,idx_beta], cluster_rows = FALSE, cluster_columns = FALSE)
dev.off()

betas_and_stderr_2 = python_like_select_rownames(summary(correlated_binom_2_fit$report), 'beta')[c(F,T),]
conf_int_2 = give_confidence_interval(betas_and_stderr_2[,1],
                                      betas_and_stderr_2[,2])

betas_and_stderr_0 = python_like_select_rownames(summary(correlated_binom_0_fit$report), 'beta')[c(F,T),]
conf_int_0 = give_confidence_interval(betas_and_stderr_0[,1],
                                      betas_and_stderr_0[,2])


ggplot(cbind.data.frame(betas_and_stderr_2, confint_lower=(conf_int_2[1,]), confint_upper=conf_int_2[2,], idx_within_dataset=1:7,
                        idx=1:7),
       aes(x=idx_within_dataset, y=(Estimate)))+
  # geom_point(data = df_beta_recovery[!(df_beta_recovery$DA_bool),], aes(x=idx_within_dataset, y=beta_true), shape=4)+
  geom_abline(slope = 0, intercept = 0, alpha=0.2)+
  geom_point()+
  geom_errorbar(aes(ymin=Estimate-`Std. Error`, ymax=Estimate+`Std. Error`),  col='blue', width=.2,
                position=position_dodge(.9))+theme_bw()+theme(legend.position = "bottom")+
  geom_errorbar(aes(ymin=confint_lower, ymax=confint_upper), width=.2,
                position=position_dodge(.9))+theme_bw()+theme(legend.position = "bottom")+
  labs(x='Signature', y='Beta slope estimate')
ggsave("../results/zeros_modelling/correlated_binomial_2_covmat_confint_betaslope.pdf", height = 4)

ggplot(cbind.data.frame(betas_and_stderr_0, confint_lower=(conf_int_0[1,]), confint_upper=conf_int_0[2,], idx_within_dataset=1:7,
                        idx=1:7),
       aes(x=idx_within_dataset, y=(Estimate)))+
  # geom_point(data = df_beta_recovery[!(df_beta_recovery$DA_bool),], aes(x=idx_within_dataset, y=beta_true), shape=4)+
  geom_abline(slope = 0, intercept = 0, alpha=0.2)+
  geom_point()+
  geom_errorbar(aes(ymin=Estimate-`Std. Error`, ymax=Estimate+`Std. Error`), col='blue', width=.2,
                position=position_dodge(.9))+theme_bw()+theme(legend.position = "bottom")+
  geom_errorbar(aes(ymin=confint_lower, ymax=confint_upper), width=.2,
                position=position_dodge(.9))+theme_bw()+theme(legend.position = "bottom")+
  labs(x='Signature', y='Beta slope estimate')
ggsave("../results/zeros_modelling/correlated_binomial_0_covmat_confint_betaslope.pdf", height = 4)


select_slope_2(i$par.fixed[which(names(i$par.fixed) == "beta")], verbatim=verbatim)
betas_and_stderr_2


fitted_2 = give_fitted_probs(correlated_binom_2_fit)
