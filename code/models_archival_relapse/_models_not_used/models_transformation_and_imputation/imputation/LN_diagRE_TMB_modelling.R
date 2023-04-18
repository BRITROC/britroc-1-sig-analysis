## Infence using stan_fit_LNM for data in the simplex
## Simple LN, no fixed or random effects

rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("../../")

set.seed(1234)

library(rstan)
library(uuid)
library(ggplot2)
library(optparse)
library(reshape2)
library(TMB)
library(pheatmap)
library(bayesplot)
library(dplyr)
library(parallel)
source("../../../../GlobalDA/code/2_inference/helper/helper_DA_stan.R")
source("../../../../GlobalDA/code/3_analysis/helper/helper_analyse_posteriors.R")
source("../../../../GlobalDA/code/2_inference_TMB/helper_TMB.R")
source("helper/functions.R")
Nits = 2500

#-------------------------------------------------------------------------------------------#
TMB::compile("tmb_RE/tmb_RE.cpp", "-std=gnu++17") ##20211105: there were problems with this file
dyn.load(dynlib("tmb_RE/tmb_RE"))
## Using the same model as Alpha trans for now (20211105)
TMB::compile("tmb_RE/tmb_RE_lm.cpp", "-std=gnu++17")
dyn.load(dynlib("tmb_RE/tmb_RE_lm"))

#-------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------#
give_short_names = function(i){
  if(grepl('TCGA', i)){
    paste0(strsplit(i, split = '-')[[1]][1:3], collapse='-')
  }else{
    i
  }
}

load("../../data/britroc_30kb_signature_data.rds")
# patient.meta
# sample_by_component
# sig_quants

exposures = t(sig_quants)
exposures = normalise_rw(exposures + 1e-4)
patient.meta = patient.meta[match(rownames(exposures), as.character(patient.meta$SAMPLE_ID)),]
sample_by_component = sample_by_component[match(rownames(exposures), rownames(sample_by_component)),]
all(rownames(exposures) == patient.meta$SAMPLE_ID)

# ------------------------------------------------------------------------------------- #

## distributions of log-ratios
pairs(as(compositions::alr(exposures), 'matrix'))
pdf("../results/LN_modelling/logratios_inputation_exposures.pdf", width = 10, height = 2)
par(mfrow=c(1,6))
sapply(1:6, function(j) hist(as(compositions::alr(exposures), 'matrix')[,j], breaks=30,
                             main=paste0('Log-ratio ', j), xlab='Value of log-ratio'))
dev.off()

ggplot(melt(as(compositions::alr(exposures), 'matrix')), aes(x=value))+geom_histogram()+facet_wrap(.~Var2, ncol=6)+
  theme_bw()
ggsave("../results/LN_modelling/logratios_inputation_exposures2.pdf",  width = 10, height = 2)

sapply(1)

ggplot(melt(as(compositions::alr(normalise_rw(exposures + 1e-3)), 'matrix')), aes(x=value))+geom_histogram()+facet_wrap(.~Var2, ncol=6)+
  theme_bw()
ggsave("../results/LN_modelling/logratios_inputation_exposures2_1emin3.pdf",  width = 10, height = 2)

par(mfrow=c(2,3))
apply(as(compositions::alr(normalise_rw(t(sig_quants) + 1e-2)), 'matrix'), 2, hist, breaks=30)

#-------------------------------------------------------------------------------------------#

x = t(cbind(1, as.numeric(factor(patient.meta$group))-1))

d = ncol(exposures) ## number of features
n = nrow(exposures) ## number of samples
num_indiv = length(unique(patient.meta$PATIENT_ID))

TMB_data = list(Y = exposures,
                num_individuals = num_indiv,
                x = t(x),
                z = sapply(unique(patient.meta$PATIENT_ID), function(i) as.numeric(patient.meta$PATIENT_ID == i)))

TMB_params = list(beta = matrix(rep(runif(1, min = -4, max = 4), 2*(d-1)),
                                              nrow = 2, byrow=TRUE),
                  u_large = matrix(rep(1, (d-1)*num_indiv), nrow=num_indiv),
                  logs_sd_RE=rep(1, d-1))

## 20211105: I haven't been able to change any of this
obj <- MakeADFun(data = TMB_data, parameters = TMB_params, DLL="tmb_RE", random = "u_large")
# obj <- MakeADFun(data = TMB_data, parameters = TMB_params, DLL="tmb_RE_lm", random = "u_large") ##20211105: changing <tmb_RE> to <tmb_RE_lm>
# obj <- MakeADFun(data = TMB_data, parameters = TMB_params, DLL="tmb_RE" )
opt <- do.call("optim", obj)
opt
opt$hessian ## <-- FD hessian from optim
rep <- sdreport(obj)
rep

beta_est = matrix(python_like_select_name(rep$par.fixed, 'beta'), nrow = 2)
u_est = matrix(python_like_select_name(rep$par.random, 'u_large'), ncol = (d-1))
thetaLR = TMB_data$x %*% beta_est + TMB_data$z %*% u_est

exposures_fitted = softmax(cbind(thetaLR, 0))
dim(exposures_fitted)
dim(exposures)


colnames(exposures_fitted) = colnames(exposures)

df_scatter = cbind.data.frame(exposures=as.vector(exposures), fitted_exposures= as.vector(exposures_fitted),
                  signature=rep(colnames(exposures), each=nrow(exposures)),
                  patient=rep(rownames(exposures), 7))
df_scatter$subtract=df_scatter$fitted_exposures - df_scatter$exposures
df_matrix = cbind.data.frame(sample=rownames(exposures), exposures=(exposures),
                             fitted_exposures= as.vector(exposures_fitted))
ggthemr::ggthemr("fresh")
ggplot(df_scatter,
     aes(x=exposures, y=fitted_exposures, col=signature))+geom_point()+facet_wrap(.~signature, nrow=1)+theme_linedraw()+
  theme(axis.text.x = element_text(angle = 30, hjust = 1))+
  geom_abline(slope = 1, intercept = 0, col='black', )
ggsave("../results/LN_modelling/RELN_TMB_scatter.png", width = 12, height = 2.5)

# levels(df_scatter$patient) = rownames(exposures)[order(df_scatter[df_scatter$signature == 's2','subtract'])]
order_patients = rownames(exposures)[order(df_scatter[df_scatter$signature == 's2','subtract'])]
## are the cases that fail from the same patients?
ggplot(df_scatter, aes(x=factor(patient, levels=order_patients),
# ggplot(df_scatter %>% filter(signature=='s2'), aes(x=factor(patient, levels=order_patients),
                                                   y=subtract, col=signature,
                                                   group=factor(patient, levels=order_patients)))+
  # geom_boxplot()+
  theme_bw()+geom_line(aes(group=signature))+facet_wrap(.~signature, ncol=1)
ggsave("../results/LN_modelling/RELN_TMB_subtract_true_fitted.png", width = 12, height = 12)

pdf("../results/LN_modelling/RELN_TMB_matrices_exposures.pdf", width = 6, height = 7)
par(mfrow=c(1,2))
image(t(exposures), main='True exposures')
image(t(exposures_fitted), main='Fitted exposures')
dev.off()

## Analysing the betas
## intercept
## to get the standard errors of estimates
beta_est[1,]
intercept_est = select_intercept(python_like_select_rownames(summary(rep), 'beta')[,1]) ## estimate
intercept_err = select_intercept(python_like_select_rownames(summary(rep), 'beta')[,2]) ## std err
ggplot(cbind.data.frame(name=paste0('LR', 1:(d-1)), intercept_est, intercept_err),
       aes(x=name, y=intercept_est))+
  geom_errorbar(aes(ymin=intercept_est-intercept_err, ymax=intercept_est+intercept_err), width=.1)+
  geom_line() + geom_point()+theme_linedraw()+labs(x='Name of beta intercept', y='Estimate for beta value')
ggsave("../results/LN_modelling/RELN_TMB_intercept_errorbar.png", width = 3, height = 2.5)

## slope

slope_est = select_slope_2(python_like_select_rownames(summary(rep), 'beta')[,1]) ## estimate
slope_err = select_slope_2(python_like_select_rownames(summary(rep), 'beta')[,2]) ## std err
ggplot(cbind.data.frame(name=paste0('LR', 1:(d-1)), slope_est, slope_err),
       aes(x=name, y=slope_est))+
  geom_errorbar(aes(ymin=slope_est-slope_err, ymax=slope_est+slope_err), width=.1)+
  geom_line() + geom_point()+theme_linedraw()+labs(x='Name of beta slope', y='Estimate for beta value')
ggsave("../results/LN_modelling/RELN_TMB_slope_errorbar.png", width = 3, height = 2.5)


## Analysing the random effects
ggthemr::ggthemr_reset()
RE_est = python_like_select_rownames(summary(rep), 'u_large')[,1] ## estimate
RE_err = python_like_select_rownames(summary(rep), 'u_large')[,2] ## std err
df_RE = cbind.data.frame(name_indiv=factor(rep(1:num_indiv, d-1)),
                         name_LR = factor(paste0('LR', 1:(d-1))[rep(1:(d-1), each=num_indiv)]),
                 RE_est, RE_err)
ggplot(df_RE,
       # aes(x=interaction(name_LR,name_indiv), y=RE_est, col=(name_LR)))+
       # aes(x=name_indiv, group = interaction(name_indiv, name_LR), y=RE_est, col=(name_LR)))+facet_wrap(.~name_indiv)+
       aes( x = name_LR, y=RE_est, col=(name_LR)))+facet_wrap(.~name_indiv, nrow=2)+
  geom_errorbar(aes(ymin=RE_est-RE_err, ymax=RE_est+RE_err), width=.1)+
  geom_point()+labs(x='Name of beta slope', y='Estimate for beta value')
ggsave("../results/LN_modelling/RELN_TMB_RE_coefs_errorbar.png", width = 20, height = 2.5)

ggplot(df_RE,
       # aes( x = factor(name_indiv, levels=rownames(exposures)[order(rowSums(abs(u_est)))]),
       aes( x = factor(name_indiv, levels=order(rowSums(abs(u_est)))),
            y=RE_est, group=name_LR))+facet_wrap(.~name_LR)+
  geom_errorbar(aes(ymin=RE_est-RE_err, ymax=RE_est+RE_err), width=.1)+
  geom_point()+geom_line()+labs(x='Ordered patients', y='Coefficient for random effects')
ggsave("../results/LN_modelling/RELN_TMB_RE_coefs_errorbar2.png", width = 20, height = 4.5)

ggplot(df_RE,
       # aes( x = factor(name_indiv, levels=rownames(exposures)[order(rowSums(abs(u_est)))]),
       aes( x = RE_est, col=name_LR))+facet_wrap(.~name_LR, nrow=1)+
  geom_density()+
  labs(x='Coefficients for random effects', y='Density')
ggsave("../results/LN_modelling/RELN_TMB_RE_coefs_densities.png", width = 10, height = 2.5)

ggplot(df_RE,
       # aes( x = factor(name_indiv, levels=rownames(exposures)[order(rowSums(abs(u_est)))]),
       aes( x = RE_est))+
  geom_density()+
  labs(x='Coefficients for random effects', y='Density')
ggsave("../results/LN_modelling/RELN_TMB_RE_coefs_density_all.png", width = 3, height = 3)
## add the normal from the estimates below!

## suspicious that it seems to be the same gaussian that all random effects coefficients come from
sds = python_like_select_rownames(summary(rep), 'logs_sd_RE')[,1]
names(sds) = paste0('LR', 1:(d-1))
sds = melt(sapply(sds, rnorm, n = 2000, mean=0)); colnames(sds) = c('name_indiv', 'name_LR', 'RE_est')
ggplot(rbind.data.frame(cbind.data.frame(sds, RE_err=NA, sim='Sim'),
                      cbind.data.frame(df_RE, sim='Obs')), aes(x=RE_est, col=sim))+geom_density()+
  facet_wrap(.~name_LR, nrow=1)
ggsave("../results/LN_modelling/RELN_TMB_RE_coefs_density_all_with_sim.png", width = 8, height = 2)

## So, is it differentially abundant?
## differential abundance
wald_TMB_wrapper(rep, verbatim = FALSE)

