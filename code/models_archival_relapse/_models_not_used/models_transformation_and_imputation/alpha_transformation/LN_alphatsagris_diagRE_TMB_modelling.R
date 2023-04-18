## Infence using stan_fit_LNM for data in the simplex
## Simple LN, no fixed or random effects

rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("../../helper/")
set.seed(1234)

source("header.R")

all(colnames(sig_quants) == patient.meta$SAMPLE_ID)

d
dim(exposures_alphatrans)

TMB_data = list(Y = exposures_alphatrans,
                num_individuals = num_indiv,
                x = t(x),
                z = sapply(unique(patient.meta$PATIENT_ID), function(i) as.numeric(patient.meta$PATIENT_ID == i)))

TMB_params = list(beta = matrix(runif( 2*(d), min = -4, max = 4),
                                nrow = 2, byrow=TRUE),
                  u_large = matrix(rep(1, (d)*num_indiv), nrow=num_indiv),
                  logs_sd_RE=rep(1, d),
                  logs_sd_FE=rep(1, d))

obj <- MakeADFun(data = TMB_data, parameters = TMB_params, DLL="tmb_RE_lm", random = "u_large")
# obj <- MakeADFun(data = TMB_data, parameters = TMB_params, DLL="tmb_RE" )
opt <- do.call("optim", obj)
opt
opt$hessian ## <-- FD hessian from optim
rep <- sdreport(obj)
rep

# lme4::glmer(exposures_alphatrans + x[2,] + (1|patient.meta$PATIENT_ID))


beta_est = matrix(python_like_select_name(rep$par.fixed, 'beta'), nrow = 2)
u_est = matrix(python_like_select_name(rep$par.random, 'u_large'), ncol = (d))
thetaLR = TMB_data$x %*% beta_est + TMB_data$z %*% u_est

## this below is not as easy as with the softmax! I need to find the inverse function
inverse_alpha_trans(thetaLR[1,], alpha = alpha_val)
exposures_fitted = t(apply(thetaLR, 1, inverse_alpha_trans, alpha=alpha_val))
dim(exposures_fitted)
dim(exposures_alphatrans)


colnames(exposures_fitted) = colnames(exposures_alphatrans)

df_scatter = cbind.data.frame(exposures=as.vector(exposures_alphatrans),
                              fitted_exposures= as.vector(thetaLR),
                              signature=rep(colnames(exposures_alphatrans),
                                            each=nrow(exposures_alphatrans)))

df_matrix = cbind.data.frame(sample=rownames(exposures_alphatrans), exposures=(exposures_alphatrans),
                             fitted_exposures= as.vector(exposures_fitted))

# df_scatter = cbind.data.frame(exposures=as.vector(exposures_alphatrans),
#                               fitted_exposures= as.vector(exposures_fitted),
#                               signature=rep(colnames(exposures_alphatrans),
#                                             each=nrow(exposures_alphatrans)))
# df_matrix = cbind.data.frame(sample=rownames(exposures), exposures=(exposures),
#                              fitted_exposures= as.vector(exposures_fitted))


ggthemr::ggthemr("fresh")
ggplot(df_scatter,
       aes(x=exposures, y=fitted_exposures, col=signature))+geom_point()+facet_wrap(.~signature, nrow=1)+theme_linedraw()+
  theme(axis.text.x = element_text(angle = 30, hjust = 1))+
  geom_abline(slope = 1, intercept = 0, col='black', )
ggsave("../results/LN_modelling/RELN_TMB_alphaTsagris_scatter.png", width = 12, height = 2.5)

pdf("../results/LN_modelling/RELN_TMB_alphaTsagris_matrices_exposures.pdf", width = 6, height = 7)
par(mfrow=c(1,2))
image(t(exposures_alphatrans), main='True transformed exposures')
image(t(thetaLR), main='Fitted transformed exposures')
dev.off()

## Analysing the betas
## intercept
## to get the standard errors of estimates
beta_est[1,]
intercept_est = select_intercept(python_like_select_rownames(summary(rep), 'beta')[,1]) ## estimate
intercept_err = select_intercept(python_like_select_rownames(summary(rep), 'beta')[,2]) ## std err
ggplot(cbind.data.frame(name=paste0('Transformation', 1:d), intercept_est, intercept_err),
       aes(x=name, y=intercept_est))+
  geom_errorbar(aes(ymin=intercept_est-intercept_err, ymax=intercept_est+intercept_err), width=.1)+
  geom_line() + geom_point()+theme_linedraw()+labs(x='Name of beta intercept', y='Estimate for beta value')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave("../results/LN_modelling/RELN_TMB_alphaTsagris_intercept_errorbar.png", width = 3, height = 2.5)

## slope

slope_est = select_slope_2(python_like_select_rownames(summary(rep), 'beta')[,1],v=F) ## estimate
slope_err = select_slope_2(python_like_select_rownames(summary(rep), 'beta')[,2],v=F) ## std err
ggplot(cbind.data.frame(name=paste0('Transformation', 1:d), slope_est, slope_err),
       aes(x=name, y=slope_est))+
  geom_errorbar(aes(ymin=slope_est-slope_err, ymax=slope_est+slope_err), width=.1)+
  geom_line() + geom_point()+theme_linedraw()+labs(x='Name of beta slope', y='Estimate for beta value')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave("../results/LN_modelling/RELN_TMB_alphaTsagris_slope_errorbar.png", width = 3, height = 2.5)


## Analysing the random effects
ggthemr::ggthemr_reset()
RE_est = python_like_select_rownames(summary(rep), 'u_large')[,1] ## estimate
RE_err = python_like_select_rownames(summary(rep), 'u_large')[,2] ## std err
df_RE = cbind.data.frame(name_indiv=factor(rep(1:num_indiv, d)),
                         name_LR = factor(paste0('Transformation', 1:(d))[rep(1:(d), each=num_indiv)]),
                         RE_est, RE_err)
# ggplot(df_RE,
#        # aes(x=interaction(name_LR,name_indiv), y=RE_est, col=(name_LR)))+
#        # aes(x=name_indiv, group = interaction(name_indiv, name_LR), y=RE_est, col=(name_LR)))+facet_wrap(.~name_indiv)+
#        aes( x = name_LR, y=RE_est, col=(name_LR)))+facet_wrap(.~name_indiv, nrow=2)+
#   geom_errorbar(aes(ymin=RE_est-RE_err, ymax=RE_est+RE_err), width=.1)+
#   geom_point()+labs(x='Name of beta slope', y='Estimate for beta value')
# ggsave("../results/LN_modelling/RELN_TMB_alphaTsagris_RE_coefs_errorbar.png", width = 20, height = 2.5)

ggplot(df_RE,
       # aes( x = factor(name_indiv, levels=rownames(exposures)[order(rowSums(abs(u_est)))]),
       aes( x = factor(name_indiv, levels=order(rowSums(abs(u_est)))),
            y=RE_est, group=name_LR))+facet_wrap(.~name_LR)+
  geom_errorbar(aes(ymin=RE_est-RE_err, ymax=RE_est+RE_err), width=.1)+
  geom_point()+geom_line()+labs(x='Ordered patients', y='Coefficient for random effects')
ggsave("../results/LN_modelling/RELN_TMB_alphaTsagris_RE_coefs_errorbar2.png", width = 20, height = 4.5)

ggplot(df_RE,
       # aes( x = factor(name_indiv, levels=rownames(exposures)[order(rowSums(abs(u_est)))]),
       aes( x = factor(name_indiv, levels=order(rowSums(abs(u_est)))),
            y=abs(RE_est), group=name_LR))+facet_wrap(.~name_LR)+
  geom_errorbar(aes(ymin=abs(RE_est)-RE_err, ymax=abs(RE_est)+RE_err), width=.1)+
  geom_point()+geom_line()+labs(x='Ordered patients', y='Coefficient for random effects')
ggsave("../results/LN_modelling/RELN_TMB_alphaTsagris_RE_coefs_errorbar2_abs.png", width = 20, height = 4.5)


ggplot(df_RE,
       # aes( x = factor(name_indiv, levels=rownames(exposures)[order(rowSums(abs(u_est)))]),
       aes( x = RE_est, col=name_LR))+facet_wrap(.~name_LR, nrow=1)+
  geom_density()+
  labs(x='Coefficients for random effects', y='Density')
ggsave("../results/LN_modelling/RELN_TMB_alphaTsagris_RE_coefs_densities.png", width = 10, height = 2.5)

ggplot(df_RE,
       # aes( x = factor(name_indiv, levels=rownames(exposures)[order(rowSums(abs(u_est)))]),
       aes( x = RE_est))+
  geom_density()+
  labs(x='Coefficients for random effects', y='Density')
ggsave("../results/LN_modelling/RELN_TMB_alphaTsagris_RE_coefs_density_all.png", width = 3, height = 3)

## suspicious that it seems to be the same gaussian that all random effects coefficients come from
sds = python_like_select_rownames(summary(rep), 'logs_sd_RE')[,1]
names(sds) = paste0('Transformation', 1:(d-1))
sds = melt(sapply(sds, rnorm, n = 2000, mean=0)); colnames(sds) = c('name_indiv', 'name_LR', 'RE_est')
ggplot(rbind.data.frame(cbind.data.frame(sds, RE_err=NA, sim='Sim'),
                        cbind.data.frame(df_RE, sim='Obs')), aes(x=RE_est, col=sim))+geom_density()+
  facet_wrap(.~name_LR , nrow=1)
ggsave("../results/LN_modelling/RELN_TMB_alphaTsagris_RE_coefs_density_all_with_sim.png", width = 8, height = 2)

# ## So, is it differentially abundant?
# ## differential abundance
wald_TMB_wrapper(rep, verbatim = FALSE)

generalisation_other_alpha = function(alpha_val){
  
  exposures_alphatrans = t(apply(sig_quants, 2, alphatrans, alpha=alpha_val))
  colnames(exposures_alphatrans)= paste0('Transformation', 1:ncol(exposures_alphatrans))
  
  TMB_data = list(Y = exposures_alphatrans,
                  num_individuals = num_indiv,
                  x = t(x),
                  z = sapply(unique(patient.meta$PATIENT_ID), function(i) as.numeric(patient.meta$PATIENT_ID == i)))
  
  TMB_params = list(beta = matrix(runif( 2*(d), min = -4, max = 4),
                                  nrow = 2, byrow=TRUE),
                    u_large = matrix(rep(1, (d)*num_indiv), nrow=num_indiv),
                    logs_sd_RE=rep(1, d),
                    logs_sd_FE=rep(1, d))
  
  obj <- MakeADFun(data = TMB_data, parameters = TMB_params, DLL="tmb_RE_lm", random = "u_large")
  # obj <- MakeADFun(data = TMB_data, parameters = TMB_params, DLL="tmb_RE" )
  opt <- do.call("optim", obj)
  opt
  opt$hessian ## <-- FD hessian from optim
  rep <- sdreport(obj)
  rep
  
  beta_est = matrix(python_like_select_name(rep$par.fixed, 'beta'), nrow = 2)
  u_est = matrix(python_like_select_name(rep$par.random, 'u_large'), ncol = (d))
  thetaLR = TMB_data$x %*% beta_est + TMB_data$z %*% u_est
  
  ## this below is not as easy as with the softmax! I need to find the inverse function
  inverse_alpha_trans(thetaLR[1,], alpha = alpha_val)
  # exposures_fitted = t(apply(thetaLR, 1, inverse_alpha_trans, alpha=alpha_val))
  # dim(exposures_fitted)
  # dim(exposures_alphatrans)
  # 
  # 
  # colnames(exposures_fitted) = colnames(exposures_alphatrans)
  
  df_scatter = cbind.data.frame(exposures=as.vector(exposures_alphatrans),
                                fitted_exposures= as.vector(thetaLR),
                                signature=rep(colnames(exposures_alphatrans),
                                              each=nrow(exposures_alphatrans)))
  
  df_matrix = cbind.data.frame(sample=rownames(exposures_alphatrans), exposures=(exposures_alphatrans),
                               fitted_exposures= as.vector(exposures_fitted))
  
    ggplot(df_scatter,
         aes(x=exposures, y=fitted_exposures, col=signature))+geom_point()+facet_wrap(.~signature, nrow=1)+theme_linedraw()+
    theme(axis.text.x = element_text(angle = 30, hjust = 1))+
    geom_abline(slope = 1, intercept = 0, col='black', )
  ggsave(paste0("../results/LN_modelling/alpha_trans/RELN_TMB_alphaTsagris", alpha_val, "_scatter.png"), width = 12, height = 2.5)
  alpha_val
  pdf(paste0("../results/LN_modelling/alpha_trans/RELN_TMB_alphaTsagris_", alpha_val, "matrices_exposures.pdf"), width = 6, height = 7)
  par(mfrow=c(1,2))
  image(t(exposures_alphatrans), main='True transformed exposures')
  image(t(thetaLR), main='Fitted transformed exposures')
  dev.off()
  
  ## Analysing the betas
  ## intercept
  ## to get the standard errors of estimates
  beta_est[1,]
  intercept_est = select_intercept(python_like_select_rownames(summary(rep), 'beta')[,1]) ## estimate
  intercept_err = select_intercept(python_like_select_rownames(summary(rep), 'beta')[,2]) ## std err
  ggplot(cbind.data.frame(name=paste0('Transformation', 1:d), intercept_est, intercept_err),
         aes(x=name, y=intercept_est))+
    geom_errorbar(aes(ymin=intercept_est-intercept_err, ymax=intercept_est+intercept_err), width=.1)+
    geom_line() + geom_point()+theme_linedraw()+labs(x='Name of beta intercept', y='Estimate for beta value')+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  ggsave(paste0("../results/LN_modelling/alpha_trans/RELN_TMB_alphaTsagris", alpha_val, "_intercept_errorbar.png"), width = 3, height = 2.5)
  
  ## slope
  
  slope_est = select_slope_2(python_like_select_rownames(summary(rep), 'beta')[,1],v=F) ## estimate
  slope_err = select_slope_2(python_like_select_rownames(summary(rep), 'beta')[,2],v=F) ## std err
  ggplot(cbind.data.frame(name=paste0('Transformation', 1:d), slope_est, slope_err),
         aes(x=name, y=slope_est))+
    geom_errorbar(aes(ymin=slope_est-slope_err, ymax=slope_est+slope_err), width=.1)+
    geom_line() + geom_point()+theme_linedraw()+labs(x='Name of beta slope', y='Estimate for beta value')+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  ggsave(paste0("../results/LN_modelling/alpha_trans/RELN_TMB_alphaTsagris", alpha_val, "_slope_errorbar.png"), width = 3, height = 2.5)
  
  
  ## Analysing the random effects
  ggthemr::ggthemr_reset()
  RE_est = python_like_select_rownames(summary(rep), 'u_large')[,1] ## estimate
  RE_err = python_like_select_rownames(summary(rep), 'u_large')[,2] ## std err
  df_RE = cbind.data.frame(name_indiv=factor(rep(1:num_indiv, d)),
                           name_LR = factor(paste0('Transformation', 1:(d))[rep(1:(d), each=num_indiv)]),
                           RE_est, RE_err)
  # ggplot(df_RE,
  #        # aes(x=interaction(name_LR,name_indiv), y=RE_est, col=(name_LR)))+
  #        # aes(x=name_indiv, group = interaction(name_indiv, name_LR), y=RE_est, col=(name_LR)))+facet_wrap(.~name_indiv)+
  #        aes( x = name_LR, y=RE_est, col=(name_LR)))+facet_wrap(.~name_indiv, nrow=2)+
  #   geom_errorbar(aes(ymin=RE_est-RE_err, ymax=RE_est+RE_err), width=.1)+
  #   geom_point()+labs(x='Name of beta slope', y='Estimate for beta value')
  # ggsave("../results/LN_modelling/RELN_TMB_alphaTsagris_RE_coefs_errorbar.png", width = 20, height = 2.5)
  
  ggplot(df_RE,
         # aes( x = factor(name_indiv, levels=rownames(exposures)[order(rowSums(abs(u_est)))]),
         aes( x = factor(name_indiv, levels=order(rowSums(abs(u_est)))),
              y=RE_est, group=name_LR))+facet_wrap(.~name_LR)+
    geom_errorbar(aes(ymin=RE_est-RE_err, ymax=RE_est+RE_err), width=.1)+
    geom_point()+geom_line()+labs(x='Ordered patients', y='Coefficient for random effects')
  ggsave(paste0("../results/LN_modelling/alpha_trans/RELN_TMB_alphaTsagris", alpha_val, "_RE_coefs_errorbar2.png"), width = 20, height = 4.5)
  
  ggplot(df_RE,
         # aes( x = factor(name_indiv, levels=rownames(exposures)[order(rowSums(abs(u_est)))]),
         aes( x = factor(name_indiv, levels=order(rowSums(abs(u_est)))),
              y=abs(RE_est), group=name_LR))+facet_wrap(.~name_LR)+
    geom_errorbar(aes(ymin=abs(RE_est)-RE_err, ymax=abs(RE_est)+RE_err), width=.1)+
    geom_point()+geom_line()+labs(x='Ordered patients', y='Coefficient for random effects')
  ggsave(paste0("../results/LN_modelling/alpha_trans/RELN_TMB_alphaTsagris", alpha_val, "_RE_coefs_errorbar2_abs.png"), width = 20, height = 4.5)
  
  
  ggplot(df_RE,
         # aes( x = factor(name_indiv, levels=rownames(exposures)[order(rowSums(abs(u_est)))]),
         aes( x = RE_est, col=name_LR))+facet_wrap(.~name_LR, nrow=1)+
    geom_density()+
    labs(x='Coefficients for random effects', y='Density')
  ggsave(paste0("../results/LN_modelling/alpha_trans/RELN_TMB_alphaTsagris", alpha_val, "_RE_coefs_densities.png"), width = 10, height = 2.5)
  
  ggplot(df_RE,
         # aes( x = factor(name_indiv, levels=rownames(exposures)[order(rowSums(abs(u_est)))]),
         aes( x = RE_est))+
    geom_density()+
    labs(x='Coefficients for random effects', y='Density')
  ggsave(paste0("../results/LN_modelling/alpha_trans/RELN_TMB_alphaTsagris", alpha_val, "_RE_coefs_density_all.png"), width = 3, height = 3)
  
  ## suspicious that it seems to be the same gaussian that all random effects coefficients come from
  sds = python_like_select_rownames(summary(rep), 'logs_sd_RE')[,1]
  names(sds) = paste0('Transformation', 1:(d-1))
  sds = melt(sapply(sds, rnorm, n = 2000, mean=0)); colnames(sds) = c('name_indiv', 'name_LR', 'RE_est')
  ggplot(rbind.data.frame(cbind.data.frame(sds, RE_err=NA, sim='Sim'),
                          cbind.data.frame(df_RE, sim='Obs')), aes(x=RE_est, col=sim))+geom_density()+
    facet_wrap(.~name_LR , nrow=1)
  ggsave(paste0("../results/LN_modelling/alpha_trans/RELN_TMB_alphaTsagris", alpha_val, "_RE_coefs_density_all_with_sim.png"), width = 8, height = 2)
  
  # ## So, is it differentially abundant?
  # ## differential abundance
  # wald_TMB_wrapper(rep, verbatim = FALSE)
  return(rep)
  
}

various_alpha_results = lapply(c(0.8, 1, 1.2, 1.5, 2), generalisation_other_alpha)
various_alpha_results = lapply(c(0.1, 0.2, 0.5), generalisation_other_alpha)

## correlation of beta values
## extremely similar betas regardless of alpha value
pairs(python_like_select_rownames(sapply(various_alpha_results, '[[', "par.fixed"), "beta"))

## statistically significant in all
sapply(various_alpha_results, wald_TMB_wrapper, verbatim = FALSE)

##----------------------------------------------------------------------------------------------------##
## removing s5

d <- 5
exposures_alphatrans = t(apply(t(normalise_cl(sig_quants[-5,])), 2, alphatrans, alpha=alpha_val))
dim(exposures_alphatrans)
x

TMB_data = list(Y = exposures_alphatrans,
                num_individuals = num_indiv,
                x = t(x),
                z = sapply(unique(patient.meta$PATIENT_ID), function(i) as.numeric(patient.meta$PATIENT_ID == i)))

sapply(TMB_data, dim)
num_indiv
dim(exposures_alphatrans)

TMB_params = list(beta = matrix(runif( 2*(d), min = -4, max = 4),
                                nrow = 2, byrow=TRUE),
                  u_large = matrix(rep(1, (d)*num_indiv), nrow=num_indiv),
                  logs_sd_RE=rep(1, d),
                  logs_sd_FE=rep(1, d))

obj <- MakeADFun(data = TMB_data, parameters = TMB_params, DLL="tmb_RE_lm", random = "u_large")
# obj <- MakeADFun(data = TMB_data, parameters = TMB_params, DLL="tmb_RE" )
opt <- do.call("optim", obj)
opt
opt$hessian ## <-- FD hessian from optim
rep <- sdreport(obj)
rep

plot_betas(rep)

wald_TMB_wrapper(rep) ## statistically signif

beta_est = matrix(python_like_select_name(rep$par.fixed, 'beta'), nrow = 2)
u_est = matrix(python_like_select_name(rep$par.random, 'u_large'), ncol = (d))
thetaLR = TMB_data$x %*% beta_est + TMB_data$z %*% u_est

colnames(exposures_alphatrans) <- paste0('Transformation ', 1:ncol(exposures_alphatrans))
df_scatter = cbind.data.frame(exposures=as.vector(exposures_alphatrans),
                              fitted_exposures= as.vector(thetaLR),
                              signature=rep(colnames(exposures_alphatrans),
                                            each=nrow(exposures_alphatrans)))

df_matrix = cbind.data.frame(sample=rownames(exposures_alphatrans), exposures=(exposures_alphatrans),
                             fitted_exposures= as.vector(exposures_fitted))


ggthemr::ggthemr("fresh")
ggplot(df_scatter,
       aes(x=exposures, y=fitted_exposures, col=signature))+geom_point()+facet_wrap(.~signature, nrow=1)+theme_linedraw()+
  theme(axis.text.x = element_text(angle = 30, hjust = 1))+
  geom_abline(slope = 1, intercept = 0, col='black')
ggsave("../../../results/_previous_results/LN_modelling/alpha_trans/RELN_TMB_alphaTsagris_scatter_nos5.png", width = 12, height = 2.5)
