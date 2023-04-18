## Infence using stan_fit_LNM for data in the simplex
## Simple LN, no fixed or random effects

rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

set.seed(1234)

library(rstan)
library(uuid)
library(ggplot2)
library(optparse)
library(reshape2)
library(pheatmap)
library(bayesplot)
source("../../../GlobalDA/code/2_inference/helper/helper_DA_stan.R")
source("../../../GlobalDA/code/3_analysis/helper/helper_analyse_posteriors.R")
source("../../../GlobalDA/code/2_inference_TMB/helper_TMB.R")
Nits = 2500

model_file_name = "../../../ProjectOvarianMultisampleTree/code/files_analysis/modelling_LN/stan_fit_LN.stan"
rstan_options(auto_write = TRUE)
stanc(model_file_name)

#-------------------------------------------------------------------------------------------#
give_short_names = function(i){
  if(grepl('TCGA', i)){
    paste0(strsplit(i, split = '-')[[1]][1:3], collapse='-')
  }else{
    i
  }
}

load("../data/britroc_30kb_signature_data.rds")
# patient.meta
# sample_by_component
# sig_quants

exposures = t(sig_quants)
exposures = normalise_rw(exposures + 1e-4)
patient.meta = patient.meta[match(rownames(exposures), as.character(patient.meta$SAMPLE_ID)),]
sample_by_component = sample_by_component[match(rownames(exposures), rownames(sample_by_component)),]

#-------------------------------------------------------------------------------------------#

x = t(cbind(1, as.numeric(factor(patient.meta$group))-1))

stan_data = list(n=nrow(exposures),
                 d = ncol(exposures),
                 W = exposures,
                 x = x)

params = c('Sigma', 'beta')

fit_stan <- stan(file = model_file_name, data = stan_data,
                 iter = Nits, chains = 4, cores = 4, thin = 1, pars = params,
                 control = list(stepsize=3, adapt_delta=.95, max_treedepth=15))

uuid_out = uuid::UUIDgenerate()
saveRDS(fit_stan, paste0("../out/inference/", "LN_FE", uuid_out, ".RDS"))

fit_stan <- readRDS(paste0("../out/inference/LN_FEe167991d-9fc4-4b2b-889a-3906eac198c4.RDS"))


# pairs(fit_stan)
# 
# ## Good convergence
# max(bayesplot::rhat(fit_stan))
posterior2 = as.matrix(fit_stan)

names_betas = colnames(posterior2)[grep("beta", colnames(posterior2))]
p <- bayesplot::mcmc_trace(posterior2,  pars = names_betas,
                           facet_args = list(nrow = 1, labeller = label_parsed))+ theme(text = element_text(size=10))
p + facet_text(size = 6)


## looking at the correlation between inferred parameters (in this case, just betas)
## there is a clear correlations between the intersect and the corresponding slope, for each of the (sub)features index by 1, ..., d-1
png(paste0("../results/LN_modelling/FELN_pairs_plot.png"), width = 15, height = 15, units = "in", res = 300)
pairs(fit_stan, pars = python_like_select(colnames(posterior2), "beta"))
dev.off()
# # ## beta coefficients for intersect
# # pairs(fit_stan, pars = names(fit_stan)[grep('beta\\[1,', names(fit_stan))], text.panel = "Coefficients for the intercept")
# # ## beta coefficients for slope
# # pairs(fit_stan, pars = names(fit_stan)[grep('beta\\[2,', names(fit_stan))], text.panel = "Coefficients for the slope")
# 
# ## Reminder that beta[1,x] corresponds to the intercept, and beta[2,x] to the slopes
# ## plot all but LP
pdf(paste0("../results/LN_modelling/FELN_pars.pdf"))
bayesplot::mcmc_areas(posterior2, pars = colnames(posterior2)[-ncol(posterior2)])+ggtitle('Slope')
dev.off()

pdf(paste0("../results/LN_modelling/FELN_parcoord.pdf"), width = 12)
bayesplot::mcmc_parcoord(posterior2)+ggtitle('Slope')
dev.off()
 
png(paste0("../results/LN_modelling/FELN_traceplot.png"), width = 14, height = 4, units = "in", res = 300)
traceplot(fit_stan, pars = c("beta", "Sigma"), inc_warmup = TRUE, nrow = 2)
dev.off()

## let's simulate under the model with the inferred parameters and compare it to the actual data
size_sim = 400
stopifnot(size_sim %% 2 == 0)
idx_posteriors = sample(1:nrow(posterior2), size = size_sim, replace = FALSE)

simulate_from_model = function(idx_posterior){
  posterior_row = posterior2[idx_posterior,]
  beta = matrix(python_like_select_name(posterior_row, 'beta'), nrow=2)
  multiv = t(apply(t(x) %*% beta, 1, mvtnorm::rmvnorm, n = 1,
                   sigma =  matrix(python_like_select_name(posterior_row, 'Sigma'),
                                   ncol = ncol(exposures)-1)))
  return(softmax(cbind(multiv, 0)))
}

sim_results = do.call('rbind', lapply(idx_posteriors, simulate_from_model))
colnames(sim_results) = colnames(exposures)
splits_group = rep(c(1,2), c(nrow(sim_results), nrow(exposures)))

png(paste0("../results/LN_modelling/FELN_comparison_with_simulation.png"), width = 10, height = 10, units = "in", res = 300)
pairs(rbind(sim_results, exposures), col=c('#080a4d', 'red')[splits_group], pch=19, cex=0.1)
dev.off()

## do this for contour plot


# sim_results_split = split(sim_results, factor(splits_group))
# min(sim_results_split[[1]])
# min(sim_results_split[[2]])
# 
