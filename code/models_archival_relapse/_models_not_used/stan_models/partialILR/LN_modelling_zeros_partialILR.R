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

model_file_name = "../../../ProjectOvarianMultisampleTree/code/files_analysis/modelling_LN/stan_fit_simple_LN.stan"
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

#-------------------------------------------------------------------------------------------#

## Transform data with partial ILR
which_zero = t(apply(exposures, 1, function(i) as.numeric((i==0)) ))

library(compositions)

irl_base_complete = compositions::ilrBase(D = 7)

give_partial_ilr_basis = function(which_zero_vector){
  irl_base_partial1 = irl_base_complete
  for(idx_missing in which_zero_vector){
    irl_base_partial1[,idx_missing-1] = 0
    irl_base_partial1[idx_missing,] = 0
    if(idx_missing <= ncol(irl_base_partial1) ){
      for(i in idx_missing:ncol(irl_base_partial1)){
        above_vals = irl_base_partial1[(1:(i-1))[-which_zero_vector],i]
        ## substitute above vals
        irl_base_partial1[(1:(i-1))[-which_zero_vector],i] = -irl_base_partial1[(i+1),i]/(i-length(above_vals))
      }
    }
  }
  return(irl_base_partial1)
}

exposures[1,]
which_zero[1,]
as.vector(compositions::ilr(exposures)[1,])
give_partial_ilr_basis(which(which_zero[1,] == 1))

exposures_partial_irl = t(apply(exposures, 1, function(rw){
  .which_zero = as.vector(which(rw==0))
  .basis_row = give_partial_ilr_basis(.which_zero)
  compositions::clr(rw) %*% .basis_row
}))

## note that the partial ILR that I have  computed is different from the ILR setting zero entries as zero ILRs
## different
compositions::clr(exposures[1,]) %*% give_partial_ilr_basis(which_zero[[1]])
compositions::clr(exposures[1,]) %*% irl_base_complete

# To show that I am computing the ILR correctly given a basis
# plot(as.vector(compositions::ilr(exposures)[1,]),
#      compositions::clr(exposures[1,]) %*% irl_base_complete)
# abline(coef = c(0,1))

#-------------------------------------------------------------------------------------------#
which_zero
exposures

which_zero_ilr = which_zero[,-1]
subset_samples = rowSums(which_zero_ilr) < (ncol(exposures)-2)
exposures_partial_irl_partial = exposures_partial_irl[subset_samples,]

stan_data = list(n=nrow(exposures_partial_irl_partial),
                 d = ncol(exposures),
                 zeros_matrix = which_zero_ilr[subset_samples,],
                 W = exposures_partial_irl_partial)

params = c('mu', 'Sigma')


model_file_name = "../../../ProjectOvarianMultisampleTree/code/files_analysis/modelling_LN/stan_fit_simple_zeros_partialILR.stan"
rstan_options(auto_write = TRUE)
stanc(model_file_name)

fit_stan <- stan(file = model_file_name, data = stan_data,
                 iter = Nits, chains = 4, cores = 4, thin = 1, pars = params,
                 control = list(stepsize=3, adapt_delta=.95, max_treedepth=15))

## Good convergence
max(bayesplot::rhat(fit_stan))
posterior2 = as.matrix(fit_stan)

# saveRDS(posterior2, file = paste0("../../../out/robj/inference/simple_", sample_name, "_posteriors.RDS"))

names_betas = colnames(posterior2)[grep("beta", colnames(posterior2))]
p <- bayesplot::mcmc_trace(posterior2,  pars = names_betas,
                           facet_args = list(nrow = 1, labeller = label_parsed))+ theme(text = element_text(size=10))
p + facet_text(size = 6)

p <- bayesplot::mcmc_trace(posterior2,  pars = python_like_select(colnames(posterior2), 'mu'),
                           facet_args = list(nrow = 1, labeller = label_parsed))+ theme(text = element_text(size=10))
p + facet_text(size = 6)


pairs(fit_stan, pars = python_like_select(colnames(posterior2), "mu"))

bayesplot::mcmc_areas(posterior2, pars = colnames(posterior2)[-ncol(posterior2)])+ggtitle('Slope')
bayesplot::mcmc_parcoord(posterior2)+ggtitle('Slope')

traceplot(fit_stan, pars = c("mu", "Sigma"), inc_warmup = TRUE, nrow = 2)

size_sim = 400
stopifnot(size_sim %% 2 == 0)
idx_posteriors = sample(1:nrow(posterior2), size = size_sim, replace = FALSE)

exposures_NAzeros = exposures
exposures_NAzeros[exposures_NAzeros==0] = NA
simulate_from_model = function(idx_posterior){
  posterior_row = posterior2[idx_posterior,]
  mu = python_like_select_name(posterior_row, 'mu')
  multiv = t(apply(t(matrix(mu)), 1, mvtnorm::rmvnorm, n = 1,
                   sigma =  matrix(python_like_select_name(posterior_row, 'Sigma'),
                                   ncol = ncol(exposures_NAzeros)-1)))
  return(softmax(c(multiv, 0)))
}

sim_results = t(sapply(idx_posteriors, simulate_from_model))
colnames(sim_results) = colnames(exposures_NAzeros)
splits_group = rep(c(1,2), c(nrow(sim_results), nrow(exposures_NAzeros)))
pairs(rbind(sim_results, exposures_NAzeros), col=c('#080a4d', 'red')[splits_group], pch=19, cex=0.1)
## simulation in dark blue; observed exposures in red
## here I am only plotting non-zero exposures

ggplot(melt(list(simulation=sim_results, observed=exposures_NAzeros)), aes(x=value, group=L1, col=L1))+facet_wrap(.~Var2)+geom_density()

