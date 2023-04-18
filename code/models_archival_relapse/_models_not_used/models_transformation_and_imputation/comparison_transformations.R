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
library(TMB)
library(pheatmap)
library(bayesplot)
library(dplyr)
library(parallel)
source("../../../GlobalDA/code/2_inference/helper/helper_DA_stan.R")
source("../../../GlobalDA/code/3_analysis/helper/helper_analyse_posteriors.R")
source("../../../GlobalDA/code/2_inference_TMB/helper_TMB.R")
source("functions.R")
Nits = 2500

#-------------------------------------------------------------------------------------------#
TMB::compile("tmb_RE/tmb_RE.cpp", "-std=gnu++17")
dyn.load(dynlib("tmb_RE/tmb_RE"))
#-------------------------------------------------------------------------------------------#

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

## compute distances
nonredundant_pairs = do.call('rbind', sapply(1:(nrow(exposures)-1), function(i) cbind(i, (i+1):nrow(exposures))))
distances_exposures = sapply(1:nrow(nonredundant_pairs),
                             function(i) dist(rbind(exposures[nonredundant_pairs[i,1],],
                                                    exposures[nonredundant_pairs[i,2],])))
zeros_imputation_vec = c(1e-1, 4e-2, 1e-2, 1e-3, 1e-4)
transformations_list = mclapply(zeros_imputation_vec, function(imputation){
  as(compositions::alr(normalise_rw(exposures + imputation)), 'matrix')
})
alphas_vec = c(0.2, 0.8, 1, 1.1, 1.2, 1.4)
alpha_transformations_list = mclapply(alphas_vec, function(alpha_val){
  t(apply(exposures, 1, alphatrans, alpha=alpha_val))
})
transformations_list = c(transformations_list, alpha_transformations_list)
distances_exposures_transformation = mclapply(transformations_list, FUN = function(transformed_dataset){
  sapply(1:nrow(nonredundant_pairs), function(i){dist(rbind(transformed_dataset[nonredundant_pairs[i,1],],
                                                            transformed_dataset[nonredundant_pairs[i,2],]))})
})

# par(mfrow=c(1,length(transformations_list)))
# sapply(1:length(transformations_list), function(j)
#   plot(distances_exposures, distances_exposures_transformation[[j]])
# )

cors = sapply(distances_exposures_transformation, function(j) cor(distances_exposures, j))
names(cors) = c(paste0('Zero inputation: ', zeros_imputation_vec), paste0('Alpha: ', alphas_vec))
cors
xtable::xtable(t(t(cors)))

## compute order of signatures and how this changes (alternative to the distance)
sapply(transformations_list)

cors_rank = sapply(transformations_list, function(j){
  cor(apply(apply(exposures, 2, order), 1, mean),
      apply(apply(j, 2, order), 1, mean))
})

plot(cors, cors_rank)

ggplot((cbind.data.frame(cors, cors_rank, type=sapply(names(cors), function(i) strsplit(i, ':')[[1]][1]),
                         val=sapply(names(cors), function(i) strsplit(i, ':')[[1]][2]))),
       aes(x=cors, y=cors_rank, col=as.numeric(val), shape=type))+geom_point()
ggsave("../results/LN_modelling/transformations_rank_correlation.pdf", width=4, height = 3)
# -------------------------

## Marginal over scaling factors
# l = 1
alpha_vec = seq(0.5, 1.5, length.out = 20)=
# i = nonredundant_pairs[l,1]
# j = nonredundant_pairs[l,2]
vec_k = seq(0.8, 4, length.out = 20)
k_it = 1
k = vec_k[k_it]

give_sum_dist_per_k = function(k){
  ## should we scale any two samples independently? probably...!!
  scaled_exposures = exposures*k
  marginal_distance_scaled = apply(nonredundant_pairs, 1, function(pair_vec){
    scaled_exposures_i = scaled_exposures[pair_vec[1],]
    scaled_exposures_j = scaled_exposures[pair_vec[2],]
    as.numeric(dist(rbind(scaled_exposures_i, scaled_exposures_j)))
  } )
  sum_marginal_distance_scaled = sum(marginal_distance_scaled)
  # marginal_distance_scaled[l]/sum_marginal_distance_scaled
  
  alpha_vec_it = alpha_vec[1]
  ## compute all the distances of the scaled composition (shared by all alphas)
  transformed_exposures =   t(apply(exposures, 1, alphatrans, alpha=alpha_vec_it))
  marginal_distance_transformed = apply(nonredundant_pairs, 1, function(pair_vec){
    scaled_exposures_i = transformed_exposures[pair_vec[1],]
    scaled_exposures_j = transformed_exposures[pair_vec[2],]
    as.numeric(dist(rbind(scaled_exposures_i, scaled_exposures_j)))
  } )
  sum_marginal_distance_transformed = sum(marginal_distance_transformed)
  
  # transformed_exposures[nonredundant_pairs[l,1],]
  # transformed_exposures[nonredundant_pairs[l,2],]
  
  sum_k = sapply(1:length(nonredundant_pairs), function(l) marginal_distance_scaled[l]/sum_marginal_distance_scaled - marginal_distance_transformed[l]/sum_marginal_distance_transformed)
  return(sum_k)
}

several_k = mclapply(vec_k, give_sum_dist_per_k)


## compute all the distacnes of the transformation

