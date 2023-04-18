#-----------------------------------------------------------------------------------------------#
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

segs = read.table("../../data/britroc_30kb_ds_absCopyNumber_segmentTable.tsv", header = T)
library(reshape2)
library(ggplot2)
library(dplyr)
library(gridExtra)
source("functions.R")
source("../functions.R")
source("../../../britroc-cnsignatures-bfb69cd72c50/main_functions.R")
source("../../../britroc-cnsignatures-bfb69cd72c50/helper_functions.R")
source("../../../../GlobalDA/code/2_inference_TMB/helper_TMB.R")
components = readRDS("../../../britroc-cnsignatures-bfb69cd72c50/data/component_parameters.rds")
feat_sig_mat = readRDS("../../../britroc-cnsignatures-bfb69cd72c50/data/feat_sig_mat.rds")
phils_exposures = t(sig_quants) ## this are phil's exposures
load("../../data/britroc_30kb_signature_data.rds") ## patient meta

# segs = split.data.frame(segs, segs$sample)
# cnfeatures = extractCopynumberFeatures(segs)

# This is not needed. This is only for signature de novo generation
# sumofpost = lapply(1:length(components), function(i){
#   calculateSumOfPosteriors(cnfeatures[[i]], components = components[[i]], name = names(components)[i])
# })
# saveRDS(sumofpost, "../../out/signature_extraction/sumofpost.RDS")

# SxCmatrix = generateSampleByComponentMatrix(cnfeatures)
# saveRDS(SxCmatrix, "../../out/signature_extraction/SxCmatrix.RDS")
SxCmatrix = readRDS("../../out/signature_extraction/SxCmatrix.RDS")


sigs = quantifySignatures_alt(SxCmatrix)
patient.meta = patient.meta[match(colnames(sigs), as.character(patient.meta$SAMPLE_ID)),]

# sigs_phil = quantifySignatures(SxCmatrix)
# phils_exposures
# all(as.vector(t(sigs_phil)) == as.vector(phils_exposures))
# sigs

## compare to the signatures I was given
# load(("../../data/britroc_30kb_signature_data.rds"))
# sigs_from_phil = sig_quants
# saveRDS(sigs_from_phil, "../../data/britroc_30kb_signatures.rds")
# 
# if(!all(match(colnames(sigs_from_phil), colnames(sigs)) == 1:ncol(sigs))){
#   stop('Wrong order')
# }
# 
# plot(as.vector(sigs_from_phil), as.vector(sigs))
# 
# sum(sigs_from_phil == 0)
# sum(sigs == 0)


sigs_from_phil = readRDS("../../data/britroc_30kb_signatures.rds")
#-----------------------------------------------------------------------------------------------#


#-----------------------------------------------------------------------------------------------#
## for each signature and for each sample, find how robust their - potentially zero - exposure is

## super naive first approach
# .SxCmatrix_mod = SxCmatrix+1
# .sigs_mod = quantifySignatures_alt(.SxCmatrix_mod)
# 
# ggplot(melt(SxCmatrix), aes(x=value))+geom_density()+facet_wrap(.~Var2, scales = "free")
# 
# replicate_add_noise = function(){
#   .SxCmatrix_mod = SxCmatrix+1
#   .sigs_mod = quantifySignatures_alt(.SxCmatrix_mod)
#   return(.sigs_mod)
# }

# noisy_replicates = replicate(100, replicate_add_noise(), simplify = F)
# melt_noisy_replicates = melt(noisy_replicates)
# melt_noisy_replicates$sig_sample = paste0(melt_noisy_replicates$Var1, '_', melt_noisy_replicates$Var2)
# average_exposures = melt_noisy_replicates %>% group_by(sig_sample) %>% summarize(mean=mean(value)) %>% arrange(mean)
# melt_noisy_replicates$sig_sample = factor(melt_noisy_replicates$sig_sample, levels=average_exposures %>% select(sig_sample) %>% unlist())
## sort by the highest and lowest

# ggplot(melt_noisy_replicates[melt_noisy_replicates$sig_sample %in% (average_exposures %>% filter(mean<0.0001) %>% select(sig_sample) %>% unlist),],
#        aes(x=sig_sample, y=value))+geom_point()


plot(SxCmatrix[1,1:10])

# ## the higher the entropy, the higher the perturbation
# perturbation = c(0,0,0,0,0,0,0,0,0,1)
# perturbation = rep(1,10)
# entropy::entropy(normalise_rw(rep(1,10)))
# entropy::entropy(normalise_rw(c(0,0,0,0,0,0,0,0,0,1)))
# install.packages('entropy')
# library(entropy)
# entropy::entropy(perturbation)
# sum(SxCmatrix[1,1:10]) * normalise_rw(normalise_rw(SxCmatrix[1,1:10]) * perturbation)
# normalise_cl

## select where to move it from
## do this procedure for some number of steps

idx_sample = 1
grouping_features = gsub('[[:digit:]]+', '', colnames(SxCmatrix))
which_features = lapply(unique(grouping_features), function(i) which(grouping_features == i))
names(which_features) = unique(grouping_features)

transition_matrices = lapply(unique(grouping_features), function(i){
  len = sum(grouping_features == i)
  transition_matrix = matrix(0, len, len)
  for(i in 1:len){
    transition_matrix[i,i-1] = 1
    if(i<len)  transition_matrix[i,i+1] = 1
  }
  transition_matrix = normalise_rw(transition_matrix)
  transition_matrix
}); names(transition_matrices) = unique(grouping_features)

mod_version = function(idx_sample, fraction_to_transport=0.02, SxCmatrix){
  mod_sop = sapply(unique(grouping_features), function(i){
  feature_vector = SxCmatrix[idx_sample,which_features[[i]]]
  
  count = 0; max_moves = 20
  ## MC for changing the proportions
  while(count < max_moves){
    ## select which proportion to move
    idx = sample(1:length(feature_vector), size = 1)
    .amount = rnorm(1, mean = sum(feature_vector)*fraction_to_transport, sd = fraction_to_transport*50)
    if(feature_vector[idx] > .amount){
      idx_to = sample(1:length(feature_vector), size = 1, prob = transition_matrices[[i]][idx,])
      feature_vector[idx_to] = feature_vector[idx_to] + .amount
      feature_vector[idx] = feature_vector[idx] - .amount
      count = count + 1
    }
  }
  return(feature_vector)
})
  names(mod_sop) = NULL
  return(unlist(mod_sop))
}

perturbed_sc = lapply(1:ncol(sigs), function(sampl) replicate(n = 20, expr = mod_version(idx_sample = sampl, SxCmatrix = SxCmatrix), simplify = F))
perturbed_sc = lapply(perturbed_sc, function(i) do.call('rbind', i))
perturbed_sc_rbind = do.call('rbind', perturbed_sc)


## compare the distributions of the simulations to the distributions of the original sums of posteriors
png("../../results/zeros_modelling/histograms_components.png", width = 10, height = 7, units = "in", res = 300)
grid.arrange(ggplot(melt(sapply(which_features, function(i) SxCmatrix[,i])), aes(x=Var2, y=value))+geom_bar(stat = "identity")+
               facet_wrap(.~L1, nrow=1, scales = "free_x")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle('Real SxC distributions'),
ggplot(melt(sapply(which_features, function(i) perturbed_sc_rbind[,i])), aes(x=Var2, y=value))+geom_bar(stat = "identity")+
  facet_wrap(.~L1, nrow=1, scales = "free_x")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle('Simulated SxC distributions')
)
dev.off()

perturbed_sc_sigs = lapply(perturbed_sc, function(i) t(quantifySignatures_alt(i)))
names(perturbed_sc_sigs) = colnames(sigs)

## See if the simulated exposures are around the value from Phil
phils_exposures
perturbed_sc_sigs

sim_in_interval = sapply(1:nrow(phils_exposures), function(sample_idx){
  sapply(1:7, function(j){
    .sim_vals = perturbed_sc_sigs[[sample_idx]][,j]
    ( phils_exposures[sample_idx,j] >= min(.sim_vals) ) & ( phils_exposures[sample_idx,j] <= max(.sim_vals) ) 
    })
})
sim_in_95_interval = sapply(1:nrow(phils_exposures), function(sample_idx){
  sapply(1:7, function(j){
    .sim_vals = perturbed_sc_sigs[[sample_idx]][,j]
    .quant = quantile(.sim_vals, c(0.05, 0.95))
    ( phils_exposures[sample_idx,j] >= .quant[1] ) & ( phils_exposures[sample_idx,j] <= .quant[2] ) 
  })
})
sum(sim_in_interval)/length(sim_in_interval)
sum(sim_in_95_interval)/length(sim_in_95_interval)
image(sim_in_interval)

## example of one sample
ggplot(melt(perturbed_sc_sigs[[2]]), aes(x=Var2, y=value))+geom_point()


zeros_df2_0005 = give_df_zeros_given_thresh(threshold = 0.005, perturbed_sc_sigs = perturbed_sc_sigs, sigs_from_phil = sigs_from_phil)

perturbed_sc_sigs_melt = melt(perturbed_sc_sigs)
perturbed_sc_sigs_melt$interaction = paste0(perturbed_sc_sigs_melt$L1, '-', perturbed_sc_sigs_melt$Var2)
perturbed_sc_sigs_melt$interaction = factor(perturbed_sc_sigs_melt$interaction, levels=unique(perturbed_sc_sigs_melt$interaction[order(perturbed_sc_sigs_melt$value)]))
perturbed_sc_sigs_melt$agreement = (zeros_df2_0005$transition %in% c("Consistent nonzero", "Consistent zero"))[match(perturbed_sc_sigs_melt$interaction, paste0(zeros_df2_0005$sample, '-', zeros_df2_0005$sig))]

ggplot(perturbed_sc_sigs_melt, aes(x=interaction, y=value, col=agreement))+geom_point(size=0.2)+facet_wrap(.~Var2)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave("../../results/zeros_modelling/robustness_zero_from_simulation_0005_ordered_scatterplot.pdf", height = 7, width = 8)

rename_columns = function(i, replacement){
  colnames(i) = replacement
  i
}

sigs_from_phil_melt = rename_columns(melt(sigs_from_phil), c('Var2', 'L1', 'value'))
sigs_from_phil_melt$interaction = paste0(sigs_from_phil_melt$L1, '-', sigs_from_phil_melt$Var2)
sigs_from_phil_melt$interaction = factor(sigs_from_phil_melt$interaction, levels=unique(sigs_from_phil_melt$interaction[order(sigs_from_phil_melt$value)]))

perturbed_sc_sigs_melt$interaction = factor(as.character(perturbed_sc_sigs_melt$interaction),
                                               levels=levels(sigs_from_phil_melt$interaction))
ggplot()+
  geom_point(data = perturbed_sc_sigs_melt, aes(x=interaction, y=value, col=agreement), size=0.2)+
  geom_line(data=sigs_from_phil_melt, col='red', aes(x=interaction, y=value, group=Var2))+
  facet_wrap(.~Var2)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave("../../results/zeros_modelling/robustness_zero_from_simulation_0005_ordered_scatterplot2.pdf", height = 7, width = 8)

exposures_0005 = sigs
for(zero_mine in 1:nrow(zeros_df2_0005)){
  if(zeros_df2_0005[zero_mine,'my_zeros']){
    exposures_0005[as.character(zeros_df2_0005[zero_mine,'sig']), as.character(zeros_df2_0005[zero_mine,'sample'])] = 0
  }
}
exposures_0005 = t(normalise_cl(exposures_0005))

plot(as.vector(exposures_0005), as.vector(sigs), col=factor(as.vector(exposures_0005 == sigs)))


pairs(cbind(as.vector(sigs), as.vector(phils_exposures), as.vector(exposures_0005)))

obj_original = list(Y = t(sigs),
                num_individuals = length(unique(patient.meta$PATIENT_ID)),
                d = 7,
                n = ncol(sigs),
                x = cbind(1, as.numeric(as.factor(patient.meta$group))-1),
                z = give_z_matrix_from_labels(lbls = patient.meta$PATIENT_ID))
obj_0005 = list(Y = t(exposures_0005),
           num_individuals = length(unique(patient.meta$PATIENT_ID)),
           d = 7,
           n = ncol(exposures_0005),
           x = cbind(1, as.numeric(as.factor(patient.meta$group))-1),
           z = give_z_matrix_from_labels(lbls = patient.meta$PATIENT_ID))
obj_original2 = list(Y = t(sigs_from_phil),
           num_individuals = length(unique(patient.meta$PATIENT_ID)),
           d = 7,
           n = ncol(sigs_from_phil),
           x = cbind(1, as.numeric(as.factor(patient.meta$group))-1),
           z = give_z_matrix_from_labels(lbls = patient.meta$PATIENT_ID))

give_LRchanges_barplot(obj_0005, nrow_facets = 1)
ggsave("../../results/zeros_modelling/paired_LR_changes_new_zeros_0005.pdf", width = 12, height = 2)
give_LRchanges_barplot(obj_original, nrow_facets = 1)

pdf("../../results/zeros_modelling/paired_LR_changes_several_methods.pdf", width = 12, height = 7)
grid.arrange(give_LRchanges_barplot(obj_0005, nrow_facets = 1)+ggtitle('LogR of exposures with zeros from simulation with threshold=0.005'),
             give_LRchanges_barplot(obj_original, nrow_facets = 1)+ggtitle('LogR of exposures from Phil, without zero threshold'),
             give_LRchanges_barplot(obj_original2, nrow_facets = 1)+ggtitle('LogR of exposures from Phil'))
dev.off()

ggplot(give_df_zeros_given_thresh(0.005, perturbed_sc_sigs, sigs_from_phil), aes(x=sig, fill=transition))+geom_bar()+
  scale_fill_manual(values = c('#a6e695', '#95e6ca', '#db728e', '#a769d6'))+
  labs(x='Signature')
ggsave("../../results/zeros_modelling/robustness_zero_from_simulation_0005.pdf", height = 3, width = 6)

plts = lapply(c(0.001, 0.002, 0.005, 0.01, 0.02), function(thresh){
  ggplot(give_df_zeros_given_thresh(thresh, perturbed_sc_sigs, sigs_from_phil),
         aes(x=sig, fill=transition))+geom_bar()+ scale_fill_manual(values = c('#a6e695', '#95e6ca', '#db728e', '#a769d6'))+
  labs(x='Signature')+ggtitle( paste0('Threshold: ', thresh) )
})

pdf("../../results/zeros_modelling/robustness_zero_from_simulation_variousthreshs.pdf", width = 7, height = 6)
do.call('grid.arrange', plts)
dev.off()

##' conclusions:
##' - there were almost no exposures for s1 in original derivation; simulation shows quite a lot
##' - signatures for which the only inconsistency is that there are more zeros in simulation: s1, s2, s4, s5
##' - in s3, there were some zeros in original which have now disappeared
##' - in s6, there are mostly consistencies, in a few inconsistencies in both directions
##' - in s7, all the zeros in the original disappear (making essentially always nonzero in simulation)
zeros_df2_0005 %>% filter(sig == 's7') %>% select(my_zeros) %>% table

## looking at examples of s7 becoming nonzero in simulation
zeros_df2_0005 %>% filter(sig == 's7') %>% filter( (phil_zeros == 'TRUE') & (my_zeros == 'FALSE') )
ggplot(melt(perturbed_sc_sigs$IM_190), aes(x=Var2, y=value))+geom_point()+
  geom_abline(slope = 0, intercept = sigs['s7','IM_190'])+
  geom_abline(slope = 0, intercept = median(perturbed_sc_sigs$IM_190[,'s7']), col='blue')
ggplot(cbind.data.frame(x=perturbed_sc_sigs$IM_190[,7]), aes(x=x))+geom_density()+
  geom_vline(xintercept = sigs['s7','IM_190'])+
  geom_vline(xintercept = median(perturbed_sc_sigs$IM_190[,'s7']), col='blue')+
  geom_vline(xintercept = 0.005, col='red')

## looking at examples of s3 becoming nonzero in simulation

# threshs = c(0.001, 0.002, 0.005, 0.01, 0.02)
# threshs = seq(from=0.0001, to=0.02, by=0.0001)
threshs = seq(from=0.0001, to=0.05, by=0.005)

strategy1_results = mclapply(threshs, give_df_zeros_given_thresh, perturbed_sc_sigs, sigs_from_phil)
names(strategy1_results) = threshs

df_grouped_strategy1 = melt(strategy1_results) %>% group_by(L1, sig) %>%
  mutate(prop_zero_original=sum(transition == 'Zero in original only')/sum(transition %in% c('Zero in original only', 'Consistent zero'))) %>%
  mutate(prop_zero_sim=sum(transition == 'Zero in simulation only')/sum(transition %in% c('Zero in simulation only', 'Consistent zero')))

# pdf("../../results/zeros_modelling/robustness_zero_from_simulation_variousthreshs_2.pdf", width = 9, height = 2)
# grid.arrange(ggplot(df_grouped_strategy1,
# aes(x=as.numeric(L1), y=prop_zero_sim, col=sig, group=sig))+geom_line()+
#   scale_color_brewer(palette="Dark2")+labs(x='Threshold', y='(Zero only in simulation)/(Zero in simulation)')+ggtitle('Simulation #1\n(Zero only in simulation)/(Zero in simulation)'),
# ggplot(df_grouped_strategy1,
#        aes(x=as.numeric(L1), y=prop_zero_original, col=sig, group=sig))+geom_line()+
#   scale_color_brewer(palette="Dark2")+labs(x='Threshold', y='(Zero in original only)/(Zero in original)')+ggtitle('Simulation #1\n(Zero only in original)/(Zero in original)'),
# ggplot(df_grouped_strategy1, aes(x=prop_zero_original, y=1-prop_zero_sim, col=sig, group=sig))+geom_point()+
#   scale_color_brewer(palette="Dark2")+ggtitle('Proportions of zeros in only original and only simulation , for various thresholds'),
# ncol=2)
# dev.off()

ggplot(df_grouped_strategy1,
       aes(x=as.numeric(L1), y=prop_zero_sim, col=sig, group=sig))+geom_line()+
         scale_color_brewer(palette="Dark2")+labs(x='Threshold', y='(Zero only in simulation)/(Zero in simulation)')+ggtitle('Simulation #1\n(Zero only in simulation)/(Zero in simulation)')
ggsave("../../results/zeros_modelling/robustness_zero_from_simulation_variousthreshs_2a.png", width = 5, height = 5)
ggplot(df_grouped_strategy1,
      aes(x=as.numeric(L1), y=prop_zero_original, col=sig, group=sig))+geom_line()+
 scale_color_brewer(palette="Dark2")+labs(x='Threshold', y='(Zero in original only)/(Zero in original)')+ggtitle('Simulation #1\n(Zero only in original)/(Zero in original)')
ggsave("../../results/zeros_modelling/robustness_zero_from_simulation_variousthreshs_2b.png", width = 5, height = 5)
ggplot(df_grouped_strategy1, aes(x=prop_zero_original, y=1-prop_zero_sim, col=sig, group=sig))+geom_point()+
 scale_color_brewer(palette="Dark2")+ggtitle('Proportions of zeros in only original and only simulation , for various thresholds')
ggsave("../../results/zeros_modelling/robustness_zero_from_simulation_variousthreshs_2c.png", width = 5, height = 5)
#-----------------------------------------------------------------------------------------------#


#-----------------------------------------------------------------------------------------------#

### Simulate exposures

# sim covariances
nsamples_sim = 250
sigma = cbind(runif(nsamples_sim, -2, 2), runif(nsamples_sim, -2, 2), runif(nsamples_sim, -2, 2), runif(nsamples_sim, -2, 2))
sigma = cbind(sigma, sigma[,1]+runif(nsamples_sim, 0, 1))
sigma = cbind(sigma, sigma[,2]+runif(nsamples_sim, -1, 1))

sigs_sim = mvtnorm::rmvnorm(nsamples_sim, mean = runif(6), sigma = cov(sigma))
sigs_sim = cbind(sigs_sim, 0)
sigs_sim = t(sweep(exp(sigs_sim), 1, rowSums(exp(sigs_sim)), '/'))
sigs_sim_counts = apply(sigs_sim, 2, function(probs_sigs) rmultinom(n = 1, size = rpois(1, 80), prob = probs_sigs))
colnames(sigs_sim) = paste0('SimSample', 1:ncol(sigs_sim))

## set randomly to zero
sigs_sim[sample(1:length(sigs_sim), size = round(length(sigs_sim)*0.1))] = 0

## simulate sum of posteriors matrix
##' for each sample, get the number of counts attributed to each sample, and select components (for each feature independently),
##' sampling with probability from the feat_sig_mat. Then add all the values for each component in a sample, and normalise feature-wise
##' so that they add to one

select_features_with_prob = function(i){
  .ret = sapply(which_features, function(idx_features){
    sampled_idx = sample(x = 1:length(idx_features), size = 1, replace = F, prob = i[idx_features])
    .x = rep(0, length(idx_features))
    names(.x) = rownames(feat_sig_mat)[idx_features]
    .x[sampled_idx] = 1
    .x
})
  names(.ret) = NULL
  unlist(.ret)
}

# normalise_per_prob = function(i){
#   .ret = sapply(which_features, function(idx_features){
#     i[idx_features] = i[idx_features]/sum(i[idx_features])
#   })
#   names(.ret) = NULL
#   unlist(.ret)
# }

SxCmatrix_sim = apply(sigs_sim_counts, 2, function(vec_counts){
  SoP_from_counts = sapply(1:length(vec_counts),
    function(idx_sig) rowSums(sapply(1:vec_counts[idx_sig], function(sig_count_it) select_features_with_prob(feat_sig_mat[,idx_sig]) )))
  # normalise_per_prob(rowSums(SoP_from_counts))
  rowSums(SoP_from_counts)
}  )

recovery = quantifySignatures(t(SxCmatrix_sim))
plot(unlist(sigs_sim), unlist(recovery), pch=8, cex=0.3)

ggplot(cbind.data.frame(true_sigs=as.vector(sigs_sim), recovery_sigs=as.vector(recovery), sig=rep(paste0('s', 1:7), ncol(sigs_sim))),
       aes(x=true_sigs, y=recovery_sigs))+facet_wrap(.~sig, nrow=1)+geom_point(size=.6)+geom_abline(slope = 1, intercept = 0, lty='dashed')+
  labs(x='True exposures (simulation)', y='Recovered exposures\n(simulation)')
ggsave("../../results/zeros_modelling/recovery_exposures.pdf", width = 10, height = 2)

amounts = c(0.005, 0.01, 0.02, 0.04, 0.1, 0.2)
perturbed_sc_sim = lapply(amounts, function(amount_it){
  print(amount_it)
  .x = lapply(1:ncol(sigs_sim), function(sampl) replicate(n = 10, expr = mod_version(sampl, amount_it, SxCmatrix = t(SxCmatrix_sim)), simplify = F))
  lapply(.x, function(i) do.call('rbind', i))
})
perturbed_sc_sigs_sim = lapply(perturbed_sc_sim, function(samp){lapply(samp, function(i) t(quantifySignatures_alt(i)))})

ggplot(melt(perturbed_sc_sigs_sim[[2]][[1]]), aes(x=Var2, y=value))+geom_point()

strategy1_results_sim = lapply(perturbed_sc_sigs_sim, function(i){
  mclapply(threshs, give_df_zeros_given_thresh, perturbed_sc_sigs=i, sigs_from_phil = sigs_sim)
  })

# perturbed_sc_sigs = i[[1]]
# threshold = threshs[1]
# sigs_from_phil = sigs_sim

names(strategy1_results_sim) = paste0('Amount:', amounts)
for(i in 1:length(strategy1_results_sim)){
  names(strategy1_results_sim[[i]]) = paste0('Thresh:', threshs)
}

df_grouped_strategy1_sim = melt(strategy1_results_sim) %>% group_by(L1, L2, sig) %>%
  mutate(prop_zero_original=sum(transition == 'Zero in original only')/sum(transition %in% c('Zero in original only', 'Consistent zero'))) %>%
  mutate(prop_zero_sim=sum(transition == 'Zero in simulation only')/sum(transition %in% c('Zero in simulation only', 'Consistent zero'))) %>%
  mutate(sensitivity=sum(transition == 'Consistent zero')/sum(transition %in% c('Zero in original only', 'Consistent zero'))) %>%
  mutate(specificity=sum(transition == 'Consistent nonzero')/sum(transition %in% c('Zero in simulation only', 'Consistent nonzero')))

df_grouped_strategy1_sim_poolsig = melt(strategy1_results_sim) %>% group_by(L1, L2) %>%
  mutate(prop_zero_original=sum(transition == 'Zero in original only')/sum(transition %in% c('Zero in original only', 'Consistent zero'))) %>%
  mutate(prop_zero_sim=sum(transition == 'Zero in simulation only')/sum(transition %in% c('Zero in simulation only', 'Consistent zero'))) %>%
  mutate(sensitivity=sum(transition == 'Consistent zero')/sum(transition %in% c('Zero in original only', 'Consistent zero'))) %>%
  mutate(specificity=sum(transition == 'Consistent nonzero')/sum(transition %in% c('Zero in simulation only', 'Consistent nonzero')))%>%
  mutate(precision=sum(transition == 'Consistent zero')/sum(transition %in% c('Zero in simulation only', 'Consistent zero'))) %>%
  mutate(recall=sum(transition == 'Consistent zero')/sum(transition %in% c('Zero in original only', 'Consistent zero')))
  
df_grouped_strategy1_sim_poolsig$F1score = (df_grouped_strategy1_sim_poolsig$precision*df_grouped_strategy1_sim_poolsig$recall)/df_grouped_strategy1_sim_poolsig$precision+df_grouped_strategy1_sim_poolsig$recall
df_grouped_strategy1_sim_poolsig = df_grouped_strategy1_sim_poolsig[!duplicated(df_grouped_strategy1_sim_poolsig[,c('sample', 'L2', 'L1')]),]

ggplot(df_grouped_strategy1_sim,
       aes(x=as.numeric(gsub('Thresh:', '', L2)), y=prop_zero_sim, col=sig, group=sig))+geom_line()+facet_wrap(.~gsub('Amount:', '', L1))+
  scale_color_brewer(palette="Dark2")+labs(x='Threshold', y='(Zero only in simulation)/(Zero in simulation)')+ggtitle('Simulation #1\n(Zero only in simulation)/(Zero in simulation)')
ggsave("../../results/zeros_modelling/robustness_zero_from_simulation_variousthreshs_2a_various_amounts_sim.png", width = 8, height = 5)

ggplot(df_grouped_strategy1_sim,
       aes(x=as.numeric(gsub('Thresh:', '', L2)), y=prop_zero_original, col=sig, group=sig))+geom_line()+facet_wrap(.~L1)+
  scale_color_brewer(palette="Dark2")+labs(x='Threshold', y='(Zero in original only)/(Zero in original)')+
  ggtitle('Simulation #1\n(Zero only in original)/(Zero in original)')

ggplot(df_grouped_strategy1_sim, aes(x=prop_zero_original, y=1-prop_zero_sim, col=as.numeric(gsub('Thresh:', '', L2)),
                                     shape= L1))+
  geom_point()+facet_wrap(.~sig, scales='free')+
  ggtitle('Proportions of zeros in only original and only simulation , for various thresholds')

ggplot(df_grouped_strategy1_sim, aes(x=1-specificity, y=sensitivity, col=as.numeric(gsub('Thresh:', '', L2)),
                                     shape= L1, group=L2))+
  geom_point()+facet_wrap(.~sig, scales='free')+geom_line()+
  ggtitle('Specificity and sensitivity, for various thresholds and amounts')

require(jcolors)
ggplot(df_grouped_strategy1_sim_poolsig, aes(x=1-specificity, y=sensitivity, col=as.numeric(gsub('Thresh:', '', L2)),
                                     shape= L1, group=(L1)))+
  geom_point()+geom_line()+
  # scale_color_continuous(name="chroma")
  scale_color_jcolors_contin("pal3", reverse = TRUE, bias = 2.25)+labs(col = "Threshold", shape="Amount")
ggsave("../../results/zeros_modelling/robustness_zero_from_simulation_variousthreshs_2a_various_amounts_sim_ROC.png", width = 5.5, height = 5)

ggplot(df_grouped_strategy1_sim_poolsig, aes(x=1-specificity, y=sensitivity, col=as.numeric(gsub('Thresh:', '', L2)),
                                             shape= L1, group=(L2)))+
  geom_point()+geom_line()+
  # scale_color_continuous(name="chroma")
  scale_color_jcolors_contin("pal3", reverse = TRUE, bias = 2.25)+labs(col = "Threshold", shape="Amount")+facet_wrap(.~L1)
ggsave("../../results/zeros_modelling/robustness_zero_from_simulation_variousthreshs_2a_various_amounts_sim_ROC_2.png", width = 8, height = 5)


ggplot(df_grouped_strategy1_sim_poolsig, aes(x=1-specificity, y=sensitivity, col=as.numeric(gsub('Amount:', '', L1)),
                                             shape= L1, group=L2))+
  geom_point()+geom_line()+
  # scale_color_continuous(name="chroma")
  scale_color_jcolors_contin("pal3", reverse = TRUE, bias = 2.25)+labs(col = "Threshold", shape="Amount")+facet_wrap(.~L1, scales = "free")

ggplot(df_grouped_strategy1_sim_poolsig, aes(x=as.numeric(gsub('Thresh:', '', L2)), y=as.numeric(gsub('Amount:', '', L1)), fill=F1score))+
  scale_y_continuous(trans='log2')+
  geom_tile()+labs(x='Threshold', y='Amount')+
  scale_fill_jcolors_contin("pal3", reverse = TRUE, bias = 2.25)
ggsave("../../results/zeros_modelling/robustness_zero_from_simulation_variousthreshs_2a_various_amounts_sim_F1.png", width = 5, height = 4)

#-----------------------------------------------------------------------------------------------#

#-----------------------------------------------------------------------------------------------#
### Real data, multiple amounts
perturbed_sc_multiple = lapply(amounts, function(amount_it){
  print(amount_it)
  .x = lapply(1:ncol(sigs), function(sampl) replicate(n = 20, expr = mod_version(sampl, amount_it), simplify = F))
  lapply(.x, function(i) do.call('rbind', i))
})
perturbed_sc_multiple_sigs = lapply(perturbed_sc_multiple, function(samp){lapply(samp, function(i) t(quantifySignatures_alt(i)))})


strategy1_results_multiple = lapply(perturbed_sc_multiple_sigs, function(i){
  mclapply(threshs, give_df_zeros_given_thresh, perturbed_sc_sigs=i, sigs_from_phil = sigs_from_phil)
})
names(strategy1_results_multiple) = paste0('Amount:', amounts)
for(i in 1:length(strategy1_results_multiple)){
  names(strategy1_results_multiple[[i]]) = paste0('Thresh:', threshs)
}
df_grouped_strategy1_multiple = melt(strategy1_results_multiple) %>% group_by(L1, L2, sig) %>%
  mutate(prop_zero_original=sum(transition == 'Zero in original only')/sum(transition %in% c('Zero in original only', 'Consistent zero'))) %>%
  mutate(prop_zero_sim=sum(transition == 'Zero in simulation only')/sum(transition %in% c('Zero in simulation only', 'Consistent zero')))
ggplot(df_grouped_strategy1_multiple,
       aes(x=as.numeric(gsub('Thresh:', '', L2)), y=prop_zero_sim, col=sig, group=sig))+geom_line()+facet_wrap(.~L1)+
  scale_color_brewer(palette="Dark2")+labs(x='Threshold', y='(Zero only in simulation)/(Zero in simulation)')+ggtitle('Simulation #1\n(Zero only in simulation)/(Zero in simulation)')
ggsave("../../results/zeros_modelling/robustness_zero_from_simulation_variousthreshs_2a_various_amounts.png", width = 8, height = 5)


