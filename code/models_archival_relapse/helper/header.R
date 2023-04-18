library(rstan)
library(uuid)
library(ggplot2)
library(optparse)
library(reshape2)
library(TMB)
library(pheatmap)
library(bayesplot)
library(dplyr)

# source("../../../../../GlobalDA/code/2_inference/helper/helper_DA_stan.R")
# source("../../../../../GlobalDA/code/3_analysis/helper/helper_analyse_posteriors.R")
# source("../../../../../GlobalDA/code/2_inference_TMB/helper_TMB.R")
source("../helper/functions.R")

#-------------------------------------------------------------------------------------------#
TMB::compile("../tmb_RE/tmb_RE_lm.cpp","-std=gnu++17")
dyn.load(dynlib("../tmb_RE/tmb_RE_lm"))
TMB::compile("../tmb_RE/tmb_correlated_multinom.cpp","-std=gnu++17")
dyn.load(dynlib("../tmb_RE/tmb_correlated_multinom"))
# TMB::compile("../tmb_RE/tmb_correlated_multinom_1.cpp")
# dyn.load(dynlib("../tmb_RE/tmb_correlated_multinom_1"))
TMB::compile("../tmb_RE/tmb_correlated_multinom_2.cpp","-std=gnu++17")
dyn.load(dynlib("../tmb_RE/tmb_correlated_multinom_2"))
TMB::compile("../tmb_RE/tmb_correlated_multinom_2_allFE.cpp","-std=gnu++17")
dyn.load(dynlib("../tmb_RE/tmb_correlated_multinom_2_allFE"))
#-------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------#
load("../../../data/britroc_30kb_signature_data.rds")
# patient.meta
# sample_by_component
# sig_quants


alpha_val = 1.2

exposures = t(sig_quants)
exposures_inputation1em4 = normalise_rw(exposures + 1e-4)
exposures_alphatrans = t(apply(sig_quants, 2, alphatrans, alpha=alpha_val))
colnames(exposures_alphatrans) = paste0('Transformed', 1:ncol(exposures_alphatrans))
patient.meta = patient.meta[match(rownames(exposures_alphatrans), as.character(patient.meta$SAMPLE_ID)),]
sample_by_component = sample_by_component[match(rownames(exposures_alphatrans), rownames(sample_by_component)),]

exposures_zeros = apply(exposures > 0, 2, as.numeric)

# par(mfrow=c(2,3))
# apply(exposures_alphatrans, 2, hist, breaks=30)

#-------------------------------------------------------------------------------------------#

x = t(cbind(1, as.numeric(factor(patient.meta$group))-1))

d = ncol(exposures_alphatrans) ## number of features
n = nrow(exposures_alphatrans) ## number of samples
num_indiv = length(unique(patient.meta$PATIENT_ID))