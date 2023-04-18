##---------------------------------------------------------------------------------##
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("../../code/models_archival_relapse/run_partialILR/")

set.seed(1234)

folder_images_out <- "../../../results/partialILRmodelling_FE_other_regression/primary_platinum_resistant/"
name_comparison <- "_primaryplatinumsensitivity_"

library(uuid)
library(ggplot2)
library(reshape2)
library(compositions)
library(gridExtra)
library(TMB)
library(xtable)

source("../../../../../CDA_in_Cancer/code/functions/meretricious/pretty_plots/prettySignatures.R")
source("../../../../../GlobalDA/code/2_inference_TMB/helper_TMB.R")
source("../../../../../other_repos/Vias_Brenton/copy_number_analysis_organoids/helper_functions.R")
source("../helper/functions.R")
source("../helper/header.R")

## Modified by PS to drop factors
primary_platinum_resistant <- read.table("../../../data/britroc_primary_platinum_resistant_patient_list.tsv",stringsAsFactors = F)

##---------------------------------------------------------------------------------##


##---------------------------------------------------------------------------------##

## partial ILR without correlations
TMB::compile("../tmb_RE/tmb_MVN_partial_ILR_FEb.cpp", "-std=gnu++17")
dyn.load(dynlib("../tmb_RE/tmb_MVN_partial_ILR_FEb"))

## partial ILR with correlations
TMB::compile("../tmb_RE/tmb_MVN_partial_ILR_FEe.cpp", "-std=gnu++17")
dyn.load(dynlib("../tmb_RE/tmb_MVN_partial_ILR_FEe"))

## MVN with correlations, without zeros
TMB::compile("../tmb_RE/tmb_MVN_ILR.cpp", "-std=gnu++17")
dyn.load(dynlib("../tmb_RE/tmb_MVN_ILR"))

## MVN without correlations, without zeros
TMB::compile("../tmb_RE/tmb_RE_20220222.cpp", "-std=gnu++17")
dyn.load(dynlib("../tmb_RE/tmb_RE_20220222"))

## partial ILR
TMB::compile("../tmb_RE/tmb_MVN_partial_ILR.cpp", "-std=gnu++17")
dyn.load(dynlib("../tmb_RE/tmb_MVN_partial_ILR"))

## partial ILR without correlations
TMB::compile("../tmb_RE/tmb_MVN_partial_ILR_notcor.cpp", "-std=gnu++17")
dyn.load(dynlib("../tmb_RE/tmb_MVN_partial_ILR_notcor"))


##---------------------------------------------------------------------------------##


##---------------------------------------------------------------------------------##
## Added by PS
load("../../../data/britroc_30kb_signature_data.rds")
exposures=t(sig_quants)

primary_platinum_resistant
primary_platinum_resistance = primary_platinum_resistant$pt_group[match(patient.meta$PATIENT_ID, rownames(primary_platinum_resistant))]
primary_platinum_resistance[is.na(primary_platinum_resistance)] <- "sensitive"
primary_platinum_resistance_char <- primary_platinum_resistance
primary_platinum_resistance <- as.numeric(factor(primary_platinum_resistance))
table(primary_platinum_resistance)
table(primary_platinum_resistance_char)

##---------------------------------------------------------------------------------##

##---------------------------------------------------------------------------------##
## INPUT

TMB_data_with_subset_for_nlminb_imput_sensitivity <- prepare_TMB_data_with_subset(exposures = impute(exposures, 1e-2), 1:7,
                                                                                  .keep_additional=!(is.na(primary_platinum_resistance)),
                                                                                  add_x_vec = primary_platinum_resistance-1)
dim(TMB_data_with_subset_for_nlminb_imput_sensitivity$Y)
TMB_data_with_subset_for_nlminb_imput_sensitivity$z <- diag(nrow(TMB_data_with_subset_for_nlminb_imput_sensitivity$z)) ## need to keep this for size of ularge
TMB_data_with_subset_for_nlminb_imput_sensitivity$num_individuals <- NULL
stopifnot(sum(is.na(TMB_data_with_subset_for_nlminb_imput_sensitivity$x[,2])) == 0)

TMB_data_with_subset_for_nlminb_imput_sensitivity_alr <- prepare_TMB_data_with_subset(exposures = impute(exposures, 1e-2), 1:7,
                                                                                      .keep_additional=!(is.na(primary_platinum_resistance)),
                                                                                      add_x_vec = primary_platinum_resistance-1, ilr_trans = F)
TMB_data_with_subset_for_nlminb_imput_sensitivity_alr$Y <- as(compositions::alr(impute(exposures , 1e-2)), 'matrix')
TMB_data_with_subset_for_nlminb_imput_sensitivity_alr$z <- diag(nrow(TMB_data_with_subset_for_nlminb_imput_sensitivity_alr$z)) ## need to keep this for size of ularge
TMB_data_with_subset_for_nlminb_imput_sensitivity_alr$num_individuals <- NULL

## some sensitivities were NA
TMB_data_with_subset_for_nlminb_imput_sensitivity_alr <- give_subset_samples_TMBobj_2(TMB_data_with_subset_for_nlminb_imput_sensitivity_alr, selected_rows_obs = which(!is.na(TMB_data_with_subset_for_nlminb_imput_sensitivity_alr$x[,2])))
stopifnot(sum(is.na(TMB_data_with_subset_for_nlminb_imput_sensitivity_alr$x[,2])) == 0)

run_model_sensitivity <- function(correlation_model_status, TMB_data=TMB_data_with_subset_for_nlminb_imput_sensitivity, add='', imput=T, FEMEbool='FE'){
  if(correlation_model_status == 'FEnocor'){
    model_arg <- 'tmb_MVN_partial_ILR_FEb'
  }else if(correlation_model_status == 'FEcor'){
    model_arg <- 'tmb_MVN_partial_ILR_FEe'
  }else if(correlation_model_status == "FEcorALR"){
    model_arg <- 'tmbMVNILR'
  }else if(correlation_model_status == "FEnocorALR"){
    model_arg <- 'tmbMVNILRnocor' ## 'tmb_RE_20220222'
  }else  if(correlation_model_status == 'partialILRcor'){
    model_arg <- 'partialILR'
  }else  if(correlation_model_status == 'partialILRnocor'){
    model_arg <- 'partialILRnocor'
  }else{
    stop("Invalid model")
  }
  
  if(imput){
    imput_add = '_imput'
  }else{
    imput_add = ''
  }
  
  cat('Model: ', model_arg, '\n')
  
  res_nlminb_Fed2_imput_sensitivity <- wrapper_run_TMB_use_nlminb(model = model_arg,
                                                                  object = TMB_data,
                                                                  use_nlminb = T)
  saveRDS(res_nlminb_Fed2_imput_sensitivity, paste0("../../../data/inference/", FEMEbool, name_comparison, correlation_model_status, add, imput_add, ".RDS"))
  
  if(ncol(TMB_data$x) == 2){
    plot_betas(res_nlminb_Fed2_imput_sensitivity)
    ggsave(paste0(folder_images_out, "betas", name_comparison, correlation_model_status, add, imput_add, ".pdf"), width = 4, height = 2.5)
    
    simulated_vals_sensitivity <- melt(sim_LNM_FE_with_sd(covariates=cbind(1, rep(c(0,1), each=200)),
                                                          TMB_res=res_nlminb_Fed2_imput_sensitivity))
    simulated_vals_sensitivity$Var1b <- as.numeric(gsub("\\..*", "", gsub("^X", "", simulated_vals_sensitivity$Var1)))
    simulated_vals_sensitivity$Var1b[simulated_vals_sensitivity$Var1b == 0] <- 'resistant'
    simulated_vals_sensitivity$Var1b[simulated_vals_sensitivity$Var1b == 1] <- 'sensitive'
    ggplot((simulated_vals_sensitivity), aes(x=(Var1),
                                             col=Var2, y=value))+geom_bar(stat = "identity")+facet_wrap(.~Var1b, scales = "free_x")+
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
            axis.ticks.x=element_blank())
    ggsave(paste0(folder_images_out, "sim", name_comparison, correlation_model_status, add, imput_add, ".pdf"), width = 7, height = 4.5)
    
    exposures_sensitivity <- as(compositions::ilrInv( TMB_data$Y), 'matrix'); rownames(exposures_sensitivity) <- make.names(TMB_data$x[,2], unique = T); colnames(exposures_sensitivity) <- paste0('s', colnames(exposures_sensitivity))
    
    fitted_vals_sensitivity <- sim_LNM_FE(covariates=TMB_data$x,
                                          TMB_res=res_nlminb_Fed2_imput_sensitivity)
    
    
    ggplot(dcast(melt(list(fitted_vals_sensitivity,
                           exposures_sensitivity)), Var1+Var2~L1, value.var = "value"),
           aes(x=`1`,`2`, group=`1`, col=as.numeric(gsub("\\..*", "", gsub("^X", "", Var1)))))+
      geom_violin()+
      geom_jitter(alpha=0.2)+
      facet_wrap(.~Var2, scales = "free")+
      geom_abline(slope = 1, intercept = 0)+theme_bw()+labs(x='Fitted values', y='Observed values', col='Sensitivity (1)')
    ggsave(paste0(folder_images_out, "fit_and_obs_scatter", name_comparison, correlation_model_status, add, imput_add, ".pdf"), width = 7, height = 7.5)
  }else{
    cat('Multiple covariates. Plots not rendered\n')
  }
}

run_model_sensitivity('FEnocor')
run_model_sensitivity('FEcor')

## with ALR, instead of ILR, transformation
give_pairs_with_mvn_wrapper(TMB_data_with_subset_for_nlminb_imput_sensitivity_alr$Y)
run_model_sensitivity('FEcorALR', TMB_data = TMB_data_with_subset_for_nlminb_imput_sensitivity_alr, add='_alr')
run_model_sensitivity('FEnocorALR', TMB_data = TMB_data_with_subset_for_nlminb_imput_sensitivity_alr, add='_alr')

python_like_grep <- function(i, j){
  i[grep(j, i)]
}


## comparing ILR and ALR
par(mfrow=c(1,1), mar=c(6,6,2,2))
plot(as.vector(TMB_data_with_subset_for_nlminb_imput_sensitivity$Y),
     as.vector(TMB_data_with_subset_for_nlminb_imput_sensitivity_alr$Y[match(rownames(TMB_data_with_subset_for_nlminb_imput_sensitivity$Y),
                                                                             rownames(TMB_data_with_subset_for_nlminb_imput_sensitivity_alr$Y)),]))
##---------------------------------------------------------------------------------##

##---------------------------------------------------------------------------------##
## Stratifying diagnosis and relapse patients

give_alr_dataset_selected <- function(group_arg){
  TMB_data_with_subset_for_nlminb_imput_sensitivity_alr <- prepare_TMB_data_with_subset(exposures = impute(exposures, 1e-2), 1:7,
                                                                                        .keep_additional= ( (!(is.na(primary_platinum_resistance))) & (patient.meta$group == group_arg)),
                                                                                        add_x_vec = primary_platinum_resistance-1, ilr_trans = F)
  TMB_data_with_subset_for_nlminb_imput_sensitivity_alr$Y <- as(compositions::alr(impute(exposures , 1e-2)), 'matrix')
  TMB_data_with_subset_for_nlminb_imput_sensitivity_alr$z <- diag(nrow(TMB_data_with_subset_for_nlminb_imput_sensitivity_alr$z)) ## need to keep this for size of ularge
  TMB_data_with_subset_for_nlminb_imput_sensitivity_alr$num_individuals <- NULL
  ## some sensitivities were NA
  TMB_data_with_subset_for_nlminb_imput_sensitivity_alr <- give_subset_samples_TMBobj_2(TMB_data_with_subset_for_nlminb_imput_sensitivity_alr, selected_rows_obs = which(!is.na(TMB_data_with_subset_for_nlminb_imput_sensitivity_alr$x[,2])))
  stopifnot(sum(is.na(TMB_data_with_subset_for_nlminb_imput_sensitivity_alr$x[,2])) == 0)
  TMB_data_with_subset_for_nlminb_imput_sensitivity_alr
}

TMB_data_with_subset_for_nlminb_imput_sensitivity_alr_rlps <- give_alr_dataset_selected('rlps')
TMB_data_with_subset_for_nlminb_imput_sensitivity_alr_arx <- give_alr_dataset_selected('arx')
dim(TMB_data_with_subset_for_nlminb_imput_sensitivity_alr_rlps$Y)
dim(TMB_data_with_subset_for_nlminb_imput_sensitivity_alr_arx$Y)

run_model_sensitivity('FEcorALR', TMB_data = TMB_data_with_subset_for_nlminb_imput_sensitivity_alr_arx, add='_alrARX')
run_model_sensitivity('FEcorALR', TMB_data = TMB_data_with_subset_for_nlminb_imput_sensitivity_alr_rlps, add='_alrRLPS')
run_model_sensitivity('FEnocorALR', TMB_data = TMB_data_with_subset_for_nlminb_imput_sensitivity_alr_arx, add='_alrARXnocor')
run_model_sensitivity('FEnocorALR', TMB_data = TMB_data_with_subset_for_nlminb_imput_sensitivity_alr_rlps, add='_alrRLPSnocor')

##---------------------------------------------------------------------------------##


##---------------------------------------------------------------------------------##

##---------------------------------------------------------------------------------##
## MIXED EFFECTS PARTIALILR WITH MULTIPLE PATIENT SAMPLES

TMB_data_with_subset_for_nlminb_sensitivity_partialILR <- prepare_TMB_data_with_subset(exposures = exposures, 1:7,
                                                                                       .keep_additional=!(is.na(primary_platinum_resistance)),
                                                                                       add_x_vec = primary_platinum_resistance-1)

image(TMB_data_with_subset_for_nlminb_sensitivity_partialILR$x)
image(TMB_data_with_subset_for_nlminb_sensitivity_partialILR$z)
run_model_sensitivity('partialILRcor',
                      TMB_data = TMB_data_with_subset_for_nlminb_sensitivity_partialILR,
                      add='_partialILR', imput = F, FEMEbool = "ME")
run_model_sensitivity('partialILRnocor',
                      TMB_data = TMB_data_with_subset_for_nlminb_sensitivity_partialILR,
                      add='_partialILRnocor', imput = F, FEMEbool = "ME")

partialILR = readRDS(paste0("../../../data/inference/ME", name_comparison, "partialILRcor", "_partialILR", ".RDS"))
partialILR

TMB_data_with_subset_for_nlminb_sensitivity_partialILR_two_groups = TMB_data_with_subset_for_nlminb_sensitivity_partialILR
TMB_data_with_subset_for_nlminb_sensitivity_partialILR_ARX = give_subset_samples_TMBobj_2(TMB_data_with_subset_for_nlminb_sensitivity_partialILR,
                                                                                          which(sapply(rownames(TMB_data_with_subset_for_nlminb_sensitivity_partialILR$Y), function(i) patient.meta[patient.meta$SAMPLE_ID == i,'group'] == "arx")))
TMB_data_with_subset_for_nlminb_sensitivity_partialILR_RLPS = give_subset_samples_TMBobj_2(TMB_data_with_subset_for_nlminb_sensitivity_partialILR,
                                                                                          which(sapply(rownames(TMB_data_with_subset_for_nlminb_sensitivity_partialILR$Y), function(i) patient.meta[patient.meta$SAMPLE_ID == i,'group'] == "rlps")))
dim(TMB_data_with_subset_for_nlminb_sensitivity_partialILR_two_groups$Y)
dim(TMB_data_with_subset_for_nlminb_sensitivity_partialILR_ARX$Y)
dim(TMB_data_with_subset_for_nlminb_sensitivity_partialILR_RLPS$Y)


TMB_data_with_subset_for_nlminb_sensitivity_partialILR_two_groups$x = cbind(TMB_data_with_subset_for_nlminb_sensitivity_partialILR_two_groups$x,
                                                                            as.numeric(factor(patient.meta$group[match(rownames(TMB_data_with_subset_for_nlminb_sensitivity_partialILR_two_groups$Y), patient.meta$SAMPLE_ID)]))-1)
image(TMB_data_with_subset_for_nlminb_sensitivity_partialILR_two_groups$x)
run_model_sensitivity(correlation_model_status = 'partialILRnocor',
                      TMB_data = TMB_data_with_subset_for_nlminb_sensitivity_partialILR_two_groups,
                      add='_partialILRnocor_group_sensitivity',
                      imput = F,
                      FEMEbool = "ME")

run_model_sensitivity(correlation_model_status = 'partialILRnocor',
                      TMB_data = TMB_data_with_subset_for_nlminb_sensitivity_partialILR_ARX,
                      add='_partialILRnocor_group_sensitivityARX',
                      imput = F,
                      FEMEbool = "ME")


run_model_sensitivity(correlation_model_status = 'partialILRnocor',
                      TMB_data = TMB_data_with_subset_for_nlminb_sensitivity_partialILR_RLPS,
                      add='_partialILRnocor_group_sensitivityRLPS',
                      imput = F,
                      FEMEbool = "ME")

# correlation_model_status = 'partialILRnocor'
# TMB_data = TMB_data_with_subset_for_nlminb_sensitivity_partialILR_two_groups
# add='_partialILRnocor_group_sensitivity'
# imput = F
# FEMEbool = "ME"
# 
# model = model_arg
# object = TMB_data
# use_nlminb = T
# initial_params = NULL

partialILR_twocovar = readRDS(paste0("../../../data/inference/ME", name_comparison, "partialILRnocor", "_partialILRnocor_group_sensitivity", ".RDS"))

plot(plot_betas(TMB_obj = partialILR_twocovar, p=3, return_plot = T,
                names_covariates=c('Intercept', 'Sensitive/resistant', 'Diagnosis/relapse')))
ggsave(paste0(folder_images_out, "betas_sensitivity_and_group_", "partialILRnocor", "_partialILRnocor_group", name_comparison, ".pdf"), width = 6, height = 2.5)

##---------------------------------------------------------------------------------##

##---------------------------------------------------------------------------------##
dev.off()
pdf(paste0(folder_images_out, "ternary_sensitivity.pdf"), width = 6, height = 2.5)
par(mfrow=c(1,2), mar=c(0,0,2,0))
plot_ternary(exposures[which(primary_platinum_resistance == 1),c(1,3,7)], main='Resistant', legend_on = F)
plot_ternary(exposures[which(primary_platinum_resistance == 2),c(1,3,7)], main='Sensitive', legend_on = F)
dev.off()


##---------------------------------------------------------------------------------##

##---------------------------------------------------------------------------------##
save.image(paste0("../../../data/inference/image_primary_platinum_resistant.RDS"))
##---------------------------------------------------------------------------------##

##---------------------------------------------------------------------------------##
files_results_model <- sort(python_like_grep(list.files(paste0("../../../data/inference/"), full.names = T), paste0(name_comparison)))
res_nlminb <- sapply(files_results_model, readRDS)
names(res_nlminb) <- gsub(".RDS", "", basename(files_results_model))

# wald_TMB_wrapper(res_nlminb[[paste0('FE', name_comparison, 'FEcor_imput')]])
# wald_TMB_wrapper(res_nlminb[[paste0('FE', name_comparison, 'FEnocor_imput')]])
# wald_TMB_wrapper(res_nlminb[[paste0('FE', name_comparison, 'FEnocorALR_alr_imput')]])

print(xtable::xtable(as(signif(sapply(res_nlminb, wald_TMB_wrapper),
                               4), 'matrix'), digits = 10))

##---------------------------------------------------------------------------------##
