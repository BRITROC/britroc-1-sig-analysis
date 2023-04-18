#-----------------------------------------------------------------------------------------------#
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source("../functions.R")

segs = read.table("../../data/britroc_30kb_ds_absCopyNumber_segmentTable.tsv", header = T)
library(reshape2)
library(ggplot2)
library(dplyr)
library(gridExtra)
source("../functions.R")
source("../../../britroc-cnsignatures-bfb69cd72c50/main_functions.R")
source("../../../britroc-cnsignatures-bfb69cd72c50/helper_functions.R")
source("../../../../GlobalDA/code/2_inference_TMB/helper_TMB.R")
components = readRDS("../../../britroc-cnsignatures-bfb69cd72c50/data/component_parameters.rds")
feat_sig_mat = readRDS("../../../britroc-cnsignatures-bfb69cd72c50/data/feat_sig_mat.rds")
# phils_exposures = t(sig_quants) ## this are phil's exposures
load("../../data/britroc_30kb_signature_data.rds") ## patient meta
SxCmatrix = readRDS("../../out/signature_extraction/SxCmatrix.RDS")
load("../../data/britroc_30kb_signature_data.rds")

library(ggplot2)
library(reshape2)
library(pheatmap)
library(gridExtra)
source("functions.R")

load("../../data/britroc_30kb_signature_data.rds")
# patient.meta
# sample_by_component
# sig_quants

exposures = t(sig_quants)
patient.meta = patient.meta[match(rownames(exposures), as.character(patient.meta$SAMPLE_ID)),]

sigs = quantifySignatures_alt(SxCmatrix)


sigs_from_phil = readRDS("../../data/britroc_30kb_signatures.rds")
#-----------------------------------------------------------------------------------------------#

obj = list(Y = exposures,
           num_individuals = length(unique(patient.meta$PATIENT_ID)),
           d = 7,
           n = nrow(exposures),
           x = cbind(1, as.numeric(as.factor(patient.meta$group))-1),
           z = give_z_matrix_from_labels(lbls = patient.meta$PATIENT_ID))
obj_original = list(Y = t(sigs),
                    num_individuals = length(unique(patient.meta$PATIENT_ID)),
                    d = 7,
                    n = ncol(sigs),
                    x = cbind(1, as.numeric(as.factor(patient.meta$group))-1),
                    z = give_z_matrix_from_labels(lbls = patient.meta$PATIENT_ID))
obj_original2 = list(Y = t(sigs_from_phil),
                     num_individuals = length(unique(patient.meta$PATIENT_ID)),
                     d = 7,
                     n = ncol(sigs_from_phil),
                     x = cbind(1, as.numeric(as.factor(patient.meta$group))-1),
                     z = give_z_matrix_from_labels(lbls = patient.meta$PATIENT_ID))

grid.arrange(give_LRchanges_barplot(obj, nrow_facets = 1),
             give_LRchanges_barplot(obj_original, nrow_facets = 1)+ggtitle('LogR of exposures from Phil, without zero threshold'),
             give_LRchanges_barplot(obj_original2, nrow_facets = 1)+ggtitle('LogR of exposures from Phil'))

give_LRchanges_barplot(obj, nrow_facets = 1)

exposures_old = exposures
obj_old = obj
####-----

# rm(list = ls())

## Infence using stan_fit_LNM for data in the simplex
## Simple LN, no fixed or random effects

# rm(list = ls())
# setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(ggplot2)
library(reshape2)
library(pheatmap)
library(gridExtra)
source("../functions.R")

load("../../data/britroc_30kb_signature_data.rds")
# patient.meta
# sample_by_component
# sig_quants

max(table(patient.meta$SAMPLE_ID))
max(table(patient.meta$PATIENT_ID))

exposures2 = t(sig_quants)
patient.meta = patient.meta[match(rownames(exposures2), as.character(patient.meta$SAMPLE_ID)),]
sample_by_component = sample_by_component[match(rownames(exposures2), rownames(sample_by_component)),]

source("../header.R")
obj = list(Y = exposures2,
           num_individuals = length(unique(patient.meta$PATIENT_ID)),
           d = 7,
           n = nrow(exposures2),
           x = cbind(1, as.numeric(as.factor(patient.meta$group))-1),
           z = give_z_matrix_from_labels(lbls = patient.meta$PATIENT_ID))

grid.arrange(give_LRchanges_barplot(obj, nrow_facets = 1),
give_LRchanges_barplot(obj_old, nrow_facets = 1))

exposures_old
exposures2
t(sigs_from_phil)

#-----------------------------------------------------------------------------------------------#


exposures_old
exposures2
t(sigs_from_phil)

obj_new = list(Y = exposures,
           num_individuals = length(unique(patient.meta$PATIENT_ID)),
           d = 7,
           n = nrow(exposures),
           x = cbind(1, as.numeric(as.factor(patient.meta$group))-1),
           z = give_z_matrix_from_labels(lbls = patient.meta$PATIENT_ID))

grid.arrange(give_LRchanges_barplot(obj, nrow_facets = 1),
             give_LRchanges_barplot(obj_new, nrow_facets = 1),
             give_LRchanges_barplot(obj_old, nrow_facets = 1))


load_one_set_of_exposures = function(){
  #-----------------------------------------------------------------------------------------------#
  rm(list = ls())
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  
  source("../functions.R")
  
  segs = read.table("../../data/britroc_30kb_ds_absCopyNumber_segmentTable.tsv", header = T)
  library(reshape2)
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
  source("../functions.R")
  source("../../../britroc-cnsignatures-bfb69cd72c50/main_functions.R")
  source("../../../britroc-cnsignatures-bfb69cd72c50/helper_functions.R")
  source("../../../../GlobalDA/code/2_inference_TMB/helper_TMB.R")
  components = readRDS("../../../britroc-cnsignatures-bfb69cd72c50/data/component_parameters.rds")
  feat_sig_mat = readRDS("../../../britroc-cnsignatures-bfb69cd72c50/data/feat_sig_mat.rds")
  # phils_exposures = t(sig_quants) ## this are phil's exposures
  load("../../data/britroc_30kb_signature_data.rds") ## patient meta
  SxCmatrix = readRDS("../../out/signature_extraction/SxCmatrix.RDS")
  load("../../data/britroc_30kb_signature_data.rds")
  return(exposures)
}


load_the_other = function(){
  setwd("/Users/morril01/Documents/PhD/other_repos/britroc-1/code/")
  library(ggplot2)
  library(reshape2)
  library(pheatmap)
  library(gridExtra)
  source("functions.R")
  
  load("../data/britroc_30kb_signature_data.rds")
  # patient.meta
  # sample_by_component
  # sig_quants
  
  max(table(patient.meta$SAMPLE_ID))
  max(table(patient.meta$PATIENT_ID))
  
  exposures = t(sig_quants)
  source("header.R")
  return(exposures)
}

xxx1 = load_one_set_of_exposures()
xxx2 = load_the_other()
xxx1
xxx2

all(xxx1 == xxx2)

all(xxx1 == exposures)
obj
obj_new
saveRDS(object = exposures, "~/Desktop/good_exposures.RDS")
saveRDS(object = patient, "~/Desktop/good_patients.RDS")