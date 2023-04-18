## FE using mixed effects code

rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("../../models/run_partialILR/")

# set.seed(1234)

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

tp53 <- read.csv("../../../data/restricted/tp53_opt_panel.csv")
meta_archival <- read.csv("../../../data/restricted/Archival_samples.csv")
clinbrit <- read.csv("../../../data/restricted/clinical/britroc_cohort_patient_data.tsv", sep = "\t")
clinbrit2 <- read.csv("../../../data/restricted/clinical/britroc_patient_treatment_data.tsv", sep = "\t")
# View(tp53)

folder_images_out <- "../../../results/partialILRmodelling_FE_other_regression/"
system(paste0('mkdir -p ', folder_images_out))

yes_TP53 <- tp53$ID[unique(tp53$ID) %in% rownames(exposures)]
no_TP53 <- tp53$ID[!(unique(tp53$ID) %in% rownames(exposures))]

yes_TP53

par(mfrow=c(1,2), mar=c(0,0,0,0))
plot_ternary(exposures[rownames(exposures) %in% yes_TP53,c(1,3,5)], main='TP53', legend_on = F)
plot_ternary(exposures[rownames(exposures) %in% no_TP53,c(1,3,5)], main='No TP53', legend_on = F)

grid.arrange(createBarplot(exposures[rownames(exposures) %in% yes_TP53,])+theme(legend.position="bottom"),
             createBarplot(exposures[rownames(exposures) %in% no_TP53,])+theme(legend.position="bottom"),
             nrow=1)

table(meta_archival$IM_ID)
locations <- unique(meta_archival$Location.of.Tumour.biopsied)
locations <- locations[locations != ""]

dev.off()
pdf("../../exploratory_analyses/signatures_by_location.pdf", width=10)
par(mfrow=c(4,4), mar=c(2,0,1,0))
for(i in locations){
  idx_samples <- exposures[rownames(exposures) %in% meta_archival$IM_ID[meta_archival$Location.of.Tumour.biopsied == i],
                           c(1,3,5)]
  if(length(idx_samples) > 3){
    plot_ternary(idx_samples,
                 main=paste0(i, ' (', nrow(idx_samples), ')'), legend_on = F)
  }
}
dev.off()

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


##---------------------------------------------------------------------------------##
## Age
patient.meta$PATIENT_ID[match(rownames(exposures), patient.meta$SAMPLE_ID)]
age_x = clinbrit$age[match(patient.meta$PATIENT_ID, paste0('BRITROC-', clinbrit$britroc_number))]
TMB_data_with_subset_for_nlminb <- prepare_TMB_data_with_subset(exposures = exposures, subset_sigs = c(1:4, 5, 6, 7), .keep_additional=!(is.na(age_x)), add_x_vec = age_x)
TMB_data_with_subset_for_nlminb_imput <- prepare_TMB_data_with_subset(exposures = impute(exposures, 1e-2), c(1:4, 5, 6, 7),
                                                                      .keep_additional=!(is.na(age_x)), add_x_vec = age_x)
TMB_data_with_subset_for_nlminb$x
TMB_data_with_subset_for_nlminb$z <- diag(nrow(TMB_data_with_subset_for_nlminb$z)) ## need to keep this for size of ularge
TMB_data_with_subset_for_nlminb_imput$z <- diag(nrow(TMB_data_with_subset_for_nlminb_imput$z)) ## need to keep this for size of ularge
TMB_data_with_subset_for_nlminb$num_individuals <- NULL
TMB_data_with_subset_for_nlminb_imput$num_individuals <- NULL

# res_nlminb_Fed <- wrapper_run_TMB_use_nlminb(model = "tmb_MVN_partial_ILR_FEd_integrateularge",
#                                                              object = TMB_data_with_subset_for_nlminb,
#                                                              use_nlminb = T)
# plot_betas(res_nlminb_Fed)
res_nlminb_Fed2 <- wrapper_run_TMB_use_nlminb(model = "tmb_MVN_partial_ILR_FEb",
                                             object = TMB_data_with_subset_for_nlminb,
                                             use_nlminb = T)
res_nlminb_Fed2cor <- wrapper_run_TMB_use_nlminb(model = "tmb_MVN_partial_ILR_FEe",
                                              object = TMB_data_with_subset_for_nlminb,
                                              use_nlminb = T)
grid.arrange(plot_betas(res_nlminb_Fed2, title='no cor'), plot_betas(res_nlminb_Fed2cor, title = 'cor'))

# res_nlminb_Fed_imput <- wrapper_run_TMB_use_nlminb(model = "tmb_MVN_partial_ILR_FEd_integrateularge",
#                                              object = TMB_data_with_subset_for_nlminb_imput,
#                                              use_nlminb = T)

res_nlminb_Fed2_imput <- wrapper_run_TMB_use_nlminb(model = "tmb_MVN_partial_ILR_FEb",
                                                   object = TMB_data_with_subset_for_nlminb_imput,
                                                   use_nlminb = T)

res_nlminb_Fed2_imputcor <- wrapper_run_TMB_use_nlminb(model = "tmb_MVN_partial_ILR_FEe",
                                                    object = TMB_data_with_subset_for_nlminb_imput,
                                                    use_nlminb = T)
grid.arrange(plot_betas(res_nlminb_Fed2_imput, title='no cor'), plot_betas(res_nlminb_Fed2_imputcor, title = 'cor'))

plot_betas(res_nlminb_Fed_imput)
plot_betas(res_nlminb_Fed2_imput, title = "Betas FEb with imput of 1e-2")
ggsave(paste0(folder_images_out, "betas_age_FEnocor_imput.pdf"), width = 7, height = 4.5)


grid.arrange(plot_betas(res_nlminb_Fed, title = 'partial ILR'),
             plot_betas(res_nlminb_Fed2, title = 'partial ILR no cor'),
             plot_betas(res_nlminb_Fed_imput, title = 'with imput'),
             plot_betas(res_nlminb_Fed2_imput, title = 'partial ILR no cor with imput')
)


fitted_vals <- sim_LNM_FE(covariates=TMB_data_with_subset_for_nlminb_imput$x,
                          TMB_res=res_nlminb_Fed2_imput)
simulated_vals <- sim_LNM_FE(covariates=cbind(1, seq(20, 100, 0.2)),
                          TMB_res=res_nlminb_Fed2_imput)
ggplot(melt(fitted_vals), aes(x=factor(Var1, levels=sort(rownames(fitted_vals))), ## nice!
                              col=Var2, y=value))+geom_bar(stat = "identity")+
  geom_vline(xintercept = factor("X80"))+
  # theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+labs(y='Exposures', fill='')

ggplot(melt(simulated_vals), aes(x=factor(Var1, levels=sort(rownames(fitted_vals))),
                              col=Var2, y=value))+geom_bar(stat = "identity")+
  geom_vline(xintercept = factor("X80"))

exposures_age <- as(compositions::ilrInv( TMB_data_with_subset_for_nlminb_imput$Y), 'matrix'); rownames(exposures_age) <- make.names(TMB_data_with_subset_for_nlminb_imput$x[,2], unique = T); colnames(exposures_age) <- paste0('s', colnames(exposures_age))
rownames(exposures_age) <- rownames(TMB_data_with_subset_for_nlminb_imput$Y)
uniq_ages <- sort(unique(c(rownames(fitted_vals), rownames(simulated_vals), rownames(exposures_age))))
ggplot((melt(list(fitted_vals=fitted_vals, simulated_vals=simulated_vals, true_exposures=exposures_age))),
       aes(x=factor(Var1, levels=uniq_ages),
                                 col=Var2, y=value))+geom_bar(stat = "identity")+
  geom_vline(xintercept = which(uniq_ages == "X40"))+
  geom_text(aes(x=which(uniq_ages == "X40"), label="40 years", y=0.5), colour="black", angle=90, vjust = 1.2)+
  geom_vline(xintercept = which(uniq_ages == "X60"))+
  geom_text(aes(x=which(uniq_ages == "X60"), label="60 years", y=0.5), colour="black", angle=90, vjust = 1.2)+
  geom_vline(xintercept = which(uniq_ages == "X80"))+
  geom_text(aes(x=which(uniq_ages == "X80"), label="80 years", y=0.5), colour="black", angle=90, vjust = 1.2)+
  facet_wrap(.~L1, ncol=1)+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+labs(y='Exposures', fill='')
ggsave(paste0(folder_images_out, "fit_and_sim_age_FEnocor_imput.pdf"), width = 7, height = 7.5)

par(mfrow=c(1,1), mar=c(4,3,2,2))
plot(TMB_data_with_subset_for_nlminb_imput$x[,2], exposures_age[,4])

dim(exposures)
dim(exposures_age)

zero_s4 <- exposures[match(rownames(exposures_age), rownames(exposures)),'s4'] == 0
table(zero_s4)
plot(TMB_data_with_subset_for_nlminb_imput$x[,2][!zero_s4], exposures_age[,4][!zero_s4])

system(paste0('open ', folder_images_out))

par(mfrow=c(1,1))
ggplot(data.frame(x=TMB_data_with_subset_for_nlminb_imput$x[,2], y=exposures_age[,'s4']), aes(x=x,y=y))+
  geom_point()+geom_smooth()

rownames(fitted_vals) == rownames(exposures_age)
ggplot(dcast(melt(list(fitted_vals,
     exposures_age)), Var1+Var2~L1, value.var = "value"), aes(x=`1`,`2`,
     col=as.numeric(gsub("\\..*", "", gsub("^X", "", Var1)))))+geom_point()+facet_wrap(.~Var2)+
  geom_abline(slope = 1, intercept = 0)+theme_bw()+labs(x='Fitted values', y='Observed values', col='Age')
ggsave(paste0(folder_images_out, "fit_and_obs_scatter_age_FEnocor_imput.pdf"), width = 7, height = 7.5)

exposures_age <- as(compositions::ilrInv( TMB_data_with_subset_for_nlminb_imput$Y), 'matrix'); rownames(exposures_age) <- make.names(TMB_data_with_subset_for_nlminb_imput$x[,2], unique = T); colnames(exposures_age) <- paste0('s', colnames(exposures_age))

ggplot(dcast(melt(list(fitted_vals,
                       exposures_age)), Var1+Var2~L1, value.var = "value"),
       aes(x=`1`,`2`, group=`1`, col=as.numeric(gsub("\\..*", "", gsub("^X", "", Var1)))))+
  geom_point(alpha=0.2)+ facet_wrap(.~Var2, scales = "free")+
  geom_abline(slope = 1, intercept = 0)+theme_bw()+
  labs(x='Fitted values', y='Observed values', col='Status')

##---------------------------------------------------------------------------------##

clinbrit$drugs
table(clinbrit$pt_sensitivity_at_reg)
all(rownames(exposures) == patient.meta$SAMPLE_ID )
factor(clinbrit$pt_sensitivity_at_reg)
as.numeric(factor(clinbrit$pt_sensitivity_at_reg))-1

##---------------------------------------------------------------------------------##
##' I guess status=0 means alive and status=1 means dead, in which case s3, s5 and s7 are good prognostic
##' markers and s1/s2, s4 and s6 are bad ones
clinbrit$status
status_x = clinbrit$status[match(patient.meta$PATIENT_ID, paste0('BRITROC-', clinbrit$britroc_number))]
TMB_data_with_subset_for_nlminb_imput_status <- prepare_TMB_data_with_subset(exposures = impute(exposures, 1e-2), 1:7,
                                                                                  .keep_additional=!(is.na(status_x)),
                                                                                  add_x_vec = status_x)
TMB_data_with_subset_for_nlminb_imput_status$z <- diag(nrow(TMB_data_with_subset_for_nlminb_imput$z)) ## need to keep this for size of ularge
TMB_data_with_subset_for_nlminb_imput_status$num_individuals <- NULL

res_nlminb_Fed2_imput_status <- wrapper_run_TMB_use_nlminb(model = "tmb_MVN_partial_ILR_FEb",
                                                                object = TMB_data_with_subset_for_nlminb_imput_status,
                                                                use_nlminb = T)
plot_betas(res_nlminb_Fed2_imput_status)
wald_TMB_wrapper(res_nlminb_Fed2_imput_status)

exposures_status <- as(compositions::ilrInv( TMB_data_with_subset_for_nlminb_imput_status$Y), 'matrix'); rownames(exposures_status) <- make.names(TMB_data_with_subset_for_nlminb_imput_status$x[,2], unique = T); colnames(exposures_status) <- paste0('s', colnames(exposures_status))

fitted_vals_status <- sim_LNM_FE(covariates=TMB_data_with_subset_for_nlminb_imput_status$x,
                                      TMB_res=res_nlminb_Fed2_imput_status)


ggplot(dcast(melt(list(fitted_vals_status,
                       exposures_status)), Var1+Var2~L1, value.var = "value"),
       aes(x=`1`,`2`, group=`1`, col=as.numeric(gsub("\\..*", "", gsub("^X", "", Var1)))))+
  geom_violin()+ geom_jitter(alpha=0.2)+ facet_wrap(.~Var2, scales = "free")+
  geom_abline(slope = 1, intercept = 0)+theme_bw()+
  labs(x='Fitted values', y='Observed values', col='Status')

##---------------------------------------------------------------------------------##
clinbrit$pfs
pfs_x = clinbrit$pfs[match(patient.meta$PATIENT_ID, paste0('BRITROC-', clinbrit$britroc_number))]
TMB_data_with_subset_for_nlminb_imput_pfs <- prepare_TMB_data_with_subset(exposures = impute(exposures, 1e-2), 1:7,
                                                                      .keep_additional=!(is.na(pfs_x)), add_x_vec = pfs_x)
TMB_data_with_subset_for_nlminb_imput_pfs$z <- diag(nrow(TMB_data_with_subset_for_nlminb_imput_pfs$z)) ## need to keep this for size of ularge
TMB_data_with_subset_for_nlminb_imput_pfs$num_individuals <- NULL

res_nlminb_Fed2_imput_pfs <- wrapper_run_TMB_use_nlminb(model = "tmb_MVN_partial_ILR_FEb",
                                              object = TMB_data_with_subset_for_nlminb_imput_pfs,
                                              use_nlminb = T)
plot_betas(res_nlminb_Fed2_imput_pfs)

fitted_vals_pfs <- melt(sim_LNM_FE(covariates=TMB_data_with_subset_for_nlminb_imput_pfs$x,
                          TMB_res=res_nlminb_Fed2_imput_pfs))
simulated_vals_pfs <- melt(sim_LNM_FE(covariates=cbind(1, seq(0, 1000, 5)),
                             TMB_res=res_nlminb_Fed2_imput_pfs))
# fitted_vals_pfs$Var1 <- as.numeric(gsub("\\..*", "", gsub("^X", "", fitted_vals_pfs$Var1)))
simulated_vals_pfs$Var1 <- as.numeric(gsub("\\..*", "", gsub("^X", "", simulated_vals_pfs$Var1)))
ggplot((fitted_vals_pfs), aes(x=Var1,
                              col=Var2, y=value))+geom_bar(stat = "identity")+
  geom_vline(xintercept = factor("X80"))+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+labs(y='Exposures', fill='') ##???


ggplot((simulated_vals_pfs), aes(x=Var1,
                                 col=Var2, y=value))+geom_bar(stat = "identity")#+
  geom_vline(xintercept = factor("X80"))

exposures_pfs <- as(compositions::ilrInv( TMB_data_with_subset_for_nlminb_imput_pfs$Y), 'matrix'); rownames(exposures_pfs) <- make.names(TMB_data_with_subset_for_nlminb_imput_pfs$x[,2], unique = T); colnames(exposures_pfs) <- paste0('s', colnames(exposures_pfs))

ggplot(dcast(melt(list(fitted_vals_pfs,
                       exposures_pfs)), Var1+Var2~L1, value.var = "value"),
       aes(x=`1`,`2`, group=`1`, col=as.numeric(gsub("\\..*", "", gsub("^X", "", Var1)))))+
  geom_point(alpha=0.2)+ facet_wrap(.~Var2, scales = "free")+
  geom_abline(slope = 1, intercept = 0)+theme_bw()+
  labs(x='Fitted values', y='Observed values', col='Status')

##-------------------------------------------------------------------------------------------##
## samples undergoing WGD from archival to relapse
undergoingWGD <- c('IM_124/JBLAB−4186', 'IM_124/JBLAB−4187', 'IM_124/JBLAB−4188',
  'IM_124/JBLAB−4189', 'IM_295/JBLAB−4960', 'IM_340/JBLAB−4965',
  'IM_341/JBLAB−4965', 'IM_383/JBLAB−19330', 'IM_395/JBLAB−19338',
  'IM_396/JBLAB−19338', 'IM_423/JBLAB−4996', 'IM_56/JBLAB−4128')

undergoingWGD <- unique(as.vector(sapply(undergoingWGD, function(i) strsplit(i, '/')[[1]])))
undergoingWGD
match(undergoingWGD, patient.meta$SAMPLE_ID)
sum(!is.na(match(undergoingWGD, patient.meta$SAMPLE_ID)))

remove_na <- function(i) i[!is.na(i)]
undergoingWGD_patients <- as.vector(unique(remove_na(patient.meta$PATIENT_ID[match(undergoingWGD, patient.meta$SAMPLE_ID)])))
undergoingWGD_patientsnum <- as.numeric(gsub("BRITROC-", "", undergoingWGD_patients))

head(clinbrit)
head(clinbrit2)

clinbrit[match(undergoingWGD_patientsnum, clinbrit$britroc_number),'Undergoing_WGD'] <- 'Arx/rel WGD'
clinbrit2[match(undergoingWGD_patientsnum, clinbrit2$britroc_number),'Undergoing_WGD'] <- 'Arx/rel WGD'

clinbrit$Undergoing_WGD[(clinbrit$Undergoing_WGD) != 'Arx/rel WGD'] <- 'No ploidy change'
clinbrit2$Undergoing_WGD[(clinbrit2$Undergoing_WGD) != 'Arx/rel WGD'] <- 'No ploidy change'

colnames(clinbrit)

plot(density(patient.meta$ploidy))

ggplot(clinbrit, aes(x=Undergoing_WGD, y=age))+geom_boxplot()+theme_bw()+geom_jitter()
ggplot(clinbrit, aes(x=Undergoing_WGD, y=status))+geom_boxplot()+theme_bw()+geom_jitter()
ggplot(clinbrit, aes(x=Undergoing_WGD, y=pt_sensitivity_at_reg))+geom_boxplot()+theme_bw()+geom_jitter()

plot(sort(patient.meta$ploidy))
threshold_WGD <- 2.6

patients_with_archival_WGD <- unique(as.numeric(gsub("BRITROC-", "",
   patient.meta$PATIENT_ID[patient.meta[patient.meta$group == 'arx','ploidy'] > threshold_WGD])))

##it looks like we don't have information for some patients in patient.meta (no sequencing I guess)
clinbrit$britroc_number[is.na(match(patients_with_archival_WGD, clinbrit$britroc_number))]
clinbrit2$britroc_number[is.na(match(patients_with_archival_WGD, clinbrit2$britroc_number))]

length(unique(clinbrit2$britroc_number))
length(unique(clinbrit$britroc_number))
length(unique(patient.meta$PATIENT_ID))

clinbrit[remove_na(match(patients_with_archival_WGD, clinbrit$britroc_number)),'Archival_WGD'] <- 'Archival WGD'
clinbrit2[remove_na(match(patients_with_archival_WGD, clinbrit2$britroc_number)),'Archival_WGD'] <- 'Archival WGD'
clinbrit$Archival_WGD[(clinbrit$Archival_WGD) != 'Archival WGD'] <- 'Archival diploid'
clinbrit2$Archival_WGD[(clinbrit2$Archival_WGD) != 'Archival WGD'] <- 'Archival diploid'


ggplot(clinbrit, aes(x=interaction(Undergoing_WGD,Archival_WGD), y=age))+geom_boxplot()+
  theme_bw()+geom_jitter()

ggplot(clinbrit, aes(x=interaction(Undergoing_WGD,Archival_WGD), y=pfs))+geom_boxplot()+
  theme_bw()+geom_jitter()


clades_WGD0 <- readRDS("/Users/morril01/Documents/PhD/other_repos/Vias_Brenton/copy_number_analysis_organoids/robjects/umap_ploidy_genes.RDS")

clades_WGD

clades_WGD <- clades_WGD0[remove_na(match(patient.meta$SAMPLE_ID, clades_WGD0$sample)),]
clades_WGD <- clades_WGD[patient.meta$group[match(clades_WGD$sample, patient.meta$SAMPLE_ID)] == 'arx',]

clades_WGD$patient <- gsub("BRITROC-", "", patient.meta$PATIENT_ID[match(clades_WGD$sample, patient.meta$SAMPLE_ID)])

patient.meta$WGD_UMAP <- NA
patient.meta$WGD_UMAP[match(rownames(clades_WGD), patient.meta$SAMPLE_ID)] <- (clades_WGD$`2` < 3)

ggplot(clades_WGD,
       aes(x=`1`, y=`2`, label=sample))+
  geom_point( alpha=0.4, aes( color=clade, size=MYC))+
  theme_bw()+
  geom_point(data = clades_WGD0, aes(x=`1`, y=`2`), alpha=0.05)

ggplot(clades_WGD,
       aes(x=`1`, y=`2`, label=sample, shape=dataset))+
  geom_point( alpha=0.4, aes( color=clade, size=MYC))+
  theme_bw()

clinbrit[remove_na(match(clades_WGD[clades_WGD$`2` > 3,]$patient, clinbrit$britroc_number)),'Archival_WGD_umap'] <- 'Archival diploid'
clinbrit[remove_na(match(clades_WGD[clades_WGD$`2` < 3,]$patient, clinbrit$britroc_number)),'Archival_WGD_umap'] <- 'Archival WGD'

# clinbrit$Archival_WGD_umap <- NA
# clinbrit$Archival_WGD_umap[remove_na(match(clinbrit$britroc_number, as.numeric(clades_WGD[clades_WGD$`2` > 3,]$patient)))] <- 'Diploid'
# clinbrit$Archival_WGD_umap[-remove_na(match(clinbrit$britroc_number, as.numeric(clades_WGD[clades_WGD$`2` > 3,]$patient)))] <- 'WGD'
clades_WGD[clades_WGD$`2` > 3,]$patient

table(ploidy=clinbrit$Archival_WGD, umap=clinbrit$Archival_WGD_umap)

ggplot(clinbrit, aes(x=interaction(Undergoing_WGD,Archival_WGD_umap), y=age))+geom_boxplot()+
  theme_bw()+geom_jitter()

#-----------------------------------------------------------------------------------------#

length(unique(yes_TP53[yes_TP53 %in% as.character(patient.meta$SAMPLE_ID)]))
length(unique(no_TP53[no_TP53 %in% as.character(patient.meta$SAMPLE_ID)]))

## has more time elapsed for the patients who undergo WGD than for those who don't?
clinbrit <- clinbrit[clinbrit$britroc_number %in% gsub("BRITROC-", "", patient.meta$PATIENT_ID),]
progression_interval <- as.Date(clinbrit$progression_date) - as.Date(clinbrit$diagnosis_date)
clinbrit$pfs
progression_interval ## same

#-----------------------------------------------------------------------------------------#
## add SVM classification of WGD
SVM_WGC <- readRDS("../../../../Vias_Brenton/copy_number_analysis_organoids/robjects/SVM_WGD.RDS")

all(rownames(exposures) == patient.meta$SAMPLE_ID)

samples_that_undergo_WGD <- c('BRITROC-23', 'BRITROC-74', 'BRITROC-209', 'BRITROC-216', 'BRITROC-241', 'BRITROC-267', 'BRITROC-274')

boxplot(list(clinbrit$pfs[match(gsub("BRITROC-", "", samples_that_undergo_WGD), clinbrit$britroc_number)],
     clinbrit$pfs[-match(gsub("BRITROC-", "", samples_that_undergo_WGD), clinbrit$britroc_number)]))

idx_undegoWGD <- match(gsub("BRITROC-", "", samples_that_undergo_WGD), clinbrit$britroc_number)
idx_not_undegoWGD <- -match(gsub("BRITROC-", "", samples_that_undergo_WGD), clinbrit$britroc_number)
idx_arxWGD <- -match(gsub("BRITROC-", "", samples_that_undergo_WGD), clinbrit$britroc_number)

patient.meta$PATIENT_ID <- as.character(patient.meta$PATIENT_ID)
patient.meta$SAMPLE_ID <- as.character(patient.meta$SAMPLE_ID)

## at least one archival, one relapse
multiple_samples_at_least_one <- patient.meta %>% group_by(PATIENT_ID) %>% summarise(tab=length(unique(group))) %>% filter(tab == 2)

patient.meta$WGD_progression <- NA
patient.meta$WGD_progression[!(patient.meta$PATIENT_ID %in% as.character(multiple_samples_at_least_one$PATIENT_ID))] <- 'Single group'


table(patient.meta$WGD_progression)
sum(is.na(patient.meta$WGD_progression))

as.character(multiple_samples_at_least_one$PATIENT_ID)

library(e1071)
patient.meta$WGD_SVM_both = predict(SVM_WGC$svmfit_both_rep[[1]]$svmfit, data.frame(as(compositions::clr(impute(exposures, 1e-2)), 'matrix')))
patient.meta$WGD_SVM_ICGC = predict(SVM_WGC$svmfit_ICGC_rep[[1]]$svmfit, data.frame(as(compositions::clr(impute(exposures, 1e-2)), 'matrix')))
patient.meta$WGD_SVM_TCGA = predict(SVM_WGC$svmfit_rep[[1]]$svmfit, data.frame(as(compositions::clr(impute(exposures, 1e-2)), 'matrix')))

boxplot(split(patient.meta$ploidy, patient.meta$WGD_SVM_both))
boxplot(split(patient.meta$ploidy, patient.meta$WGD_SVM_ICGC))
boxplot(split(patient.meta$ploidy, patient.meta$WGD_SVM_TCGA))

## several WGD samples seem to be categorised as non-WGD
table(TCGA=patient.meta$WGD_SVM_TCGA, ICGC=patient.meta$WGD_SVM_ICGC)

## alternative: using the ICGC SVM, which seems to be more conservative in calling non-WGD
multiple_samples_at_least_one$archival_is_WGD = sapply(multiple_samples_at_least_one$PATIENT_ID, function(i) paste0(unique(patient.meta$WGD_SVM_both[(patient.meta$PATIENT_ID == i) & (patient.meta$group == "arx")]  )))
multiple_samples_at_least_one$archival_is_WGD_ICGC = sapply(multiple_samples_at_least_one$PATIENT_ID, function(i) paste0(unique(patient.meta$WGD_SVM_ICGC[(patient.meta$PATIENT_ID == i) & (patient.meta$group == "arx")]  )))
multiple_samples_at_least_one$archival_is_WGD_TCGA = sapply(multiple_samples_at_least_one$PATIENT_ID, function(i) paste0(unique(patient.meta$WGD_SVM_TCGA[(patient.meta$PATIENT_ID == i) & (patient.meta$group == "arx")]  )))
table(multiple_samples_at_least_one$archival_is_WGD)
table(multiple_samples_at_least_one$archival_is_WGD_ICGC)
table(multiple_samples_at_least_one$archival_is_WGD_TCGA)

patient.meta$WGD_progression[!(patient.meta$PATIENT_ID %in% names(which(multiple_samples_at_least_one$archival_is_WGD == "FALSE")))] <- 'Diploid archival'
patient.meta$WGD_progression[!(patient.meta$PATIENT_ID %in% names(which(multiple_samples_at_least_one$archival_is_WGD == "TRUE")))] <- 'WGD archival'
patient.meta$WGD_progression[patient.meta$PATIENT_ID %in% samples_that_undergo_WGD] <- 'Undergoing WGD'
patient.meta$WGD_progression[patient.meta$PATIENT_ID == "BRITROC-67"] <- 'Revertal WGD'

table(patient.meta$WGD_progression)
sum(is.na(patient.meta$WGD_progression)) ## all have been classified

patient.meta_single_obs_per_patient <- patient.meta[!duplicated(patient.meta$PATIENT_ID),]
dim(patient.meta_single_obs_per_patient)
dim(patient.meta)

patient.meta_single_obs_per_patient <- cbind(patient.meta_single_obs_per_patient,
                                             clinbrit[match(gsub("BRITROC-", "", patient.meta_single_obs_per_patient$PATIENT_ID), clinbrit$britroc_number), c('pfs', 'status', 'os', 'age')])

library(survival)
library(survminer)
library(gg)
patient.meta_single_obs_per_patient[!is.na(patient.meta_single_obs_per_patient$status),'WGD_progression']
patient.meta_single_obs_per_patient[!is.na(patient.meta_single_obs_per_patient$status),'status']
patient.meta_single_obs_per_patient[!is.na(patient.meta_single_obs_per_patient$status),'os']
patient.meta_single_obs_per_patient$WGD_progression <- factor(patient.meta_single_obs_per_patient$WGD_progression)
patient.meta_single_obs_per_patient <- patient.meta_single_obs_per_patient[!is.na(patient.meta_single_obs_per_patient$status),]

# resCox= coxph(Surv(status, os) ~ WGD_progression,
resCox= coxph(Surv(status) ~ WGD_progression+age,
              data =  patient.meta_single_obs_per_patient)
resCox
pCox1 = ggforest(resCox, data = patient.meta_single_obs_per_patient,
                              main = NULL)
pCox1
## undegoing WGD is ss

resCox= coxph(Surv(os) ~ WGD_progression+age,
              data =  patient.meta_single_obs_per_patient)
resCox
pCox1 = ggforest(resCox, data = patient.meta_single_obs_per_patient,
                 main = NULL)
pCox1

resCox= coxph(Surv(os, status) ~ WGD_progression+age,
              data =  patient.meta_single_obs_per_patient)
resCox
pCox1 = ggforest(resCox, data = patient.meta_single_obs_per_patient,
                 main = NULL)
pCox1
## undegoing WGD is ss

resCox= coxph(Surv(pfs) ~ WGD_progression+age,
              data =  patient.meta_single_obs_per_patient)
resCox
pCox1 = ggforest(resCox, data = patient.meta_single_obs_per_patient,
                 main = NULL)
pCox1

resCox= coxph(Surv(pfs) ~ WGD_progression+age,
              data =  patient.meta_single_obs_per_patient)
resCox
pCox1 = ggforest(resCox, data = patient.meta_single_obs_per_patient,
                 main = NULL)
pCox1

kmFit_groups= survfit(Surv(os, status) ~ WGD_progression, data = patient.meta_single_obs_per_patient)
ggsurvplot( kmFit_groups, data = patient.meta_single_obs_per_patient,
                          risk.table = TRUE, pval = TRUE, conf.int = FALSE,
                          xlim = c(0,2000), break.time.by = 365, 
                          surv.median.line = "hv", risk.table.y.text.col = TRUE, risk.table.y.text = FALSE,
                          ggtheme = theme_bw() ) + 
  ggtitle("Survival in HGSOC by s3/s4 clusterm, in high s3 group")

#-----------------------------------------------------------------------------------------#

#-----------------------------------------------------------------------------------------#
table(SVM=patient.meta$WGD_SVM_both, UMAP=patient.meta$WGD_UMAP) ## several cases where SVM says they are diploid and UMAP that they are WGD

t.test(clinbrit$pfs[idx_undegoWGD],
       clinbrit$pfs[idx_not_undegoWGD])
mean(clinbrit$pfs[idx_undegoWGD])
sd(clinbrit$pfs[idx_undegoWGD])

mean(clinbrit$pfs[idx_not_undegoWGD], na.rm = T)
sd(clinbrit$pfs[idx_not_undegoWGD], na.rm = T)


