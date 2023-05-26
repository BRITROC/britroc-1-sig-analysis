#-------------------------------------------------------------------------------------------#
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
#-------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------#
folder_out_RDS <- "../../../out/inference/partialILR/"
source("../helper/functions.R")
source("../helper/header.R")
source("../helper/helper_functions.R")
source("../helper/helper_TMB.R")
source("../../../../britroc-1-cn-analysis/colour_palettes.R")
library(fmsb) ## radar plots
library(Ternary) ## ternary plots
library(ggradar) ## Added by PS
library(dplyr) ## Added by PS
folder_images_out <- "../../../results/partialILRmodelling_ME/"
#-------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------#
# ext func

#-------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------#
## created in LN_modelling_zeros_partialILR_TMB.R
res_nlminb_nocoroutsidesd_only_matched_allsigs = readRDS(paste0(folder_out_RDS, "res_nlminb_nocoroutsidesd_only_matched_allsigs.RDS"))
#-------------------------------------------------------------------------------------------#
## Added for paper figures - PS
path="../../../data/britroc_30kb_signature_data.rds"
load(path)
exposures=t(sig_quants)
source("../../../../britroc-1-cn-analysis/colour_palettes.R")
#-------------------------------------------------------------------------------------------#
## radar/spider plot

archival_fit <- python_like_select_name(res_nlminb_nocoroutsidesd_only_matched_allsigs$par.fixed, 'beta')[c(T,F)]
slope_fit <- python_like_select_name(res_nlminb_nocoroutsidesd_only_matched_allsigs$par.fixed, 'beta')[c(F,T)]

all(rownames(exposures) == patient.meta$SAMPLE_ID)

# pdf("../../../results/partialILRmodelling_ME/radarplot_with_res_nlminb_nocoroutsidesd_only_matched_allsigs_fit.pdf", width = 6, height = 2.5)
tikzDevice::tikz("../../../results/partialILRmodelling_ME/radarplot_with_res_nlminb_nocoroutsidesd_only_matched_allsigs_fit.tex", width = 6, height = 2.5)
#pdf("../../../../BriTROC/plots/radarplot_with_res_nlminb_nocoroutsidesd_only_matched_allsigs_fit.pdf")
par(mfrow=c(1,2), mar=c(0,0,2,0))

radarchart(as.data.frame(rbind(1, 0, exposures[patient.meta$group == 'arx',], as.numeric(compositions::ilrInv(archival_fit))), ncol=7),
           cglty = 1,
           plty = 1,
           pfcol = rgb(0, 0.4, 1, 0.01),
           pcol=c(rep(NA, sum(patient.meta$group == "arx")), 1),
           title = "Archival")

radarchart(as.data.frame(rbind(1, 0, exposures[patient.meta$group == 'rlps',], as.numeric(compositions::ilrInv(slope_fit+archival_fit))), ncol=7),
           cglty = 1,
           plty = 1,
           pfcol = rgb(0, 0.4, 1, 0.01),
           pcol=c(rep(NA, sum(patient.meta$group == "rlps")), 1),
           title = "Relapse")
dev.off()

tikzDevice::tikz("../../../results/partialILRmodelling_ME/radarplot_with_res_nlminb_nocoroutsidesd_only_matched_allsigs_fit2.tex", width = 6, height = 2.5)
par(mfrow=c(1,2), mar=c(0,0,2,0))
radarchart(as.data.frame(rbind(1, 0, exposures[patient.meta$group == 'arx',], as.numeric(compositions::ilrInv(archival_fit))), ncol=7), cglty = 1,
           plty = 1, pfcol = NA, pty=c(rep(NA, sum(patient.meta$group == "arx")), 19),
           lty=c(rep(0.1, sum(patient.meta$group == "arx")), 1),
           pcol=c(rep(rgb(0, 0.4, 1, 0.1), sum(patient.meta$group == "arx")), 1), title = "Diagnosis")
radarchart(as.data.frame(rbind(1, 0, exposures[patient.meta$group == 'rlps',], as.numeric(compositions::ilrInv(slope_fit+archival_fit))), ncol=7), cglty = 1,
           plty = 1, pfcol = NA, pty=c(rep(NA, sum(patient.meta$group == "rlps")), 19),
           lty=c(rep(0.1, sum(patient.meta$group == "rlps")), 1),
           pcol=c(rep(rgb(0, 0.4, 1, 0.1), sum(patient.meta$group == "rlps")), 1), title = "Relapse")
dev.off()


## ggplot radar plots and additional visualisation - PS
arx_data <- as.data.frame(rbind(exposures[patient.meta$group == 'arx',]), ncol=7) %>%
              tibble::rownames_to_column(var = "id") %>%
              tidyr::pivot_longer(cols = 2:8,names_to = "signature") %>%
              mutate(group = rep("diagnosis",nrow(.)))

rlps_data <- as.data.frame(rbind(exposures[patient.meta$group == 'rlps',]), ncol=7) %>%
  tibble::rownames_to_column(var = "id") %>%
  tidyr::pivot_longer(cols = 2:8,names_to = "signature") %>%
  mutate(group = rep("relapse",nrow(.)))

arx_beta=as.numeric(compositions::ilrInv(archival_fit))
names(arx_beta) <- paste0("s",1:7)

rlps_beta <- as.numeric(compositions::ilrInv(slope_fit+archival_fit))
names(rlps_beta) <- paste0("s",1:7)

beta_data <- as.data.frame(rbind(arx_beta,rlps_beta),row.names = c("diagnosis","relapse")) %>%
              tibble::rownames_to_column(var = "group") %>%
              tidyr::pivot_longer(cols = 2:8,names_to = "signature")

combined_radar <- rbind(arx_data,rlps_data)

## coord_radar function from
## https://stackoverflow.com/questions/62462681/ggradar-highlight-top-values-in-radar
coord_radar <- function (theta = "x", start = - pi / 2, direction = 1) {
  theta <- match.arg(theta, c("x", "y"))
  r <- if (theta == "x") "y" else "x"
  ggproto(
    "CordRadar", CoordPolar, theta = theta, r = r, start = start,
    direction = sign(direction),
    is_linear = function(coord) TRUE)
}

arx_relapse_radar_beta <- ggplot() +
  geom_polygon(data = combined_radar,
               aes(signature,value,fill=group,group=id),
               alpha=0.03) +
  geom_polygon(data = beta_data,
             aes(signature,value,colour=group,group=group),size=1,fill=NA) +
  geom_point(data = beta_data,
               aes(signature,value,fill=group),colour="gray15",pch=21,size=4) +
  scale_y_continuous(limits = c(0,max(beta_data$value)),breaks = seq.int(0,1,0.1)) +
  scale_fill_manual(values = colour_palettes$diagnosis_relapse) +
  scale_colour_manual(values = colour_palettes$diagnosis_relapse) +
  coord_radar() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(#axis.line = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_line(color="grey80"),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        plot.margin = margin(0,0,0,0))

arx_relapse_bar_beta <- ggplot() +
  geom_boxplot(data = combined_radar,
               aes(signature,value,fill=group),
               alpha=0.2) +
  # geom_line(data = beta_data,
  #              aes(signature,value,colour=group,group=group),fill=NA) +
  geom_point(data = beta_data,
             aes(signature,value,fill=group),
             colour="gray15",alpha=0.5,
             size=4,pch=21) +
  scale_y_continuous(limits = c(0,1)) +
  scale_fill_manual(values = colour_palettes$diagnosis_relapse) +
  scale_colour_manual(values = colour_palettes$diagnosis_relapse) +
  theme_bw() +
  theme(legend.position = "none")
             
# Export source data
write.table(x = combined_radar,file = "../../../../britroc-1-cn-analysis/source_data/figure_5C_layer_1.tsv",
            quote = F,sep = "\t",row.names = F,col.names = T)

write.table(x = beta_data,file = "../../../../britroc-1-cn-analysis/source_data/figure_5C_layer_2.tsv",
            quote = F,sep = "\t",row.names = F,col.names = T)
             
saveRDS(arx_relapse_radar_beta,file = "../../../../britroc-1-cn-analysis/plots/arx_relapse_radar_partialILR_beta.RDS")
ggsave(filename = "../../../../britroc-1-cn-analysis/plots/arx_relapse_radar_partialILR_beta.png",
       plot = arx_relapse_radar_beta,device = "png",width = 8,height = 8,units = "in",dpi = 300)
ggsave(filename = "../../../../britroc-1-cn-analysis/plots/arx_relapse_radar_partialILR_beta.pdf",
       plot = arx_relapse_radar_beta,device = "pdf",width = 8,height = 8,units = "in",dpi = 300)

saveRDS(arx_relapse_bar_beta,file = "../../../../britroc-1-cn-analysis/plots/arx_relapse_bar_partialILR_beta.RDS")
ggsave(filename = "../../../../britroc-1-cn-analysis/plots/arx_relapse_bar_partialILR_beta.png",
       plot = arx_relapse_bar_beta,device = "png",width = 10,height = 8,units = "in",dpi = 300)
ggsave(filename = "../../../../britroc-1-cn-analysis/plots/arx_relapse_bar_partialILR_beta.pdf",
       plot = arx_relapse_bar_beta,device = "pdf",width = 10,height = 8,units = "in",dpi = 300)

## all combinations of differences
# df_paired_differences <- give_LRchanges_barplot(obj_it = res_nlminb_nocoroutsidesd_only_matched_allsigs_TMBdata,
#                                                 return_df = T, nrow_facets = 3)

# df_paired_differences
both_arx_rlps <- sapply(unique(patient.meta$PATIENT_ID), function(i) all(c('rlps', 'arx') %in% 
                                                                           patient.meta$group[patient.meta$PATIENT_ID == i]))
both_arx_rlps <- rownames(exposures) %in% patient.meta$SAMPLE_ID[patient.meta$PATIENT_ID %in% patient.meta$PATIENT_ID[both_arx_rlps]]
both_arx_rlps <- as.character(patient.meta$PATIENT_ID[both_arx_rlps])
both_arx_rlps

paired_subtractions <- give_paired_subtractions(patients_list = both_arx_rlps, exposures = exposures)

## scaled beta slope coefficients
tikzDevice::tikz("../../../results/partialILRmodelling_ME/radarplot_with_res_nlminb_nocoroutsidesd_only_matched_difference.tex",
                 width = 2.5, height = 2.5)
par(mfrow=c(1,1), mar=c(0,0,2,0))
radarchart(add_colnames(as.data.frame(rbind(1, 0, as.numeric(compositions::ilrInv(slope_fit))/max(as.numeric(compositions::ilrInv(slope_fit)))), ncol=7),
                        colnames_arg = paste0('s', 1:7)), cglty = 3,
           plty = 1, pfcol = NA, pty=c(rep(NA, sum(patient.meta$group == "rlps")), 19),
           lty=c(rep(0.1, sum(patient.meta$group == "rlps")), 1),
           pcol=c(rep(rgb(0, 0.4, 1, 1.0), sum(patient.meta$group == "rlps")), 1), title = "Difference")
dev.off()

## this subtraction is NOT compositional!
tikzDevice::tikz("../../../results/partialILRmodelling_ME/radarplot_with_res_nlminb_nocoroutsidesd_only_matched_difference_subtractions.tex",
                 width = 2.5, height = 2.5)
par(mfrow=c(1,1), mar=c(0,0,2,0))
radarchart(as.data.frame(rbind(1, 0, paired_subtractions, as.numeric(compositions::ilrInv(slope_fit))), ncol=7), cglty = 1,
           plty = 1, pfcol = NA, pty=c(rep(NA, nrow(paired_subtractions)), 19),
           lty=c(rep(0.1, nrow(paired_subtractions)), 1),
           pcol=c(rep(rgb(0, 0.4, 1, 0.1), nrow(paired_subtractions)), 1), title = "Difference")
dev.off()
#-------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------#
## ternary plot with connections between matching archival and relapse
exposures <- exposures[rownames(exposures) %in% patient.meta$SAMPLE_ID,]
stopifnot(rownames(exposures) %in% patient.meta$SAMPLE_ID)

levels(factor(patient.meta$group))
pdf("../../../results/exploratory/simplex_arx.pdf", width = 3, height = 3)
give_ternary_v2(remove_all_NA(normalise_rw(exposures[patient.meta$group == "arx",c(3,4,5)])), legend_off = T) ## arx
dev.off()
pdf("../../../results/exploratory/simplex_rlps.pdf", width = 3, height = 3)
give_ternary_v2(remove_all_NA(normalise_rw(exposures[patient.meta$group == "rlps",c(3,4,5)])), legend_off = T) ## rlps
dev.off()

d = ncol(exposures)
TMB_data_all_samples = list(Y = exposures,
                            num_individuals = ncol(give_z_matrix_from_labels(patient.meta$PATIENT_ID)),
                            d = d,
                            n = nrow(exposures),
                            x = cbind(1, as.numeric(as.factor(patient.meta$group))-1),
                            z = give_z_matrix_from_labels(patient.meta$PATIENT_ID))

pdf(paste0(folder_images_out, "ternary_all.pdf"), width = 3, height = 3)
par(mfrow=c(1,1))
give_ternary_v3(TMB_data_all_samples, exposures, legend_off = T,col = colour_palettes$diagnosis_relapse) ## arx
dev.off()

combinations_k3 <- combn(7, 3)
apply(combinations_k3, 2, function(combination_sigs){
  #pdf(paste0(folder_images_out, "ternary_all_combinations/ternary_all_",
             #paste0(combination_sigs, collapse = '+'), ".pdf"),
      #width = 3, height = 3)
  pdf(paste0("../../../../britroc-1-cn-analysis/copy_number_signatures/plots/ternary_all_combinations/ternary_all_",
           paste0(combination_sigs, collapse = '+'), ".pdf"),
    width = 3, height = 3)
  par(mfrow=c(1,1))
  give_ternary_v3(TMB_data_all_samples, exposures, legend_off = T, selected_sigs = combination_sigs,col = colour_palettes$diagnosis_relapse) ## arx
  dev.off()
})
system("open ../../../results/partialILRmodelling_ME/ternary_all_combinations/")
#-------------------------------------------------------------------------------------------#



#-------------------------------------------------------------------------------------------#
