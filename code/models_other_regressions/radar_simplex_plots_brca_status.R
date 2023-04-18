#-------------------------------------------------------------------------------------------#
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
#-------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------#
folder_out_RDS <- "../../data/inference/"
source("../models_archival_relapse/helper/functions.R")
source("../models_archival_relapse/helper/header.R")
source("../models_archival_relapse/helper/helper_TMB.R")
library(fmsb) ## radar plots
library(Ternary) ## ternary plots
library(ggradar) ## Added by PS
library(dplyr)
library(ggplot2)
folder_images_out <- "../../../results/partialILRmodelling_ME_brca/"
#-------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------#
## created in LN_modelling_brca_resistance_status.R
ME_brcastatus_partialILRnocor_partialILRnocor = readRDS(
  paste0(folder_out_RDS,"ME_brcastatus_partialILRnocor_partialILRnocor.RDS"))
#-------------------------------------------------------------------------------------------#
## Added for paper figures - PS
path="../../data/britroc_30kb_signature_data.rds"
load(path)
exposures=t(sig_quants)
source("../../../britroc-cn-analysis/colour_palettes.R")


#-------------------------------------------------------------------------------------------#
## radar/spider plot

archival_fit <- python_like_select_name(ME_brcastatus_partialILRnocor_partialILRnocor$par.fixed, 'beta')[c(T,F)]
slope_fit <- python_like_select_name(ME_brcastatus_partialILRnocor_partialILRnocor$par.fixed, 'beta')[c(F,T)]

all(rownames(exposures) %in% patient.meta$SAMPLE_ID)

brca_status <- read.table("../../data/britroc_brca_patients_list.tsv",header = T,sep = "\t")
patient.meta$brca_status <- factor(ifelse(patient.meta$PATIENT_ID %in% brca_status$PATIENT_ID,"brca","non-brca"),
                                         levels = c("brca","non-brca"))

## ggplot radar plots and additional visualisation - PS
brca_data <- as.data.frame(rbind(exposures[patient.meta$brca_status == 'brca',]), ncol=7) %>%
              tibble::rownames_to_column(var = "id") %>%
              tidyr::pivot_longer(cols = 2:8,names_to = "signature") %>%
              mutate(group = rep("BRCA",nrow(.)))

nonbrca_data <- as.data.frame(rbind(exposures[patient.meta$brca_status == 'non-brca',]), ncol=7) %>%
  tibble::rownames_to_column(var = "id") %>%
  tidyr::pivot_longer(cols = 2:8,names_to = "signature") %>%
  mutate(group = rep("non-BRCA",nrow(.)))

arx_beta=as.numeric(compositions::ilrInv(archival_fit))
names(arx_beta) <- paste0("s",1:7)

rlps_beta <- as.numeric(compositions::ilrInv(slope_fit+archival_fit))
names(rlps_beta) <- paste0("s",1:7)

beta_data <- as.data.frame(rbind(arx_beta,rlps_beta),row.names = c("BRCA","non-BRCA")) %>%
              tibble::rownames_to_column(var = "group") %>%
              tidyr::pivot_longer(cols = 2:8,names_to = "signature")
beta_data$group <- factor(beta_data$group,levels = c("BRCA","non-BRCA"))

combined_radar <- rbind(brca_data,nonbrca_data)
combined_radar$group <- factor(combined_radar$group,levels = c("BRCA","non-BRCA"))

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

brca_radar_beta <- ggplot() +
  geom_polygon(data = combined_radar,
               aes(signature,value,fill=group,group=id),
               alpha=0.01) +
  geom_polygon(data = beta_data,
             aes(signature,value,colour=group,group=group),size=1,fill=NA) +
  geom_point(data = beta_data,
               aes(signature,value,fill=group),colour="gray15",pch=21,size=4) +
  scale_y_continuous(limits = c(0,max(beta_data$value)),breaks = seq.int(0,1,0.1)) +
  scale_fill_manual(values = colour_palettes$brca_carriers) +
  scale_colour_manual(values = colour_palettes$brca_carriers) +
  coord_radar() +
  theme_bw() +
  theme(legend.position = "bottom") +
  theme(#axis.line = element_blank(),
    axis.text.y = element_blank(),
    panel.grid = element_line(color="grey80"),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    plot.margin = margin(0,0,0,0))

brca_radar_beta

brca_bar_beta <- ggplot() +
  geom_point(data = combined_radar,
              aes(signature,value,color=group,group=group),position = position_jitterdodge(),
              alpha=0.2) +
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
  scale_fill_manual(values = colour_palettes$brca_carriers) +
  scale_colour_manual(values = colour_palettes$brca_carriers) +
  theme_bw() +
  theme(legend.position = "bottom")
brca_bar_beta

saveRDS(brca_radar_beta,file = "../../../britroc-cn-analysis/plots/brca_status_radar_partialILR_beta.RDS")
ggsave(filename = "../../../britroc-cn-analysis/plots/brca_status_radar_partialILR_beta.png",
       plot = brca_radar_beta,device = "png",width = 8,height = 8,units = "in",dpi = 300)
ggsave(filename = "../../../britroc-cn-analysis/plots/brca_status_radar_partialILR_beta.pdf",
       plot = brca_radar_beta,device = "pdf",width = 8,height = 8,units = "in",dpi = 300)

saveRDS(brca_bar_beta,file = "../../../britroc-cn-analysis/plots/brca_status_bar_partialILR_beta.RDS")
ggsave(filename = "../../../britroc-cn-analysis/plots/brca_status_bar_partialILR_beta.png",
       plot = brca_bar_beta,device = "png",width = 10,height = 8,units = "in",dpi = 300)
ggsave(filename = "../../../britroc-cn-analysis/plots/brca_status_bar_partialILR_beta.pdf",
       plot = brca_bar_beta,device = "pdf",width = 10,height = 8,units = "in",dpi = 300)
#-------------------------------------------------------------------------------------------#
