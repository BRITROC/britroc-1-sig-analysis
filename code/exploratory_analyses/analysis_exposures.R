## Infence using stan_fit_LNM for data in the simplex
## Simple LN, no fixed or random effects

rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(ggplot2)
library(reshape2)
library(pheatmap)
library(gridExtra)
library(dplyr)
source("../models/helper/functions.R")
source("../../../../GlobalDA/code/2_inference_TMB/helper_TMB.R")

load("../../data/britroc_30kb_signature_data.rds")
# patient.meta
# sample_by_component
# sig_quants

max(table(patient.meta$SAMPLE_ID))
max(table(patient.meta$PATIENT_ID))

exposures = t(sig_quants)
patient.meta = patient.meta[match(rownames(exposures), as.character(patient.meta$SAMPLE_ID)),]
sample_by_component = sample_by_component[match(rownames(exposures), rownames(sample_by_component)),]

all(rownames(exposures) == patient.meta$SAMPLE_ID)
all(rownames(exposures) == rownames(sample_by_component))

table(patient.meta$use)
table(patient.meta$paired)
pdf("../../results/exploratory/BriTROC2_pheatmap_exposures.pdf", height = 6, width = 6)
print(pheatmap(exposures,
         annotation_row = data.frame(group=patient.meta[,c('group')], row.names = as.character(patient.meta$SAMPLE_ID)),
         show_rownames = FALSE))
dev.off()

exposures_df <- (melt(exposures))
exposures_df_nos5 <- (melt(normalise_rw(exposures[,-5])))
exposures_df$group <- patient.meta$group[match(exposures_df$Var1, patient.meta$SAMPLE_ID)]
exposures_df$patient <- patient.meta$PATIENT_ID[match(exposures_df$Var1, patient.meta$SAMPLE_ID)]
exposures_df_nos5$group <- patient.meta$group[match(exposures_df_nos5$Var1, patient.meta$SAMPLE_ID)]
exposures_df_nos5$patient <- patient.meta$PATIENT_ID[match(exposures_df_nos5$Var1, patient.meta$SAMPLE_ID)]

exposures_df$group[exposures_df$group == 'arx'] <- 'Archival'
exposures_df$group[exposures_df$group == 'rlps'] <- 'Relapse'
exposures_df_nos5$group[exposures_df_nos5$group == 'arx'] <- 'Archival'
exposures_df_nos5$group[exposures_df_nos5$group == 'rlps'] <- 'Relapse'

ggplot(exposures_df, aes(x=Var1, y=value, fill=Var2))+geom_bar(stat = "identity")+facet_wrap(.~group, scales = "free_x")+
  scale_fill_brewer(palette="Dark2")+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+labs(y='Exposures', fill='')
ggsave("../../results/exploratory/BriTROC2_pheatmap_exposures_barplot.pdf", height = 3, width = 6)
ggplot(exposures_df, aes(x=Var2, y=value, col=group, group=interaction(Var2,group)))+
  geom_boxplot()+ geom_point(position=position_jitterdodge(),
                             alpha=0.2, size=0.9, width = 0.2)+
  theme_bw()+labs(x='Signature', y='Exposure', col='Group')
ggsave("../../results/exploratory/BriTROC2_pheatmap_exposures_boxplot.pdf", height = 3, width = 6)

ggplot(exposures_df_nos5, aes(x=Var2, y=value, col=group, group=interaction(Var2,group)))+
  geom_boxplot()+ geom_point(position=position_jitterdodge(),
                             alpha=0.2, size=0.9, width = 0.2)+
  theme_bw()+labs(x='Signature', y='Exposure', col='Group')
ggsave("../../results/exploratory/BriTROC2_pheatmap_exposures_boxplot_nos5.pdf", height = 3, width = 6)


exposures_df$group = factor(exposures_df$group, levels = unique(exposures_df$group))

give_changes_with_direction <- function(exposures_df){
  lapply(unique(exposures_df$patient), function(patient_it){
  cur_exposures <- exposures_df[exposures_df$patient == patient_it,]
  ## for each sample
  cur_samples <- unique(cur_exposures$Var1)
  cur_groups <- cur_exposures$group[sapply(cur_samples, function(i) which(cur_exposures$Var1 == i)[1] )]
  
  if(length(unique(cur_groups)) == 2){
    ## at least one archival and one relapse
    ## for all combinations of archival and relapse, compute the direction of change
    arc_patient <- cur_samples[cur_groups == 'Archival']
    rel_patient <- cur_samples[cur_groups == 'Relapse']
    return(do.call('rbind', lapply(arc_patient, function(arc_patient_it){
      do.call('rbind', lapply(rel_patient, function(rel_patient_it){
        exp_rel <- cur_exposures[cur_exposures$Var1 == rel_patient_it,]
        exp_arx <- cur_exposures[cur_exposures$Var1 == arc_patient_it,]
        .a <- rbind.data.frame(cbind.data.frame(exp_rel,
                                          direction_change=(exp_rel['value'] -exp_arx['value'])),
                         cbind.data.frame(exp_arx,
                                          direction_change=(exp_rel['value'] -exp_arx['value'])))
        colnames(.a)[ncol(.a)] <- 'direction_change'
        .a
        
      }))
    })))

  }else{
    return(cbind.data.frame(cur_exposures, direction_change=NA))
  }
  
})
}
exposures_df_with_changes <- give_changes_with_direction(exposures_df)
exposures_df_with_changes <- do.call('rbind', exposures_df_with_changes)
exposures_df_with_changes$direction_change_bool <- exposures_df_with_changes$direction_change > 0
exposures_df_with_changes_nos5 <- give_changes_with_direction(exposures_df_nos5)
exposures_df_with_changes_nos5 <- do.call('rbind', exposures_df_with_changes_nos5)
exposures_df_with_changes_nos5$direction_change_bool <- exposures_df_with_changes_nos5$direction_change > 0

dim(exposures_df_with_changes)
dim(exposures_df)

# exposures_df_paired <- exposures_df[exposures_df$patient %in% names(table(exposures_df$patient)[table(exposures_df$patient) == 14]),]
# unique(exposures_df_paired$patient)
# exposures_df_paired_list <- list(exposures_df_paired[exposures_df_paired$group == 'Archival',],
#                             exposures_df_paired[exposures_df_paired$group == 'Relapse',])
# exposures_df_paired_list <- lapply(exposures_df_paired_list, function(j) j[with(exposures_df_paired_list[[1]], order(Var2, patient)), ])
# 
# direction_change <- exposures_df_paired_list[[2]]$value - exposures_df_paired_list[[1]]$value
# names(direction_change) <- exposures_df_paired_list[[1]]$patient
# lapply(unique(exposures_df_paired$patient), function(i) exposures_df_paired[(exposures_df_paired$patient == i) & 
#                                                                               (exposures_df_paired$))
# exposures_df$direction_change <-  (direction_change[match(exposures_df$patient, names(direction_change))] > 0 )

ggplot(exposures_df_with_changes, aes(x=interaction(Var2, group), y=value, group=patient))+
  geom_point(size=0.5)+geom_line(aes(col=direction_change_bool))+
  theme_bw()+labs(x='Signature', y='Exposure', col='Direction of')+facet_wrap(.~Var2, scales = "free_x", nrow=1)+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave("../../results/exploratory/BriTROC2_pheatmap_exposures_lines.pdf", height = 3, width = 10)

ggplot(exposures_df_with_changes[!is.na(exposures_df_with_changes$direction_change_bool),], aes(x=(Var2), fill=direction_change_bool))+
  geom_bar()+
  theme_bw()+labs(x='Signature', y='Exposure', col='Direction of')+facet_wrap(.~Var2, scales = "free_x", nrow=1)+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave("../../results/exploratory/BriTROC2_pheatmap_exposures_directionofchangebarplot.pdf", height = 3, width = 6)

ggplot(exposures_df_with_changes_nos5[!is.na(exposures_df_with_changes_nos5$direction_change_bool),], aes(x=(Var2), fill=direction_change_bool))+
  geom_bar()+
  theme_bw()+labs(x='Signature', y='Exposure', col='Direction of')+facet_wrap(.~Var2, scales = "free_x", nrow=1)+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave("../../results/exploratory/BriTROC2_pheatmap_exposures_directionofchangebarplot_nos5.pdf", height = 3, width = 6)

ggplot(exposures_df_with_changes_nos5, aes(x=interaction(Var2, group), y=value, group=patient))+
  geom_point(size=0.5)+geom_line(aes(col=direction_change_bool))+
  theme_bw()+labs(x='Signature', y='Exposure', col='Direction of')+facet_wrap(.~Var2, scales = "free_x", nrow=1)+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave("../../results/exploratory/BriTROC2_pheatmap_exposures_lines_nos5.pdf", height = 3, width = 10)


exposures_zeros = apply(exposures > 0, 2, as.numeric)
pdf("../results/exploratory/BriTROC2_pheatmap_exposures_zeros.pdf", height = 6, width = 6)
print(pheatmap(data.frame(exposures_zeros, row.names=rownames(exposures)),
               annotation_row = data.frame(group=patient.meta[,c('group')], row.names = as.character(patient.meta$SAMPLE_ID)),
               show_rownames = FALSE))
dev.off()

clean_inf_na = function(mat_A){
  mat_A[is.infinite(mat_A)] = NA
  mat_A
}

change_colnames = function(df, new_colnames){
  colnames(df) = new_colnames
  df
}

pdf("../results/exploratory/BriTROC2_pheatmap_exposures_logratios.pdf", height = 6, width = 6)
print(pheatmap(data.frame(change_colnames(clean_inf_na(log(exposures/exposures[,7])[,1:6]), paste0('LR', 1:6)), row.names=rownames(exposures)),
               annotation_row = data.frame(group=patient.meta[,c('group')], row.names = as.character(patient.meta$SAMPLE_ID)),
               show_rownames = FALSE, cluster_cols = FALSE, cluster_rows = FALSE, cellheight=1.2,cellwidth=40))
dev.off()

pca = prcomp(exposures, scale. = TRUE)
pca$sdev
eigs <- pca$sdev^2
eigs[1] / sum(eigs)
df = cbind.data.frame(pca$x[,1:4],
                      group=patient.meta$group, patient=patient.meta$PATIENT_ID,
                      exposures)
df$patient2 = df$patient

## basic PCA
ggplot()+
  geom_point(data=df[,c("PC1","PC2","group","patient")], aes(x=PC1, y=PC2, alpha=0.99))+
  scale_alpha(range = c(0.1, .2))

## PCA faceting by patient, and showing all exposures as background (takes a bit long)
# ggplot()+
#   geom_point(data=df[,c("PC1","PC2","group","patient")], aes(x=PC1, y=PC2, alpha=0.8))+
#   scale_alpha(range = c(0.1, .2))+
#   geom_point(data=df, aes(x=PC1, y=PC2, col='red'))+
#   guides(col=FALSE, alpha=FALSE)+
#   facet_wrap(.~patient2)
# ggsave("../results/exploratory/PCA_with_background.png", height = 12, width = 12)

## PCA with colouring of groups, and joining samples from the same patient with a line
eigs <- pca$sdev^2
eigs[1] / sum(eigs)
ggplot(df,
       aes(x=PC1, y=PC2))+geom_point(aes(col=group))+geom_line(aes(group=patient))+
  labs(x=paste0('PC1 (', round(eigs[1] / sum(eigs)*100,2), '% variance explained)'),
       y=paste0('PC2 (', round(eigs[2] / sum(eigs)*100,2), '% variance explained)'))
# ggsave("../results/exploratory/PCA_with_links.png", height = 5)

## how many samples do we have for each patient?
## options:
#' a) a single or multiple arx sample, no rlps
#' b) a single or multiple rlps sample, no arx
#' c) a single pair of arx-rlps
#' d1) a pair of arx and rlps, and more arx
#' d2) a pair of arx and rlps, and more rlps
classification_patients = sapply(unique(df$patient), function(i){
  if('arx' %in% df$group[df$patient == i]){
    ## a, c, or d
    if('rlps' %in% df$group[df$patient == i]){
      if(length(df$group[df$patient == i]) == 2){
        return('A single pair of arx-rlps')
      }else{
        if(sum(df$group[df$patient == i] == 'rlps') > 1){
          return('Pair with multiple rlps')
        }else{
          return('Pair with multiple arx')
        }
      }
    }else{
      return('Only arx')
    }
  }else{
    return('Only rlps')
  }
})
names(classification_patients) = unique(df$patient)

xtable::xtable(table(classification_patients))

df_only_paired = df[df$patient %in% names(classification_patients)[classification_patients %in% c('A single pair of arx-rlps', 'Pair with multiple arx',  'Pair with multiple rlps')],]
pca_only_paired = pca
pca_only_paired$x = pca_only_paired$x[df$patient %in% names(classification_patients)[classification_patients %in% c('A single pair of arx-rlps', 'Pair with multiple arx',  'Pair with multiple rlps')],]

ggplot(df_only_paired,
       aes(x=PC1, y=PC2))+geom_point(aes(col=group))+geom_line(aes(group=patient))+
  labs(x=paste0('PC1 (', round(eigs[1] / sum(eigs)*100,2), '% variance explained)'),
       y=paste0('PC2 (', round(eigs[2] / sum(eigs)*100,2), '% variance explained)'))
# ggsave("../results/exploratory/PCA_with_links_onlypaired.png", height = 5)

## it looks like there is a somewhat coordinated shift in PC2
## for which patients the PC2 of arx higher than the PC2 of rlps?
## there is the problem with weird samples so 
df_only_paired$col_coordinated_pc2 = sapply((df_only_paired$patient), function(i){
  any(df[(df$patient == i) & (df$group == 'arx'),'PC2'] < df[(df$patient == i) & (df$group == 'rlps'),'PC2'])
})
df_only_paired$col_coordinated_pc1 = sapply((df_only_paired$patient), function(i){
  any(df[(df$patient == i) & (df$group == 'arx'),'PC1'] < df[(df$patient == i) & (df$group == 'rlps'),'PC1'])
})
df_only_paired$col_coordinated_pc3 = sapply((df_only_paired$patient), function(i){
  any(df[(df$patient == i) & (df$group == 'arx'),'PC3'] < df[(df$patient == i) & (df$group == 'rlps'),'PC3'])
})
df_only_paired$col_coordinated_pc4 = sapply((df_only_paired$patient), function(i){
  any(df[(df$patient == i) & (df$group == 'arx'),'PC4'] < df[(df$patient == i) & (df$group == 'rlps'),'PC4'])
})

ggplot(df_only_paired,
       aes(x=PC1, y=PC2, alpha=col_coordinated_pc2))+geom_point(aes(col=group))+geom_line(aes(group=patient))+
  labs(x=paste0('PC1 (', round(eigs[1] / sum(eigs)*100,2), '% variance explained)'),
       y=paste0('PC2 (', round(eigs[2] / sum(eigs)*100,2), '% variance explained)'))

ggplot()+
  geom_point(data=df_only_paired[,1L:(ncol(df_only_paired)-2)], aes(x=PC1, y=PC2), alpha=0.1)+
  geom_point(data=df_only_paired,
             aes(x=PC1, y=PC2, col=group, group=))+
  geom_line(data=df_only_paired, aes(x=PC1, y=PC2, group=patient))+
  facet_wrap(.~interaction(col_coordinated_pc2, col_coordinated_pc1))+
  labs(x=paste0('PC1 (', round(eigs[1] / sum(eigs)*100,2), '% variance explained)'),
       y=paste0('PC2 (', round(eigs[2] / sum(eigs)*100,2), '% variance explained)'))+
  theme_bw()

ggplot()+
  geom_point(data=df_only_paired[,1L:(ncol(df_only_paired)-2)], aes(x=PC1, y=PC3), alpha=0.1)+
  geom_point(data=df_only_paired,
             aes(x=PC1, y=PC3, col=group, group=))+
  geom_line(data=df_only_paired, aes(x=PC1, y=PC3, group=patient))+
  facet_wrap(.~interaction(col_coordinated_pc2, col_coordinated_pc1))+
  labs(x=paste0('PC1 (', round(eigs[1] / sum(eigs)*100,2), '% variance explained)'),
       y=paste0('PC3 (', round(eigs[3] / sum(eigs)*100,2), '% variance explained)'))+
  theme_bw()
# ggsave("../results/exploratory/PCA_with_links_onlypaired_facets_shifts.png", height = 5, width = 5)

ggplot()+
  geom_point(data=df_only_paired[,1L:(ncol(df_only_paired)-2)], aes(x=PC3, y=PC4), alpha=0.1)+
  geom_point(data=df_only_paired,
             aes(x=PC3, y=PC4, col=group, group=))+
  geom_line(data=df_only_paired, aes(x=PC3, y=PC4, group=patient))+
  facet_wrap(.~interaction(col_coordinated_pc3, col_coordinated_pc4))+
  labs(x=paste0('PC3 (', round(eigs[3] / sum(eigs)*100,2), '% variance explained)'),
       y=paste0('PC4 (', round(eigs[4] / sum(eigs)*100,2), '% variance explained)'))+
  theme_bw()
# ggsave("../results/exploratory/PCA_with_links_onlypaired_facets_shifts_pc3_pc4.png", height = 5, width = 5)

## PC4 seems extremely indicative!

## rotation

Compositional::bivt.contour(pca)
pca$rotation

biplot(pca, choices = c(3,4), scale = FALSE, pc.biplot = TRUE)

ggplot()+
  geom_point(data=cbind.data.frame((pca$x), group=df$group), aes(x=PC3, y=PC4, col=group), alpha=.5)+
  geom_segment(data=cbind.data.frame(sig=rownames(pca$rotation), pca$rotation),
               aes(x=0, y=0, xend=PC3*4, yend=PC4*4), arrow = arrow(length = unit(0.03, "npc")), col='red', size=1)+
  geom_label(data=cbind.data.frame(sig=rownames(pca$rotation), pca$rotation),
             aes(x=PC3*3, y=PC4*3, label=sig)) +theme_bw()+ theme(legend.position = "bottom")
ggsave("../results/exploratory/PCA_loadings_pc3_pc4.png", height = 5, width = 5)

## the same but only keeping the paired datapoints
ggplot()+
  geom_point(data=cbind.data.frame((pca_only_paired$x), group=df_only_paired$group), aes(x=PC3, y=PC4, col=group), alpha=.5)+
  geom_segment(data=cbind.data.frame(sig=rownames(pca_only_paired$rotation), pca_only_paired$rotation),
               aes(x=0, y=0, xend=PC3*4, yend=PC4*4), arrow = arrow(length = unit(0.03, "npc")), col='red', size=1)+
  geom_label(data=cbind.data.frame(sig=rownames(pca_only_paired$rotation), pca_only_paired$rotation),
             aes(x=PC3*3, y=PC4*3, label=sig)) +theme_bw()+ theme(legend.position = "bottom")
ggsave("../results/exploratory/PCA_loadings_pc3_pc4_onlypaired.png", height = 5, width = 5)

## add directionality to arrows, from arx to rlps
df_only_paired_segments = do.call('rbind.data.frame', lapply(df_only_paired$patient, function(i){
  do.call('rbind.data.frame', apply(df_only_paired[ (df_only_paired$patient == i) & (df_only_paired$group == 'arx'),], 1, function(rw_arx){
    do.call('rbind.data.frame', apply(df_only_paired[ (df_only_paired$patient == i) & (df_only_paired$group == 'rlps'),], 1, function(rw_rlps){
      cbind.data.frame(x=rw_arx['PC3'], y=rw_arx['PC4'], xend=rw_rlps['PC3'], yend=rw_rlps['PC4'])
    }))
  }))
  ## now do the same from rlps to arx
  })
)
df_only_paired_segments = df_only_paired_segments[!duplicated(df_only_paired_segments),]
df_only_paired_segments_names_samples = rownames(df_only_paired_segments)
df_only_paired_segments = data.frame(apply(df_only_paired_segments, 2, as.numeric))
df_only_paired_segments$sample = df_only_paired_segments_names_samples
ggplot()+
  geom_point(data=cbind.data.frame((pca_only_paired$x), group=df_only_paired$group), aes(x=PC3, y=PC4, col=group), alpha=.5)+
  geom_segment(data=cbind.data.frame(sig=rownames(pca_only_paired$rotation), pca_only_paired$rotation),
               aes(x=0, y=0, xend=PC3*4, yend=PC4*4), arrow = arrow(length = unit(0.03, "npc")), col='red', size=1)+
  geom_segment(data=data.frame(df_only_paired_segments), aes(x=x, y=y, xend=xend, yend=yend),
               arrow = arrow(length = unit(0.02, "npc")), alpha=.2)+
  geom_label(data=cbind.data.frame(sig=rownames(pca_only_paired$rotation), pca_only_paired$rotation),
             aes(x=PC3*3, y=PC4*3, label=sig)) +theme_bw()+ theme(legend.position = "bottom")
ggsave("../results/exploratory/PCA_loadings_pc3_pc4_loadings.png", height = 5, width = 5)

ggplot()+
  geom_point(data=cbind.data.frame((pca_only_paired$x), group=df_only_paired$group), aes(x=PC3, y=PC4, col=group), alpha=.5)

df_only_paired_segments$col_coordinated_pc3 = (df_only_paired_segments$x-df_only_paired_segments$xend)<0
df_only_paired_segments$col_coordinated_pc4 = (df_only_paired_segments$y-df_only_paired_segments$yend)<0

## This plot is not perfect, probably because in the col_coordinated_pc3 I have used the "any" strategy (grouping them by patient)
ggplot()+
  geom_point(data=df_only_paired[,1L:(ncol(df_only_paired)-2)], aes(x=PC3, y=PC4), alpha=0.1)+
  geom_point(data=df_only_paired,
             aes(x=PC3, y=PC4, col=group, group=group))+
  geom_segment(data=data.frame(df_only_paired_segments), aes(x=x, y=y, xend=xend, yend=yend),
               arrow = arrow(length = unit(0.01, "npc")))+
  facet_wrap(.~interaction(col_coordinated_pc3, col_coordinated_pc4))+
  labs(x=paste0('PC3 (', round(eigs[3] / sum(eigs)*100,2), '% variance explained)'),
       y=paste0('PC4 (', round(eigs[4] / sum(eigs)*100,2), '% variance explained)'))+
  theme_bw()

ggplot()+
  geom_bar(data=(df_only_paired_segments),
           aes(x=factor(sample, levels=df_only_paired_segments$sample[order(df_only_paired_segments$y - df_only_paired_segments$yend)]),
               y=y-yend), stat = 'identity')


do.call(grid.arrange, c(ncol=7, lapply(colnames(exposures),
       function(sig_it) ggplot(melt(split(df[,sig_it], df$group)), aes(x=L1, y=value))+geom_boxplot()+ggtitle(sig_it))))

remove_organoids = function(df){
  df[!grepl('Sample ', rownames(df)),]
}

give_zero = function(a) a[a == 0]
give_infty = function(a) a[is.infinite(a)]
give_nonzero = function(a) a[(a > 0) & !is.infinite(a)]

exposures_archival = exposures[rownames(exposures) %in% as.character(patient.meta$SAMPLE_ID[patient.meta$group == "arx"]),]
exposures_relapse = exposures[rownames(exposures) %in% as.character(patient.meta$SAMPLE_ID[patient.meta$group == "rlps"]),]

pdf("../results/exploratory/logR_plots.pdf", height = 2)
for(j in 1:7){
  .archival = exposures_archival[,j] / rowSums(exposures_archival[,(1:7)[!c(1:7 %in% j)]])
  .archival_zero = give_zero(.archival)
  .archival_infty = give_infty(.archival)
  .archival_nonzero = give_nonzero(.archival)
  .relapse = exposures_relapse[,j] / rowSums(exposures_relapse[,(1:7)[!c(1:7 %in% j)]])
  .relapse_zero = give_zero(.relapse)
  .relapse_infty = give_infty(.relapse)
  .relapse_nonzero = give_nonzero(.relapse)
  
  grid.arrange(
    ggplot(melt(cbind.data.frame(zero=c(length(.archival_zero), length(.relapse_zero)),
                                 total=c(length(c(.archival_infty,.archival_nonzero)), length(c(.relapse_infty,.relapse_nonzero))),
                                 cohort=c('Archival', ' Relapse')), id.vars = "cohort"), aes(x=cohort, fill=variable, y=value))+geom_bar(stat = "identity")+
      facet_wrap(.~cohort, scales = "free")+theme(legend.position = "bottom", legend.title = element_blank(), legend.key = element_blank(), legend.key.width = unit(.01, "cm"), legend.key.height = unit(.1, "cm"), strip.background = element_blank(), strip.text.x = element_blank())+labs(x='', y=''),
    ggplot(cbind.data.frame(x=c(.relapse_nonzero, .archival_nonzero), y=c(rep("Relapse", length(.relapse_nonzero)), rep("Archival", length(.archival_nonzero)))))+
      geom_density(aes(x=log(x), col=y))+theme(legend.position = "bottom", legend.title = element_blank(), legend.key.width = unit(.3, "cm"), legend.key.height = unit(.1, "cm"))+scale_color_manual(values=c("#56B4E9", "#E69F00")),
    ggplot(melt(cbind.data.frame(zero=c(length(.archival_infty), length(.relapse_infty)),
                                 total=c(length(c(.archival_zero,.archival_nonzero)), length(c(.relapse_zero,.relapse_nonzero))),
                                 cohort=c('Archival', ' Relapse')), id.vars = "cohort"), aes(x=cohort, fill=variable, y=value))+geom_bar(stat = "identity")+
      facet_wrap(.~cohort, scales = "free")+theme(legend.position = "bottom", legend.title = element_blank(), legend.key = element_blank(), legend.key.width = unit(.01, "cm"), legend.key.height = unit(.1, "cm"), strip.background = element_blank(), strip.text.x = element_blank())+labs(x='', y=''),
    ncol=3, widths=c(0.2, 0.6, 0.2), top=paste0('Signature ', j))
}

dev.off()

remove_inf = function(i){
  i = i[!is.infinite(i)]
  i[!is.nan(i)]
}

source("header.R")
obj = list(Y = exposures,
           num_individuals = length(unique(patient.meta$PATIENT_ID)),
           d = 7,
           n = nrow(exposures),
           x = cbind(1, as.numeric(as.factor(patient.meta$group))-1),
           z = give_z_matrix_from_labels(lbls = patient.meta$PATIENT_ID))



give_LRchanges_barplot(obj, nrow_facets = 1)
ggsave("../results/exploratory/paired_LR_changes.pdf", width = 12, height = 2)

good_exposures = readRDS("~/Desktop/good_exposures.RDS")
all(exposures == good_exposures)

## hclust
source("../../../Vias_Brenton/copy_number_analysis_organoids/helper_functions.R")
require(dendextend)
require(ggdendro)
require(ggrepel)
extra_expand = 0.02
extra_expand_v2 = 0.02
plt_hierach <- give_dendrogram_from_imputation(impute_VALUE = 1e-2, exposures = exposures, plot = T, return_grob = T)
plot(plt_hierach)

plt_hierach <- give_dendrogram_from_imputation(impute_VALUE = 4e-2, exposures = exposures, plot = T, return_grob = T)
plot(plt_hierach)

# exposures
# 
# archival <- read.csv("../../data/restricted/Archival_samples.csv")
# relapse <- read.csv("../../data/restricted/Biopsies_relapse.csv")
# brca <- read.csv("../../data/restricted/brca.csv")
# 
# cbind(archival$IM_ID[match(brca$TRIALNO, archival$BritRoc.No)],
#       relapse$JBLAB_ID[match(brca$TRIALNO, relapse$BritRoc.No)])
# brca
# 
# patient.meta


