rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

segs = read.table("../../data/britroc_30kb_ds_absCopyNumber_segmentTable.tsv", header = T)
library(reshape2)
library(ggplot2)
library(dplyr)
library(gridExtra)
source("functions.R")
source("../../../britroc-cnsignatures-bfb69cd72c50/main_functions.R")
source("../../../britroc-cnsignatures-bfb69cd72c50/helper_functions.R")
source("../../../../GlobalDA/code/2_inference_TMB/helper_TMB.R")
components = readRDS("../../../britroc-cnsignatures-bfb69cd72c50/data/component_parameters.rds")

SxCmatrix = readRDS("../../out/signature_extraction/SxCmatrix.RDS")

sigs = quantifySignatures_alt(SxCmatrix)

## remove samples with zero s5, as nothing will change in those
sigs = sigs[, !(sigs['s5',] == 0)]

sigs_nos5 = quantifySignatures_alt(SxCmatrix, subsample_sigs = c(1:4, 6:7))
sigs_nos5 = sigs_nos5[,colnames(sigs_nos5) %in% colnames(sigs)]
# sigs_nos5 = rbind(sigs_nos5[1:4,], NA, sigs_nos5[5:6,])
# sigs_nos5

sigs_subcomp = t(normalise_cl(sigs[-5,]))

dim(sigs)
dim(sigs_nos5)
dim(sigs_subcomp)

ggplot(melt(list(sigs, sigs_nos5)), aes(x=interaction(Var1,Var2)))

# ## in some samples they are exactly the same
# dim(alrs7_subcomp[rowSums(alrs7_subcomp == alrs7_nos5),])
# dim(alrs7_subcomp)
# 
# subset_samples_all_same = rownames(alrs7_subcomp[rowSums(alrs7_subcomp == alrs7_nos5) == ncol(alrs7_nos5),])
# subset_samples_all_same = subset_samples_all_same[!is.na(subset_samples_all_same)]
# 
# ## go back to the original signatures. they are the ones with zero s5 exposure
# sigs['s5',subset_samples_all_same]
# 
# ## remove those samples as they are not informative
# 
# sigs = sigs[,!(colnames(sigs) %in% subset_samples_all_same)]
# alrs7_subcomp = alrs7_subcomp[!(rownames(alrs7_subcomp) %in% subset_samples_all_same),]


ggplot(cbind.data.frame(sigs=melt(sigs_subcomp),
                 sigs_nos5=melt(sigs_nos5)), aes(x=sigs.value, y=sigs_nos5.value))+geom_point()+facet_wrap(.~sigs.Var1, nrow=1)+
labs(x='Subcompositional of exposures with all signatures', y='Exposures with extraction without s5')+
  geom_abline(slope = 1, intercept = 0, lty='dashed')
ggsave("../../results/exploratory/bleeding_signature_s5/scatterplot_extractionwithouts5.pdf", width = 12, height = 3)


## Compare all log-ratios in the two conditions
sixchoose2 = combinat::combn(x = 1:6, 2)
all_logR_sigs_subcomp = apply(sigs_subcomp, 2, function(i) apply(sixchoose2, 2, function(j) log(i[j[1]]/i[j[2]])  ))
all_logR_sigs_nos5 = apply(sigs_nos5, 2, function(i) apply(sixchoose2, 2, function(j) log(i[j[1]]/i[j[2]])  ))

ggplot(cbind.data.frame(logR_subcomp=as.vector(all_logR_sigs_subcomp), logR_nosig5=as.vector(all_logR_sigs_nos5),
                        col=rep(apply(sixchoose2, 2, paste0, collapse='_'), ncol(all_logR_sigs_nos5))),
        aes(x=logR_subcomp, y=logR_nosig5, col=col))+geom_point()+geom_abline(intercept = 0, slope = 1, lty='dashed')+
  facet_wrap(.~col)
ggsave("../../results/exploratory/bleeding_signature_s5/all_logR_change.pdf", width = 12, height = 12)

## compute the perturbations
as(compositions::alr(t(sigs_subcomp[6:1,])), 'matrix')
as(compositions::alr(t(sigs_nos5[6:1,])), 'matrix')

perturbations = t(sigs_nos5)/t(sigs_subcomp)
# perturbations_norm = normalise_rw(perturbations)
normalise_rw(t(sigs_subcomp) * perturbations)
t(sigs_nos5)

apply(perturbations, 2, boxplot)
boxplot(log(perturbations))

ggplot(melt(log(perturbations)), aes(x=Var2, y=value))+geom_density()+geom_jitter()+ggtitle('Normalised perturbation from\nincluding to excluding s5 ')+
  labs(x='Signature', y='Normalised perturbation')+theme_bw()
ggplot(melt(log(perturbations)), aes(x=Var2, y=value))+geom_boxplot()+ggtitle('Normalised perturbation from\nincluding to excluding s5 ')+
  labs(x='Signature', y='Normalised perturbation')+theme_bw()
ggsave("../../results/exploratory/bleeding_signature_s5/boxplot_extractionwithouts5.pdf", width = 3, height = 3)
ggplot(melt(log(perturbations)), aes(x=Var2, y=value))+geom_boxplot()+
  labs(x='Signature', y='Normalised perturbation')+theme_bw()
ggsave("../../results/exploratory/bleeding_signature_s5/boxplot_extractionwithouts5_notitle.pdf", width = 3, height = 3)

pairs(log(perturbations))

alrs7_subcomp = as(compositions::alr(t(sigs_subcomp)[,6:1]), 'matrix')
alrs7_nos5 = as(compositions::alr(t(sigs_nos5)[,6:1]), 'matrix')
pairs(alrs7_nos5 - alrs7_subcomp )

ggplot(cbind.data.frame(alrs7_original_subcomposition=alrs7_subcomp[,1], alrs7_without_s5=alrs7_nos5[,1],
                        perturbation_s1=perturbations[,1]),
       aes(x=alrs7_original_subcomposition, y=alrs7_without_s5, col=1-perturbation_s1))+geom_point()+
  labs(x="log(s7/s1) including s5", y="log(s7/s1) excluding s5")+
  geom_abline(slope = 1, intercept = 0, lty='dashed')+theme_bw()+guides(col=FALSE)
ggsave("../../results/exploratory/bleeding_signature_s5/extractionwithouts5_alr_s7s1.pdf", width = 3, height = 3)

##' The alr of s7 wrt s1 is constant.
##' (a) neither s1 nor s7 are affected by s5 (most likely)
##'     perform analyses wothout s1, s7
##' (b) s1 and s7 get a proportional bleed from s5

sigs_nos5[c(2:5),]
sigs_subcomp[c(2:5),]
pairs(t(sigs_nos5[c(2:5),]))
pairs(t(sigs_subcomp[c(2:5),]))
pairs(t(sigs_nos5[c(2:5),]-sigs_subcomp[c(2:5),]))

boxplot(alrs7_nos5[,5:2] - alrs7_subcomp[,5:2])
pairs(alrs7_nos5[,5:2] - alrs7_subcomp[,5:2])


ggplot(cbind.data.frame(alrs_original_subcomposition=log(sigs_subcomp['s6',]/sigs_subcomp['s4',]),
                        alrs_without_s5=log(sigs_nos5['s6',]/sigs_nos5['s4',])),
       aes(x=alrs_original_subcomposition, y=alrs_without_s5))+geom_point()+
  labs(x="log(s6/s4) including s5", y="log(s6/s4) excluding s5")+theme_bw()+
  geom_abline(slope = 1, intercept = 0, lty='dashed')
ggsave("../../results/exploratory/bleeding_signature_s5/extractionwithouts5_alr_s6s4.pdf", width = 3, height = 3)

ggplot(cbind.data.frame(alrs_original_subcomposition=log(sigs_subcomp['s4',]/sigs_subcomp['s6',]),
                        alrs_without_s5=log(sigs_nos5['s4',]/sigs_nos5['s6',])),
       aes(x=alrs_original_subcomposition, y=alrs_without_s5))+geom_point()+
  labs(x="log(s4/s6) including s5", y="log(s4/s6) excluding s5")+theme_bw()+
  geom_abline(slope = 1, intercept = 0, lty='dashed')
ggsave("../../results/exploratory/bleeding_signature_s5/extractionwithouts5_alr_s4s6.pdf", width = 3, height = 3)

ggplot(cbind.data.frame(alrs_original_subcomposition=log(sigs_subcomp['s6',]/(sigs_subcomp['s1',]+sigs_subcomp['s7',])),
                        alrs_without_s5=log(sigs_nos5['s6',]/(sigs_nos5['s1',]+sigs_nos5['s7',]))),
       aes(x=alrs_original_subcomposition, y=alrs_without_s5))+geom_point()+
  labs(x="log(s6/(s1+s7)) including s5", y="log(s6/(s1+s7)) excluding s5")+theme_bw()+
  geom_abline(slope = 1, intercept = 0, lty='dashed')


ggplot(cbind.data.frame(alrs_original_subcomposition=log(sigs_subcomp['s3',]/(sigs_subcomp['s1',]+sigs_subcomp['s7',])),
                        alrs_without_s5=log(sigs_nos5['s3',]/(sigs_nos5['s1',]+sigs_nos5['s7',]))),
       aes(x=alrs_original_subcomposition, y=alrs_without_s5))+geom_point()+
  labs(x="log(s3/(s1+s7)) including s5", y="log(s3/(s1+s7)) excluding s5")+theme_bw()+
  geom_abline(slope = 1, intercept = 0, lty='dashed')

ggplot(cbind.data.frame(alrs_original_subcomposition=log(sigs_subcomp['s3',]/(sigs_subcomp['s1',])),
                        alrs_without_s5=log(sigs_nos5['s3',]/(sigs_nos5['s1',]))),
       aes(x=alrs_original_subcomposition, y=alrs_without_s5))+geom_point()+
  labs(x="log(s3/(s1)) including s5", y="log(s3/(s1)) excluding s5")+theme_bw()+
  geom_abline(slope = 1, intercept = 0, lty='dashed')

ggplot(cbind.data.frame(alrs_original_subcomposition=log(sigs_subcomp['s3',]/colSums(sigs_subcomp[rownames(sigs_subcomp)[!(rownames(sigs_subcomp) %in% 's3')],])),
                        alrs_without_s5=log(sigs_nos5['s3',]/colSums(sigs_nos5[rownames(sigs_nos5)[!(rownames(sigs_nos5) %in% 's3')],]))),
       aes(x=alrs_original_subcomposition, y=alrs_without_s5))+geom_point()+
  labs(x="log(s3/(s1+s7)) including s5", y="log(s3/(s1+s7)) excluding s5")+theme_bw()+
  geom_abline(slope = 1, intercept = 0, lty='dashed')

ggplot(cbind.data.frame(alrs_original_subcomposition=log(sigs_subcomp['s2',]/colSums(sigs_subcomp[rownames(sigs_subcomp)[!(rownames(sigs_subcomp) %in% 's3')],])),
                        alrs_without_s5=log(sigs_nos5['s3',]/colSums(sigs_nos5[rownames(sigs_nos5)[!(rownames(sigs_nos5) %in% 's3')],]))),
       aes(x=alrs_original_subcomposition, y=alrs_without_s5))+geom_point()+
  labs(x="log(s3/(s1+s7)) including s5", y="log(s3/(s1+s7)) excluding s5")+theme_bw()+
  geom_abline(slope = 1, intercept = 0, lty='dashed')

do.call('grid.arrange', lapply(rownames(sigs_subcomp), function(i){
  ggplot(cbind.data.frame(alrs_original_subcomposition=log(sigs_subcomp[i,]/colSums(sigs_subcomp[rownames(sigs_subcomp)[!(rownames(sigs_subcomp) %in% i)],])),
                          alrs_without_s5=log(sigs_nos5[i,]/colSums(sigs_nos5[rownames(sigs_nos5)[!(rownames(sigs_nos5) %in% i)],]))),
         aes(x=alrs_original_subcomposition, y=alrs_without_s5))+geom_point()+
    labs(x="log(i/(s1+s7)) including s5", y="log(s3/(s1+s7)) excluding s5")+theme_bw()+
    geom_abline(slope = 1, intercept = 0, lty='dashed')
}))

df_coordinated_change = droplevels(cbind.data.frame(difference_alr_s4s6 = log(abs(log(sigs_subcomp['s6',]/sigs_subcomp['s4',])-log(sigs_nos5['s6',]/sigs_nos5['s4',]))),
                            difference_alr_s7s1= log(abs(alrs7_subcomp[,1]- alrs7_nos5[,1]))))

## check if the samples that line up perfectly for (s1, s7) are the same than for (s4,s6)
## they are not. The plot is all over the place
ggplot(df_coordinated_change,
                        aes(x=difference_alr_s4s6, y=difference_alr_s7s1))+geom_point()+
  theme_bw()+
  geom_abline(slope = 1, intercept = 0, lty='dashed')


ggplot(melt(log(perturbations[,c('s2','s3')])), aes(x=Var2, y=value))+geom_violin()+ggtitle('Normalised perturbation from\nincluding to excluding s5 ')+
  labs(x='Signature', y='Normalised perturbation')+theme_bw()

ggplot(cbind.data.frame((perturbations[,c('s2','s3')])), aes(x=s2, y=s3))+geom_point()+ggtitle('Normalised perturbation from\nincluding to excluding s5 ')+
  labs(x='Signature', y='Normalised perturbation')+theme_bw()

## see if the alr of s4 or s6 in the nonsig5 group depends on the quantity of s5
plot(alrs7_nos5[,c('s4')]-log(sigs['s5',]/sigs['s1',]))
plot(alrs7_nos5[,c('s6')]-log(sigs['s5',]/sigs['s1',]))
plot(alrs7_nos5[,c('s2')]-log(sigs['s5',]/sigs['s1',]))
plot(alrs7_nos5[,c('s3')]-log(sigs['s5',]/sigs['s1',])-alrs7_subcomp[,'s3'])

## If s1 doesn't change, then the alr wrt s1 tells you the change betweeen no including and excluding s5
boxplot(alrs7_subcomp)
ggplot(melt(alrs7_subcomp), aes(x=Var2, y=value))+geom_boxplot()+geom_jitter(alpha=0.3)+geom_hline(yintercept = 0, col='red')

perturbations_renormalised = sweep(perturbations, 1, perturbations[,1], '/')
perturbations_renormalised_no_outliers = perturbations_renormalised
perturbations_renormalised_no_outliers[perturbations_renormalised_no_outliers > 4] = NA
perturbations_renormalised_no_outliers[perturbations_renormalised_no_outliers < -4] = NA
pairs(perturbations_renormalised_no_outliers)

## see to which the amount increased in s4+s6 is due
plot(log(perturbations[,c('s4', 's6')]))

pairs(cbind.data.frame(perturbations_renormalised[,'s1'],
perturbations_renormalised[,'s3'],
perturbations_renormalised[,'s7']))

ggplot(cbind.data.frame(pert_s3=perturbations_renormalised[,'s3'],
                 pert_s7=perturbations_renormalised[,'s7']), aes(x=pert_s3, y=pert_s7))+geom_point()

boxplot(log(perturbations_renormalised))

ggplot(melt((alrs7_nos5-alrs7_subcomp)[,c('s2', 's3', 's4', 's6')]), aes(x=Var1, y=value, fill=Var2))+
  geom_bar(stat = 'identity')+facet_wrap(.~Var2)

normalised_alrchange_subset = sweep((alrs7_nos5-alrs7_subcomp)[,c('s2', 's3', 's4', 's6')], 1,
                                    rowSums((alrs7_nos5-alrs7_subcomp)[,c('s2', 's3', 's4', 's6')]),
                                    '/')
normalised_alrchange_subset = normalised_alrchange_subset[!apply(normalised_alrchange_subset, 1, function(i) any(is.na(i))),]
ggplot(melt(normalised_alrchange_subset), aes(x=Var1, y=value, fill=Var2))+
  geom_bar(stat = 'identity')+theme_minimal()
ggsave("../../results/exploratory/bleeding_signature_s5/extractionwithouts5_relative_allocation.pdf", width = 14, height = 2)

# plot(normalised_alrchange_subset[,'s4'],
#      normalised_alrchange_subset[,'s6'])
# plot(normalised_alrchange_subset[,'s3'],
#      normalised_alrchange_subset[,'s4'])
# 
# 
# plot(normalised_alrchange_subset[,'s2'],
#      normalised_alrchange_subset[,'s4'])
# 

## looking only at: sum of s4+s6, s2, and s3
tripartite_analysis_df_includings5 = cbind.data.frame(s2=sigs['s2',], s3=sigs['s3',], s4s6 = (sigs['s4',]+sigs['s6',]))
tripartite_analysis_df_excluding5 = cbind.data.frame(s2=sigs_nos5['s2',], s3=sigs_nos5['s3',], s4s6 = (sigs_nos5['s4',]+sigs_nos5['s6',]))

pairs(tripartite_analysis_df_includings5)
pairs(normalise_rw(tripartite_analysis_df_includings5))
pairs(tripartite_analysis_df_excluding5)
pairs(normalise_rw(tripartite_analysis_df_excluding5))


plot(normalise_rw(tripartite_analysis_df_excluding5)[,'s3'],
normalise_rw(tripartite_analysis_df_excluding5)[,'s4s6'])


