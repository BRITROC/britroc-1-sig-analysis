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
component_by_signature<-readRDS(paste(this_path,"data/feat_sig_mat.rds",sep="/"))


sigs = quantifySignatures_alt(SxCmatrix)

## remove samples with zero s5, as nothing will change in those
sigs = sigs[, !(sigs['s5',] == 0)]


## simulate SxC
## simulate exposures with zero s5
exposures_true_sim = rbind(runif(100),runif(100),runif(100),runif(100),0,runif(100),runif(100))
SxC_sim = component_by_signature %*% exposures_true_sim
exposures_sim = quantifySignatures_alt(t(SxC_sim))

ggplot(cbind.data.frame(inferred_sig=as.vector(exposures_sim),
                        true_sig=as.vector(t(normalise_cl(exposures_true_sim))),
                        sig=rep(paste0('s', 1:7), ncol(exposures_sim)),
                        sample=rep(1:ncol(exposures_sim), each=nrow(exposures_sim))),
       aes(x=inferred_sig, y=true_sig))+geom_point()+facet_wrap(.~sig, nrow=1)+
  geom_abline(slope = 1, intercept = 0, lty='dashed')
ggsave("../../results/exploratory/bleeding_signature_s5/extractionwithouts5_recovery_simulation.pdf",
       width = 10, height = 2)


plts = lapply(1:7, function(sig_it){
  exposures_true_sim = rbind(runif(100),runif(100),runif(100),runif(100),runif(100), runif(100), runif(100))
  exposures_true_sim[sig_it,] = 0
  SxC_sim = component_by_signature %*% exposures_true_sim
  exposures_sim = quantifySignatures_alt(t(SxC_sim))
  
  ggplot(cbind.data.frame(inferred_sig=as.vector(exposures_sim),
                          true_sig=as.vector(t(normalise_cl(exposures_true_sim))),
                          sig=rep(paste0('s', 1:7), ncol(exposures_sim)),
                          sample=rep(1:ncol(exposures_sim), each=nrow(exposures_sim))),
         aes(x=inferred_sig, y=true_sig))+geom_point()+facet_wrap(.~sig, nrow=1)+
    geom_abline(slope = 1, intercept = 0, lty='dashed')
})

pdf("../../results/exploratory/bleeding_signature_s5/extractionwithoutsigs_recovery_simulation.pdf", height=1.5)
for(i in 1:7) print(plts[i])
dev.off()



### With matrices
exposures_true_sim[-5,]
exposures_sim[-5,]

solve(exposures_true_sim[-5,])
