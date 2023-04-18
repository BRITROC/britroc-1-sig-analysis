rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(transport)
library(dplyr)
library(ggplot2)
library(reshape2)
library(qgraph)
source("../../../britroc-cnsignatures-bfb69cd72c50/main_functions.R")
source("../../../britroc-cnsignatures-bfb69cd72c50/helper_functions.R")

SxCmatrix = readRDS("../../out/signature_extraction/SxCmatrix.RDS")

quantifySignatures_alt = function(sample_by_component,component_by_signature=NULL)
{
  if(is.null(component_by_signature))
  {
    component_by_signature<-readRDS(paste(this_path,"data/feat_sig_mat.rds",sep="/"))
  }
  signature_by_sample<-YAPSA::LCD(t(sample_by_component),
                                  YAPSA:::normalize_df_per_dim(component_by_signature,2))
  signature_by_sample<-normaliseMatrix(signature_by_sample, sig_thresh = 0)
  signature_by_sample
}

sigs = quantifySignatures_alt(SxCmatrix)
sigs


## super naive first approach
.SxCmatrix_mod = SxCmatrix+1
.sigs_mod = quantifySignatures_alt(.SxCmatrix_mod)

plot(as.vector(.sigs_mod), as.vector(sigs))
ggplot(melt(.sigs_mod/ sigs), aes(x=log(value)))+geom_density()+facet_wrap(.~Var1, scales = "free")
logR_change = log(.sigs_mod/ sigs)
## Positive value: this signature gains exposure. Negative value: this signature loses exposure
apply(logR_change, 1, function(i) log(sum(i>0, na.rm = T)/sum(i<0, na.rm = T)))

sigs = sigs[,1:3]
.sigs_mod = .sigs_mod[,1:3]
cost_matrix = cost_matrix_numbers = matrix(1e8, length(as.vector(sigs)), length(as.vector(sigs)))
for(i in 1:ncol(sigs)){
  first_entry = 1+((i-1)*nrow(sigs))
  cost_matrix[ first_entry:(first_entry+6), first_entry:(first_entry+6)] = 1
  cost_matrix_numbers[ first_entry:(first_entry+6), first_entry:(first_entry+6)] = i
}
# cost_matrix = t(cost_matrix)
# cost_matrix_numbers = t(cost_matrix_numbers)

vec1 = as.vector(sigs); names(vec1) = as.vector(outer(rownames(sigs), colnames(sigs), Vectorize(function(i,j) paste0(i, '_', j))))
vec2 = as.vector(.sigs_mod); names(vec2) = names(vec1)
colnames(cost_matrix) = rownames(cost_matrix) = names(vec1)
opttrans = transport(a = vec1, b = vec2,
                     costm = cost_matrix)
# opttrans$from_sig = ( (opttrans$from) %% (nrow(sigs))); opttrans$from_sig[opttrans$from_sig == 0] = 7
# opttrans$to_sig = ( (opttrans$to) %% (nrow(sigs))); opttrans$to_sig[opttrans$to_sig == 0] = 7
opttrans$from_sig = rep(rownames(sigs), ncol(sigs))[opttrans$from]
opttrans$to_sig = rep(rownames(sigs), ncol(sigs))[opttrans$to]
opttrans$sample = rep(colnames(sigs), nrow(sigs))[opttrans$from]
opttrans$sample2 = rep(colnames(sigs), nrow(sigs))[opttrans$to]
# opttrans$sample = 1+(opttrans$from %/% nrow(sigs))
# opttrans$sample2 = 1+(opttrans$from %/% nrow(sigs))
head(opttrans)
View(opttrans)

## doesn't make sense. it should split in components
qgraph(data.frame(opttrans$from, opttrans$to))
qgraph(data.frame(paste0(gsub('IM_', '', opttrans$sample), '_', opttrans$from_sig),
                  paste0(gsub('IM_', '', opttrans$sample), '_', opttrans$to_sig)), label.cex=2)

## is this a coincidence?
max(c(opttrans$from, opttrans$to))
sum(sigs > 0)

fakecost_matrix = matrix(0, length(as.vector(sigs)), length(as.vector(sigs)))
fakecost_matrix[c(opttrans$to, opttrans$from)] = 1
image(fakecost_matrix)

sum(opttrans$mass)
##' for each signature, see where the signal goes. we have to do the average to each of the remaining signatures,
##' once we have normalised them by the total mass of the signature
##' Within each patient, the probabilities add up to one
optttrans_per_sig_mass = opttrans %>% group_by(from_sig) %>% mutate(relative_mass=mass/sum(mass)) %>% group_by(from_sig,to_sig) %>% summarise(average_mass=mean(relative_mass))
optttrans_per_sig_mass2 = opttrans %>% group_by(from_sig) %>% group_by(from_sig,to_sig) %>% summarise(average_mass=mean(mass))
opttrans
pheatmap(opttrans)

give_names= function(m){
  rownames(m) = colnames(m) = paste0('s', 1:7); return(m)
}

pheatmap::pheatmap(give_names(matrix(unlist(optttrans_per_sig_mass[,3]), nrow=7)), cluster_rows = F, cluster_cols = F, scale = "row")
pheatmap::pheatmap(give_names(matrix(unlist(optttrans_per_sig_mass2[,3]), nrow=7)), cluster_rows = F, cluster_cols = F, scale = "row")

transition_matrix = give_names(matrix(unlist(optttrans_per_sig_mass2[,3]), nrow=7))
# transition_matrix = sweep(transition_matrix, 1, rowSums(transition_matrix), '/') ## normalise by from node
transition_matrix = sweep(transition_matrix, 2, rowSums(transition_matrix), '/') ## normalise by to node

qgraph(data.frame(to=optttrans_per_sig_mass2$to_sig, from=optttrans_per_sig_mass2$from_sig, thickness=(optttrans_per_sig_mass2$average_mass)),
       esize=10, theme="gray")



### independent optimal transport

opttrans = lapply((1:ncol(sigs)), function(i){
  try(transport(a = sigs[,i], b = .sigs_mod[,i],
                     costm = matrix(0,7,7)))
})
names(opttrans) = colnames(sigs)
opttrans = opttrans[sapply(opttrans, typeof) == "list"]

opttrans_melt = (melt(opttrans, id.vars=c('mass', 'to', 'from')))

ggplot(opttrans_melt)+#geom_segment(aes(x=from,y=from,xend=to, yend=to, size=mass))
  geom_jitter(aes(y=from,x=to,size=mass), alpha=0.2, col='#5334e0')+theme_bw()+guides(size=FALSE)
ggsave("../../results/exploratory/bleeding_signature_s5/independent_optimaltransport.pdf", width = 5, height = 5)
