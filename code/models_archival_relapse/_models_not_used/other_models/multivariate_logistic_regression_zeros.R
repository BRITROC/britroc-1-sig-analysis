# output <- glm(sta ~ sex + ser, data=icu1.dat, family=binomial)
# logistic.regression.or.ci(output)

# https://datasciencebeginners.com/2018/12/20/multinomial-logistic-regression-using-r/

# library(rattle.data)
# Loading the wine data
# data(wine)
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


source("header.R")
library(dfidx)

# Using sample_frac to create 70 - 30 slipt into test and train
# train <- sample_frac(wine, 0.7)
# sample_id <- as.numeric(rownames(train)) # rownames() returns character so as.numeric
# test <- wine[-sample_id,]
# 
# # Setting the basline 
# train$Type <- relevel(train$Type, ref = "3")
# 
# table(train$Type)
# 
# 
# require(nnet)
# # Training the multinomial model
# multinom.fit <- multinom(Type ~ Alcohol + Color -1, data = train)

# Checking the model
# summary(multinom.fit)

data_multinom = cbind.data.frame(exp_bool=as.vector(exposures_zeros),
sig=rep(colnames(exposures), each=nrow(exposures)), sample=rep(rownames(exposures), ncol(exposures)))
data_multinom$sig = as.numeric(as.factor(data_multinom$sig))
# data_multinom <- data_multinom[data_multinom$exp_bool == 0, c('sig', 'sample')]
data_multinom <- data_multinom[data_multinom$exp_bool == 1, c('sig', 'sample')]
data_multinom[,'group'] = patient.meta$group[match(data_multinom$sample, patient.meta$SAMPLE_ID)]
data_multinom[,'patient'] = patient.meta$PATIENT_ID[match(data_multinom$sample, patient.meta$SAMPLE_ID)]

multinom.fit2 <- nnet::multinom(sig ~ 1, data = data_multinom)
multinom.fit2$fitted.values

multinom.fitFE <- nnet::multinom(sig ~ group, data = data_multinom)
multinom.fitFE$fitted.values
coefficients(multinom.fitFE)
summary(multinom.fitFE)


## ???
# head(Fishing)

## And now add random effects
## I thunk for this I need some other package
colnames(data_multinom) = c('sig', 'sample.sample', 'group.group', 'patient.patient')
Fish <- dfidx(data_multinom,
              # varying = c('sample', 'group', 'patient'),
              varying = 3, 
              shape = "wide", choice = "sig")

multinom.fitFE <- mlogit::mlogit(sig ~ group + (1|patient), data = data_multinom)
multinom.fitFE <- mlogit::mlogit(sig ~ 1, data = data_multinom)

library(TMB)
source("../../../GlobalDA/code/2_inference_TMB/helper_TMB.R")
TMB::compile("../../../GlobalDA/code/2_inference_TMB/mm_multinomial/fullRE_ME_multinomial.cpp")
dyn.load(dynlib("../../../GlobalDA/code/2_inference_TMB/mm_multinomial/fullRE_ME_multinomial"))
TMB::compile("../../../GlobalDA/code/2_inference_TMB/mm_multinomial/diagRE_ME_multinomial.cpp")
dyn.load(dynlib("../../../GlobalDA/code/2_inference_TMB/mm_multinomial/diagRE_ME_multinomial"))


## below: unnecessary, we can keep all the samples
## select onyl those with two samples, one in each group
## patient.meta = droplevels(patient.meta[patient.meta$paired,])
## table(patient.meta$PATIENT_ID)
## if there is more than one sample, choose one of each group at random
## unique(patient.meta$PATIENT_ID)

data <- list()
data$Y = exposures_zeros
data$x = cbind(1, as.numeric(as.factor(patient.meta$group))-1)
data$z = give_z_matrix_from_labels(patient.meta$PATIENT_ID)

d <- ncol(data$Y)
n <- ncol(data$z) ## number of INDIVIDUALS, not samples
data$num_individuals = n

parameters <- list(
  beta = (matrix(rep(runif(1, min = -4, max = 4), 2*(d-1)),
                 nrow = 2, byrow=TRUE)),
  u_large = matrix(rep(1, (d-1)*n), nrow=n),
  logs_sd_RE=rep(1, d-1),
  cov_par_RE = rep(1, ((d-1)*(d-1)-(d-1))/2)
)

# obj <- MakeADFun(data, parameters, DLL="fullRE_ME_multinomial", random = "u_large")
## I think perhaps the full RE multinomial doesn't work because we have too many parameters
## to infer. The diagonal RE do work.
obj <- MakeADFun(data, parameters, DLL="diagRE_ME_multinomial", random = "u_large")
obj$hessian <- TRUE
opt <- do.call("optim", obj)
opt
opt$hessian ## <-- FD hessian from optim
report = sdreport(obj)

select_slope_2(python_like_select_name(report$par.fixed, 'beta'),v=F)

fitted_multinom = simulate_from_M_TMB(report, full_RE = T, x_matrix=data$x, z_matrix=data$z)
image(t(fitted_multinom))
dim(fitted_multinom)[1] == dim(data$z)[1] ## it should be the same

colnames(fitted_multinom) = colnames(data$Y)

ggplot(melt(fitted_multinom))+
  geom_raster( aes( x = Var2, y = Var1, fill = 1-value )  )+
  scale_fill_gradient(low="blue", high="red")+
  geom_point(data=melt((data$Y)) %>% filter(value == 0),
             aes( x = Var2, y = Var1 ), col='white', shape=18, size=.3)+
  labs(x='Signature', y='Sample')+
  ggtitle('Probabilities of zero exposure and observed zero exposures')
ggsave("../results/zeros_modelling/ME_multinomial_simulation_and_true.png",
       width = 6, height = 4)
# ------------------------------------------------------------ #

## I don't think this above is correct at all. We should do an analysis for each
## signature independently. This below
rowSums(fitted_multinom)
## doesn't make any sense. zeros are not competing for each other

give_binom_per_sig = function(signature_idx){
  require(lme4)
  m <- glmer(data$Y[,signature_idx] ~ data$x[,2] + (1 | patient.meta$PATIENT_ID),
             family = binomial, control = glmerControl(optimizer = "bobyqa"),
             nAGQ = 10)
  return(m)
}

give_binom_per_sig_noRE = function(signature_idx){
  m <- glm(data$Y[,signature_idx] ~ data$x[,2] ,
             family = binomial)
  return(m)
}

res_binom_per_sig = sapply(1:7, give_binom_per_sig)
res_binom_per_sig_noRE = lapply(1:7, give_binom_per_sig_noRE)
fitted_multiple_binom = sapply(res_binom_per_sig, function(i) fitted(i))
colnames(fitted_multiple_binom) = colnames(data$Y)

ggplot(melt(fitted_multiple_binom))+
  geom_raster( aes( x = Var2, y = Var1, fill = 1-value )  )+
  scale_fill_gradient(low="blue", high="red")+
  geom_point(data=melt((data$Y)) %>% filter(value == 0),
             aes( x = Var2, y = Var1 ), col='white', shape=18, size=.3)+
  labs(x='Signature', y='Sample')+
  ggtitle('Probabilities of zero exposure and observed zero exposures')
ggsave("../results/zeros_modelling/ME_multiple_binom_simulation_and_true.png",
       width = 6, height = 4)

## data is sorted by group, but I haven't included this information in the
## plots
patient.meta$group

lapply(res_binom_per_sig, function(i) drop1(i,test="Chisq"))
lapply(res_binom_per_sig_noRE, function(i) drop1(i,test="Chisq"))

## same but without random effects

####
# Using each of the signatures as random effects

x = t(cbind(1, as.numeric(factor(patient.meta$group))-1))

d = ncol(exposures) ## number of features
n = nrow(exposures) ## number of samples
num_indiv = length(unique(patient.meta$PATIENT_ID))

##' change the layout of the exposures in such a way that we have two columns; one
##' for zeros and one for non-zeros, and all samples and signatures and groups have been
##' concatenated
