rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# set.seed(1234)

library(uuid)
library(ggplot2)
library(reshape2)
library(compositions)
library(TMB)
library(gridExtra)
library(Ternary)

source("../../../../../GlobalDA/code/2_inference_TMB/helper_TMB.R")
source("../helper/functions.R")
source("../helper/header.R")
source("../../../../Vias_Brenton/copy_number_analysis_organoids/helper_functions.R")

## simulate data where it is only one signature that changes

d <- 5
n <- 100
cov <- diag(rep(1, d-1))
beta <- matrix(c(runif(d-1), c(0,0,0,1)), byrow=T, nrow=2)

## the third category changes

x <- cbind(1, c(rep(0, n/2), rep(1, n/2)))

W <- x%*%beta
W <- as(compositions::alrInv(W+runif(n*(d-1))), 'matrix')
rowSums(W) ## observed data in the simplex

Wilr <- as(ilr(W), 'matrix')
Wilr

Walr <- as(alr(W), 'matrix')
Walr

image(Wilr)
image(Walr)
ilrBase(D = d)

TMB_data = list(Y = Wilr,
                d = d-1,
                n = n,
                x = x,
                z=give_z_matrix(100))

TMB_data_ALR = list(Y = Walr,
                d = d-1,
                n = n,
                x = x,
                z=give_z_matrix(100))
TMB_data_ALRv2 <- TMB_data_ALR
TMB_data_ALRv2$Y <- TMB_data_ALRv2$Y[,ncol(TMB_data_ALRv2$Y):1]

## the model, in cpp
folder_of_TMB_model <- "../tmb_RE/"
TMB::compile(paste0(folder_of_TMB_model, "tmb_MVN_partial_ILR_FEe.cpp"), "-std=gnu++17")
dyn.load(dynlib(paste0(folder_of_TMB_model, "../tmb_RE/tmb_MVN_partial_ILR_FEe")))
TMB::compile(paste0(folder_of_TMB_model, "tmb_MVN_partial_ILR_FEb.cpp"), "-std=gnu++17")
dyn.load(dynlib(paste0(folder_of_TMB_model, "../tmb_RE/tmb_MVN_partial_ILR_FEb")))

# res_FE <- wrapper_run_TMB_use_nlminb(model = "tmb_MVN_partial_ILR_FEe", object = TMB_data)
res_FE <- wrapper_run_TMB_use_nlminb(model = "tmb_MVN_partial_ILR_FEb", object = TMB_data)
res_FE_alr <- wrapper_run_TMB_use_nlminb(model = "tmb_MVN_partial_ILR_FEb", object = TMB_data_ALR)
res_FE_alrv2 <- wrapper_run_TMB_use_nlminb(model = "tmb_MVN_partial_ILR_FEb", object = TMB_data_ALR)
res_FE

grid.arrange(
plot_betas(res_FE),
plot_betas(res_FE_alr))#,
# plot_betas(res_FE_alrv2))


