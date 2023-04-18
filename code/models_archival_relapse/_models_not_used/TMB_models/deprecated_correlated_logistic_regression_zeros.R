rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source("header.R")

exposures_transformed = t(sapply(as.vector(exposures_zeros), function(i) table(factor(i, levels=c(0,1)))))
d = ncol(exposures_transformed)

V_mat = matrix(0, nrow=dim(x)[2]*ncol(exposures), ncol=ncol(exposures))
V_mat[1:(dim(x)[2]), 1] = 1
for(i in 2:ncol(exposures)) V_mat[((i-1)*dim(x)[2]): ((i)*dim(x)[2]), i] = 1

TMB_data = list(Y = exposures_transformed,
                num_individuals = num_indiv,
                num_sigs = ncol(exposures),
                x = t(sapply(rep(1:ncol(exposures), each=nrow(exposures_transformed)/ncol(exposures)), function(i) table(factor(i, levels=1:ncol(exposures))))),
                z = do.call('rbind', lapply(1:ncol(exposures), function(dummy){
                  sapply(unique(patient.meta$PATIENT_ID), function(i) as.numeric(patient.meta$PATIENT_ID == i))
                })),
                V = V_mat)

TMB_params = list(u_large = matrix(rep(1, (d-1)*(TMB_data$num_individuals), nrow=TMB_data$num_individuals)),
                  shared_intercept = matrix(rep(1, nrow(exposures_transformed))),
                  logs_sd_RE=rep(1, d-1),
                  mu_beta = matrix(rep(1, (TMB_data$num_sigs))),
                  logs_sd_cov_beta=rep(1, TMB_data$num_sigs),
                  logs_var_cov_beta=rep(1, ((TMB_data$num_sigs)*(TMB_data$num_sigs)-(TMB_data$num_sigs))/2)
)

obj <- MakeADFun(data = TMB_data, parameters = TMB_params, DLL="tmb_RE_betasignaturecovmat", random = "u_large")
obj$hessian <- TRUE
opt <- do.call("optim", obj)
opt
opt$hessian ## <-- FD hessian from optim
report = sdreport(obj)


unique(names(opt$par))
python_like_select_name(opt$par, "shared_intercept")
python_like_select_name(opt$par, "logs_sd_RE")




## simulate data
par(mfrow=c(1,2))
cov_mat_fitted = give_UNSTRUCTURED_CORR_t_matrix(vec = python_like_select_name(opt$par, "logs_var_cov_beta"), dim_mat = 7)
image(cov_mat_fitted)
diag(cov_mat_fitted)= python_like_select_name(opt$par, "logs_sd_cov_beta")
image(cov_mat_fitted)

python_like_select_name(opt$par, "mu_beta")


report$par.fixed

LR_fitted = matrix(python_like_select_name(report$par.fixed, "shared_intercept"))+
  TMB_data$z %*% matrix(report$par.random)+TMB_data$x %*% matrix(python_like_select_name(report$par.fixed,  "mu_beta"))
thetas = softmax(cbind(LR_fitted, 0))
image(thetas)

molten_thetas = melt(matrix(thetas, ncol=ncol(exposures)))

ggplot(melt(matrix(LR_fitted, ncol=ncol(exposures))))+
  geom_raster( aes( x = Var2, y = Var1, fill = 1-value )  )+
  scale_fill_gradient(low="blue", high="red")+
  geom_point(data=melt(matrix(TMB_data$Y, ncol=ncol(exposures))) %>% filter(value == 0),
             aes( x = Var2, y = Var1 ), col='white', shape=18, size=.3)+
  labs(x='Signature', y='Sample')+
  ggtitle('Probabilities of zero exposure and observed zero exposures')

ggplot(molten_thetas)+
  geom_raster( aes( x = Var2, y = Var1, fill = 1-value )  )+
  scale_fill_gradient(low="blue", high="red")+
  geom_point(data=melt(matrix(TMB_data$Y, ncol=ncol(exposures))) %>% filter(value == 0),
             aes( x = Var2, y = Var1 ), col='white', shape=18, size=.3)+
  labs(x='Signature', y='Sample')+
  ggtitle('Probabilities of zero exposure and observed zero exposures')



