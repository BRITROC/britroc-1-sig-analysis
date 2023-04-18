## originally in LN_modelling_zeros_partialILR_TMB.R

dummy_data <- TMB_data
nrow(dummy_data$Y)

dummy_alr <- mvtnorm::rmvnorm(n = nrow(dummy_data$Y), mean = c(1, 0.2, 0.8, 1.3, 0.2, 1.1))
dummy_ilr <- compositions::ilr(softmax(cbind(dummy_alr, 0)))
dummy_partialilr <- give_partial_irl(softmax(cbind(dummy_alr, 0)))

par(mfrow=c(1,1))
plot(as.vector(dummy_ilr), as.vector(dummy_partialilr))

dummy_data$Y <- dummy_partialilr


res_nlminb_dummy <- wrapper_run_TMB_use_nlminb(model = "partialILR", object = dummy_data,
                                               use_nlminb = T, iter.max=1000)
res_nlminb_dummy
plot_betas(res_nlminb_dummy) ## with simulated data it also doesn't converge

give_pairs_with_mvn_wrapper(matrix(res_nlminb_dummy$par.random, ncol=dummy_data$d), common_lims = T) ## very correlated
## but the data were not created as correlated with the ALR!! they were not correlated at all, in fact
pairs(dummy_alr)
pairs(dummy_ilr)
res_nlminb_dummyc <- wrapper_run_TMB_use_nlminb(model = "mvnbetacor", object = dummy_data,
                                                use_nlminb = T, iter.max=1000) ## correlated, no zero, no RE
pheatmap(L_to_cov(python_like_select_name(res_nlminb_dummyc$par.fixed, 'cov_RE'),
                  d =  ncol(dummy_data$Y))) ## no correlations

## adding more dummy data
dummy_data2 <- dummy_data
dummy_data2$Y <- rbind(dummy_data2$Y, dummy_data2$Y)
dummy_data2$x <- rbind(dummy_data2$x, dummy_data2$x)
dummy_data2$z <- rbind(dummy_data2$z, dummy_data2$z)
dummy_data2$n <- dummy_data2$n*2
give_pairs_with_mvn_wrapper(dummy_data2$Y)


## does it not converge because of the number of patients? fewer patients
## we actually want more, and not fewer, patients!
# dummy_data2$z <- give_z_matrix_from_labels(sample(1:4, size = nrow(dummy_data$Y), replace = T))
# dummy_data2$num_individuals <- ncol(dummy_data2$z)
## also doesn't converge
res_nlminb_dummy2 <- wrapper_run_TMB_use_nlminb(model = "partialILR", object = dummy_data2,
                                                use_nlminb = T, iter.max=1000) ## correlated, zeros, RE
res_nlminb_dummy2b <- wrapper_run_TMB_use_nlminb(model = "tmb_MVN_ILR", object = dummy_data2,
                                                 use_nlminb = T, iter.max=1000) ## correlated, no zeros, RE
res_nlminb_dummy2c <- wrapper_run_TMB_use_nlminb(model = "mvnbetacor", object = dummy_data2,
                                                 use_nlminb = T, iter.max=1000) ## correlated, no zero, no RE
give_pairs_with_mvn_wrapper(matrix(res_nlminb_dummy2$par.random, ncol=dummy_data2$d), common_lims = T) ## very correlated
give_pairs_with_mvn_wrapper(matrix(res_nlminb_dummy2b$par.random, ncol=dummy_data2$d), common_lims = T) ## very correlated

grid.arrange(plot_betas(res_nlminb_dummy2), plot_betas(res_nlminb_dummy2b),
             plot_betas(res_nlminb_dummy2c))


table(factor(apply(dummy_data2$z, 1, function(i) which(i == 1)), levels=1:(max(apply(dummy_data2$z, 1, function(i) which(i == 1))))),
      factor(dummy_data2$x[,2] == 1, levels=c(T,F)))

## is it a problem of some patients having many samples?
dummy_data3 <- give_subset_samples_TMBobj(dummy_data2, 1) ## removing one sample to make the number even
dummy_data3$z <- give_z_matrix(n_times_2 = nrow(dummy_data3$Y))
dummy_data3$num_individuals <- ncol(dummy_data3$z)
dummy_data3$x[,2] <- rep(c(0, 1), each=nrow(dummy_data3$Y)/2)
dim(dummy_data3$Y)
dummy_data3$num_individuals
dummy_data3$d
dummy_data3$n
dim(dummy_data3$x)
dim(dummy_data3$z)
pheatmap(cov(dummy_data3$Y))

res_nlminb_dummy3 <- wrapper_run_TMB_use_nlminb(model = "partialILR", object = dummy_data3,
                                                use_nlminb = T, iter.max=1000) ## correlated, zeros, RE
res_nlminb_dummy3

## changes
dummy_data3$Y[dummy_data3$x[,2] == 1,] <- dummy_data3$Y[dummy_data3$x[,2] == 1,][,ncol(dummy_data3$Y):1]
res_nlminb_dummy3 <- wrapper_run_TMB_use_nlminb(model = "partialILR", object = dummy_data3,
                                                use_nlminb = T, iter.max=1000) ## correlated, zeros, RE
plot_betas(res_nlminb_dummy3)


give_pairs_with_mvn_wrapper(matrix(res_nlminb_dummy3$par.random, ncol=dummy_data3$d), common_lims = T) ## very correlated

## is it a problem of th covariance matrix (i.e. no correlations?)

correlated_obs <- cbind(runif(200),runif(200) ); correlated_obs <- cbind(correlated_obs, correlated_obs[,1]+runif(200, max = 0.5), correlated_obs[,2]+runif(200, max = 0.5),
                                                                         correlated_obs[,1]*0.1+correlated_obs[,2]*0.2+runif(200, max = 1), correlated_obs[,1]+runif(200, max = 0.25))
pheatmap(cov(correlated_obs)*40, cluster_rows = F, cluster_cols = F)
dummy_alr <- mvtnorm::rmvnorm(n = nrow(dummy_data3$x), mean = c(1, 0.2, 0.8, 1.3, 0.2, 1.1),
                              sigma = cov(correlated_obs)*40)
pheatmap(cov(correlated_obs)*40, cluster_rows = F, cluster_cols = F)
dummy_ilr <- compositions::ilr(softmax(cbind(dummy_alr, 0)))
dummy_data3$Y <- dummy_ilr
dummy_data3$Y <- dummy_alr

dim(dummy_data3$Y)
dummy_data3$n
dummy_data3$d
dummy_data3$num_individuals
dim(dummy_data3$x)
dim(dummy_data3$z)

dummy_data3$num_individuals <- ncol(dummy_data3$z)

give_pairs_with_mvn_wrapper(dummy_data3$Y)
##  the matrix is either rank-deficient or indefinite, if perfect correlation
res_nlminb_dummy3 <- wrapper_run_TMB_use_nlminb(model = "partialILR", object = dummy_data3,
                                                use_nlminb = T, iter.max=1000) ## correlated, zeros, RE
plot_betas(res_nlminb_dummy3)
res_nlminb_dummy3

give_pairs_with_mvn_wrapper(matrix(res_nlminb_dummy3$par.random, ncol=6)) ## random effects are highly correlated
give_pairs_with_mvn_wrapper(dummy_data3$Y)
give_pairs_with_mvn_wrapper(dummy_data3$Y[dummy_data3$x[,2] == 0,])

## in general, the random intercepts for observed data are not too correlated
give_pairs_with_mvn_wrapper(matrix(rep$par.random, ncol=5), common_lims = T) ## correlated
give_pairs_with_mvn_wrapper(matrix(res_nlminb$par.random, ncol=6)) ## ILR1 and ILR6 quite correlated
give_pairs_with_mvn_wrapper(matrix(res_nlminb_amalgamation1$par.random, ncol=6))
give_pairs_with_mvn_wrapper(matrix(res_nlminb_amalgamation1_nos5$par.random, ncol=6))
give_pairs_with_mvn_wrapper(matrix(res_nlminb0$par.random, ncol=6))

rep_nlminb <- wrapper_run_TMB_use_nlminb(model = "partialILR", object = TMB_data_with_subset,
                                         use_nlminb = T, iter.max=1000)
d_xx <- TMB_data_with_subset$d
plot_betas(rep_nlminb)
dim(TMB_data_with_subset$z) ## 262 158
TMB_data_with_subset$num_individuals
## but random effects only have 132 rows

pdf("~/Desktop/partialilrproblems/scatter_random.pdf")
give_pairs_with_mvn_wrapper(matrix(rep_nlminb$par.random, ncol=d_xx), common_lims = T) ## looks good
dev.off()

grid.arrange(plot_betas(rep), plot_betas(rep_nlminb))
rep_nlminb

### the problem is the extreme values in estimates of covariances
.cov_est <- L_to_cov(python_like_select_name(rep_nlminb$par.fixed, 'cov_RE'),
                     d =  ncol(TMB_data_with_subset$Y))
# diag(.cov_est) <- exp(python_like_select_name(rep_nlminb$par.fixed, 'logs_sd_RE'))**2
# .cov_est <- .cov_est / exp(python_like_select_name(rep_nlminb$par.fixed, 'logs_sd_RE')) %*% t(exp(python_like_select_name(rep_nlminb$par.fixed, 'logs_sd_RE')))
pdf("~/Desktop/partialilrproblems/cov_from_estimates.pdf")
pheatmap(.cov_est, cluster_rows = F, cluster_cols = F)
dev.off()

pdf("~/Desktop/partialilrproblems/cov_from_RE_intercepts.pdf")
pheatmap(cov(matrix(rep_nlminb$par.random, ncol=d_xx)), cluster_rows = F, cluster_cols = F)
dev.off()

pdf("~/Desktop/partialilrproblems/cor_from_RE_intercepts.pdf")
pheatmap(cor(matrix(rep_nlminb$par.random, ncol=d_xx)), cluster_rows = F, cluster_cols = F)
dev.off()

plot(matrix(rep_nlminb$par.random, ncol=d_xx)[,c(1,5)]) ## the very high correlation between ILR1 and ILR5

.cor_est <- outer(1:ncol(.cov_est), 1:ncol(.cov_est), Vectorize(function(i,j) .cov_est[i,j]/(sqrt(diag(.cov_est))[i]*sqrt(diag(.cov_est))[j])))
.cor_est ## there are values more extreme than 1, so there is a problem here

sim_u_large <- mvtnorm::rmvnorm(dummy_data3$num_individuals, mean = c(1, 0.2, 0.8, 1.3, 0.2, 1.1), sigma = diag(1, 6))
sim_logR <- dummy_data3$x %*% matrix(python_like_select_name(res_nlminb_dummy3$par.fixed, 'beta'), nrow=2)+
  dummy_data3$z %*% sim_u_large

dummy_data_sim <- dummy_data3
dummy_data_sim$Y <- sim_logR
rep_nlminb_sim <- wrapper_run_TMB_use_nlminb(model = "partialILR", object = dummy_data_sim,
                                             use_nlminb = T, iter.max=1000)
rep_nlminb_sim
### THIS LOOKS GOOD!!!!
give_pairs_with_mvn_wrapper(matrix(rep_nlminb_sim$par.random, ncol=d_xx), common_lims = T) ## looks good
give_pairs_with_mvn_wrapper(matrix(rep_nlminb$par.random, ncol=d_xx), common_lims = T) ## looks good

image(dummy_data_sim$Y)
image(dummy_data$Y)
image(dummy_data$Y)
