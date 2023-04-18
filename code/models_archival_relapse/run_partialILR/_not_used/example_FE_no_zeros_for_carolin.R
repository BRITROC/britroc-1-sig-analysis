
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(compositions)
library(TMB)
library(gridExtra)
library(ggplot2)

#-------------------------------------------------------------------------------------------#
## the model, in cpp
folder_of_TMB_model <- "../tmb_RE/" ## YOU CHANGE THIS! RELATIVE PATH TO TMB MODEL (IN CPP)
TMB::compile(paste0(folder_of_TMB_model, "mvn_beta_no_cor.cpp"), "-std=gnu++17")
dyn.load(dynlib(paste0(folder_of_TMB_model, "../tmb_RE/mvn_beta_no_cor")))
#-------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------#
## Functions

impute = function(mat, inputation_value){
  mat[mat == 0] = inputation_value
  normalise_rw(mat)
}

normalise_rw <- function(x){
  if(is.null(dim(x))){
    x/sum(x)
  }else{
    ## normalise row-wise
    sweep(x, 1, rowSums(x), '/')
  }
}

python_like_select_name = function(vector, grep_substring){
  vector[grepl(pattern = grep_substring, x = names(vector))]
}

plot_betas <- function(TMB_obj, names_cats=NULL, rotate_axis=T, theme_bw=T, remove_SBS=T, only_slope=F, return_df=F, plot=T,
                       line_zero=T, add_confint=F, return_plot=T, return_ggplot=F, title=NULL, add_median=F){
  if(typeof(TMB_obj) == 'character'){
    .summary_betas <- NA
    if(theme_bw){
      plt <- ggplot()+theme_bw()
      if(plot) print(plt)
    }else{
      plt <- ggplot()
      if(plot) print(plt)
    }
  }else{
    .summary_betas <- summary(TMB_obj)
    .summary_betas <- cbind.data.frame(python_like_select_rownames(.summary_betas, 'beta'),
                                       type_beta=rep(c('Intercept', 'Slope')),
                                       LogR=rep(1:(nrow(python_like_select_rownames(.summary_betas, 'beta'))/2), each=2))
    if(only_slope){
      .summary_betas <- .summary_betas[.summary_betas$type_beta == 'Slope',]
    }
    
    if(!is.null(names_cats)){
      if(remove_SBS){
        names_cats <- gsub("SBS", "", names_cats) 
      }
      if(length(unique(.summary_betas$LogR)) != length(names_cats)){
        stop('Number of beta slope/intercept pairs should be the same as the length of the name of the categories')
      }
      .summary_betas$LogR = names_cats[.summary_betas$LogR]
    }
    plt <- ggplot(.summary_betas, aes(x=LogR, y=`Estimate`))
    
    if(line_zero) plt <- plt + geom_hline(yintercept = 0, lty='dashed', col='blue')
    if(add_median) { plt <- plt +
      geom_hline(yintercept = median(c(0,.summary_betas$Estimate[.summary_betas$type_beta == 'Slope'])),
                 lty='dashed', col='red') }
    
    plt <- plt +
      geom_point()+
      geom_errorbar(aes(ymin=`Estimate`-`Std. Error`, ymax=`Estimate`+`Std. Error`), width=.1)+
      ggtitle('Slopes')+facet_wrap(.~type_beta, scales = "free")
    
    if(theme_bw){
      plt <- plt + theme_bw()
    }
    
    if(add_confint){
      confints <- cbind(.summary_betas, confint=t(give_confidence_interval(.summary_betas[,'Estimate'], .summary_betas[,'Std. Error'])))
      plt <- plt+
        geom_errorbar(data = confints, aes(ymin=confint.1, ymax=confint.2), width=.1 ,col='blue', alpha=0.6)
      
    }
    
    if(rotate_axis){
      plt <- plt + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    }
    
    if(!is.null(title)){
      plt <- plt + ggtitle(title)
    }
    
    if(!TMB_obj$pdHess){
      plt <- plt + annotate("text", x = -Inf, y=Inf, label="not PD", vjust=1, hjust=-.2)+geom_point(col='red')
      if(plot) print(plt)
    }else{
      if(plot) print(plt)
    }
  }
  
  if(return_df){
    .summary_betas
  }else{
    if(return_plot & return_df){stop('<return_plot=T> and <return_df=T> are incompatible')}
    plot_list <- list(plt)
    class(plot_list) <- c("quiet_list", class(plot_list))
    if(return_plot){
      return(cowplot::as_grob(plt))
    }else if(return_ggplot){
      return(plt)
    }
  }
}

softmax = function(x){
  if(is.null(dim(x))){
    ## vector
    sum_x = sum(exp(x))
    exp(x)/sum_x
  }else{
    ## matrix
    sum_x = rowSums(exp(x))
    sweep(exp(x), 1, sum_x, '/')
  }
}

python_like_select_rownames = function(matrix, grep_substring){
  matrix[grepl(pattern = grep_substring, x = rownames(matrix)),]
}

select_slope_2 = function(i, verbatim=TRUE){
  if(is.null(dim(i))){
    i[c(F,T)]
  }else{
    i[,c(F,T)]
  }
}

#-------------------------------------------------------------------------------------------#

## Creating data. We'll make only one signature change in absolute terms
n <- 200
d <- 7

## six signatures have the same abundance in the two groups
absolute_exposures <- sapply(c(5, 6, 2, 10, 4, 4), function(lambda_it) rpois(n, lambda_it))
## the signature that we will put in the first row has different abundance between the two groups
absolute_exposures <- cbind(as.vector(sapply(c(10, 3), function(lambda_it) rpois(n/2, lambda_it))),
                            absolute_exposures)
image(t(absolute_exposures)) ## absolute exposures. Very clearly, s1 differs between the groups

exposures <- sweep(absolute_exposures, 1, rowSums(absolute_exposures), '/')
colnames(exposures) <- paste0('s', 1:d)
rownames(exposures) <- paste0('sample_', 1:n)
image(t(exposures)) ## normalised, relative, exposures

## Add imputation, if needed, Alternative ways of imputating data: zCompositions::multLN(), and others
exposures <- impute(exposures, 1e-2)
which_zero = t(apply(exposures, 1, function(i) as.numeric((i==0)) ))

ALR_transformed_data = as(compositions::alr(exposures), 'matrix')
image(t(ALR_transformed_data)) ## absolute exposures

dmin1 = d-1

#-------------------------------------------------------------------------------------------#
TMB_data = list(Y = ALR_transformed_data,
                d = dmin1,
                n = n,
                x = cbind(1, rep(c(0,1), each=n/2)))

TMB_params = list(logs_sd=runif(n = dmin1, min = 0, max = 2),
                  beta = (matrix(runif(dmin1*2, min = -4, max = 4),
                                 nrow = 2, byrow=TRUE))
) ## some initial parameters. The parameter that we are most interested in is the second row of beta (the coefficients for the
## change between the two conditions)

#-------------------------------------------------------------------------------------------#
## Run model

obj <- MakeADFun(data = TMB_data, parameters = TMB_params, DLL="mvn_beta_no_cor")
opt <- do.call("optim", obj)
opt
opt$hessian ## <-- FD hessian from optim
rep <- sdreport(obj)
rep

plot_betas(rep) ## one very clear negative coefficient (corresponding to the logR involving, as the numerator,
## the signature that we have simutated to change, s1)


#-------------------------------------------------------------------------------------------#
## If you wanted to run regression with some continuous value, instead of 
TMB_data_continuous = list(Y = ALR_transformed_data,
                           d = dmin1,
                           n = n,
                           x = cbind(1, runif(n)))
# obj <- MakeADFun(data = TMB_data_continuous, parameters = TMB_params, DLL="tmb_MVN_partial_ILR_FEb")
obj <- MakeADFun(data = TMB_data_continuous, parameters = TMB_params, DLL="mvn_beta_no_cor")
opt <- do.call("optim", obj)
opt
opt$hessian ## <-- FD hessian from optim
rep_cont <- sdreport(obj)
rep_cont
plot_betas(rep_cont) ## x values were generated randomly, so there shouldn't be any change
#-------------------------------------------------------------------------------------------#

