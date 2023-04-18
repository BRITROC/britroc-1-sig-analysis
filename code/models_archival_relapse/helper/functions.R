alphatrans = function(compositional_vec, alpha){
  # stopifnot(all(rowSums(compositional_vec) == 1))
  if(! (sum(compositional_vec) == 1) ){
    compositional_vec = compositional_vec/sum(compositional_vec)
  }
  
  D = length(compositional_vec)
  ua = (compositional_vec**alpha)
  ua = matrix(ua/sum(ua))
  H = Compositional::helm(D)
  return(H %*% (D * (ua-1) / alpha))
}

inverse_alpha_trans = function(transformed_vec, alpha){
  D = length(transformed_vec) + 1
  H = Compositional::helm(D)
  .part = (alpha*t(H) %*% transformed_vec + 1)
  .part = .part**(1/alpha)
  .part = .part/sum(.part)
  return(.part)
}

add_diag <- function(m, d){
  stopifnot(ncol(m) == nrow(m))
  stopifnot(ncol(m) == length(d))
  diag(m) <- d
  m
}

give_LRchanges_barplot = function(obj_it, nrow_facets=NULL, return_df=F){
  ## find the matches. It can be that there are multiple archival, or multiple relapse, or both
  
  ## we're only going to keep samples in which we have at least one sample in each group
  mat_LRchanges0 = sapply(1:obj_it$num_individuals, function(l){
    # patient.meta[patient.meta$PATIENT_ID == unique(patient.meta$PATIENT_ID)[i],]
    subset_meta = patient.meta[patient.meta$PATIENT_ID == (unique(patient.meta$PATIENT_ID))[l],]
    if(all(c('arx', 'rlps') %in% subset_meta$group)){
      ## if at least one from either group
      ## do all possible combinations of log-ratios between arx and rlps
      x = do.call('rbind', lapply(as.character(subset_meta[subset_meta$group == 'arx','SAMPLE_ID']), function(i){
        do.call('rbind', lapply(as.character(subset_meta[subset_meta$group == 'rlps','SAMPLE_ID']), function(j){
          log (obj_it$Y[j,] / obj_it$Y[i,])
        }))
      }))
      rownames(x) = paste0(unique(patient.meta$PATIENT_ID)[l], '_', 1:nrow(x))
      return(x)
    }
  })
  mat_LRchanges0 = do.call('rbind', mat_LRchanges0)
  mat_LRchanges = melt(mat_LRchanges0)
  mat_LRchanges$col = 'NonInf'
  mat_LRchanges$col[is.infinite(mat_LRchanges$value)] = 'Inf'
  mat_LRchanges$value[is.infinite(mat_LRchanges$value)] = min( remove_inf(mat_LRchanges$value))
  mat_LRchanges$col[is.nan(mat_LRchanges$value)] = 'Inf'
  mat_LRchanges$value[is.nan(mat_LRchanges$value)] = max( remove_inf(mat_LRchanges$value))
  
  if(return_df){
    return(value)
  }else{
    a = ggplot(mat_LRchanges,
               aes(x=interaction(Var1, Var2), y=value, fill=col))+geom_bar(stat = 'identity')+
      labs(x='Archival/relapse combination', y='logR change')+
      scale_fill_manual(values=c('Inf'='blue', 'NonInf'='black'))+
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(), legend.title = element_blank())+
      ggtitle(paste0('Log ratio changes'))
    if(is.null(nrow_facets)){
      a = a + facet_wrap(.~Var2, scales = "free_x")
    }else{
      a = a + facet_wrap(.~Var2, scales = "free_x", nrow=nrow_facets)
    }
    return(a)
  }
}

give_paired_subtractions = function(patients_list, exposures, nrow_facets=NULL){
  ## find the matches. It can be that there are multiple archival, or multiple relapse, or both
  
  ## we're only going to keep samples in which we have at least one sample in each group
  mat_LRchanges0 = sapply(patients_list, function(l){
    # patient.meta[patient.meta$PATIENT_ID == unique(patient.meta$PATIENT_ID)[i],]
    subset_meta = patient.meta[patient.meta$PATIENT_ID == l,]
    if(all(c('arx', 'rlps') %in% subset_meta$group)){
      ## if at least one from either group
      ## do all possible combinations of log-ratios between arx and rlps
      x = do.call('rbind', lapply(as.character(subset_meta[subset_meta$group == 'arx','SAMPLE_ID']), function(i){
        do.call('rbind', lapply(as.character(subset_meta[subset_meta$group == 'rlps','SAMPLE_ID']), function(j){
           (exposures[j,] - exposures[i,])
        }))
      }))
      rownames(x) = paste0(unique(patient.meta$PATIENT_ID)[l], '_', 1:nrow(x))
      return(x)
    }
  })
  mat_LRchanges0 = do.call('rbind', mat_LRchanges0)
  # mat_LRchanges = melt(mat_LRchanges0)
  # mat_LRchanges$col = 'NonInf'
  # mat_LRchanges$col[is.infinite(mat_LRchanges$value)] = 'Inf'
  # mat_LRchanges$value[is.infinite(mat_LRchanges$value)] = min( remove_inf(mat_LRchanges$value))
  # mat_LRchanges$col[is.nan(mat_LRchanges$value)] = 'Inf'
  # mat_LRchanges$value[is.nan(mat_LRchanges$value)] = max( remove_inf(mat_LRchanges$value))

  return(mat_LRchanges0)
}

give_short_names = function(i){
  if(grepl('TCGA', i)){
    paste0(strsplit(i, split = '-')[[1]][1:3], collapse='-')
  }else{
    i
  }
}

give_z_matrix_from_labels = function(lbls){
  t(apply(sapply(lbls, function(i) i == unique(lbls)), 2, as.numeric))
}

dropcolnames = function(df){
  colnames(df) = NULL
  df
}

remove_inf = function(i){
  i = i[!is.infinite(i)]
  i[!is.nan(i)]
}

quantifySignatures_alt = function(sample_by_component,component_by_signature=NULL, subsample_sigs=NULL)
{
  if(is.null(component_by_signature))
  {
    component_by_signature<-readRDS(paste(this_path,"data/feat_sig_mat.rds",sep="/"))
  }
  if(!is.null(subsample_sigs)){
    component_by_signature = component_by_signature[,subsample_sigs]
  }
  signature_by_sample<-YAPSA::LCD(t(sample_by_component),
                                  YAPSA:::normalize_df_per_dim(component_by_signature,2))
  signature_by_sample<-normaliseMatrix(signature_by_sample, sig_thresh = 0)
  signature_by_sample
}

give_df_zeros_given_thresh = function(threshold, perturbed_sc_sigs, sigs_from_phil){
  zeros_simulation = do.call('rbind', lapply(perturbed_sc_sigs, function(k) apply(k, 2, median) < threshold))
  ## some super naive approach: if the median is lower than 0.05, then say it's a zero
  
  # table(phil=as.vector(t(sigs_from_phil) == 0), mine=as.vector(zeros_simulation))
  
  zeros_df = cbind.data.frame(phil_zeros=as.vector(t(sigs_from_phil) == 0),
                              my_zeros=as.vector(zeros_simulation),
                              sig=rep(paste0('s', 1:7), each=nrow(zeros_simulation)),
                              sample=colnames(sigs_from_phil))
  
  ## with my simulation we gain zeros of s1, s2, and s4, and gain zeros of s7
  ggplot((melt(zeros_df, id.vars = c('sig', 'sample'))), aes(x=variable, y=value, group=sample))+geom_line()+facet_wrap(.~sig, nrow=1)
  
  zeros_df2 = zeros_df; zeros_df2$transition=paste0(zeros_df2$phil_zeros, '_', zeros_df2$my_zeros)
  head(zeros_df2)
  modify_names = function(i){
    sapply(i, function(j){
      if(j == 'TRUE_TRUE'){
        'Consistent zero'
      }else if(j == 'FALSE_FALSE'){
        'Consistent nonzero'
      }else if(j == 'TRUE_FALSE'){
        'Zero in original only' 
      }else if(j == 'FALSE_TRUE'){
        'Zero in simulation only' 
      }
    })
  }
  zeros_df2$transition = modify_names(zeros_df2$transition)
  zeros_df2
}


give_partial_ilr_basis = function(which_zero_vector, d, previous_implementation=F){
  warning('Warning - this gave an ILR basis which is not orthogonal, but not orthonormal. 20211215: changed to that it is')
  if(length(which_zero_vector) == 0){
    ## no zeros
    return(compositions::ilrBase(D = d))
  }else{
    if(previous_implementation){
      irl_base_complete = compositions::ilrBase(D = d)
      irl_base_partial1 = irl_base_complete
      for(idx_missing in which_zero_vector){
        irl_base_partial1[,idx_missing-1] = 0
        irl_base_partial1[idx_missing,] = 0
        if(idx_missing <= ncol(irl_base_partial1) ){
          for(i in idx_missing:ncol(irl_base_partial1)){
            above_vals = irl_base_partial1[(1:(i))[-which_zero_vector],i]
            ## substitute above vals
            irl_base_partial1[(1:(i))[-which_zero_vector],i] = -irl_base_partial1[(i+1),i]/(length(above_vals))
          }
        }
      }
      return(irl_base_partial1)
    }else{
      ## new, to preserve orthonormality (20211217)
      if(d-length(which_zero_vector) == 1){
        return(rep(0, d))
      }else{
        irl_base_complete_sub <- compositions::ilrBase(D = d-length(which_zero_vector))
        
        irl_base_partial1 <- matrix(NA, ncol=d-1, nrow=d)
        irl_base_partial1[which_zero_vector,] = 0
        irl_base_complete_sub
        for(i in which_zero_vector){
          irl_base_partial1[,max(1, i-1)] = 0 ## max for exception of column 1 when i==1
        }
        irl_base_partial1[is.na(irl_base_partial1)] <-  irl_base_complete_sub
        return(irl_base_partial1)
      }
    }
  }
}



give_partial_irl = function(arg_exposures, ...){
  d = ncol(arg_exposures)
  t(apply(arg_exposures, 1, function(rw){
    .which_zero = as.vector(which(rw==0))
    if(length(.which_zero) >= (d-1)){
      rep(0, d-1)
    } else{
      .basis_row = give_partial_ilr_basis(.which_zero, d, ...)
      compositions::clr(rw) %*% .basis_row
    }
  }))
}


give_short_names = function(i){
  if(grepl('TCGA', i)){
    paste0(strsplit(i, split = '-')[[1]][1:3], collapse='-')
  }else{
    i
  }
}

dropcolnames = function(df){
  colnames(df) = NULL
  df
}

give_short_names = function(i){
  if(grepl('TCGA', i)){
    paste0(strsplit(i, split = '-')[[1]][1:3], collapse='-')
  }else{
    i
  }
}

give_TMB_params = function(arg_d, arg_num_individuals){
  list(beta = (matrix(runif(arg_d*2, min = -4, max = 4),
                      nrow = 2, byrow=TRUE)),
       u_large = matrix(rep(1, (arg_d)*(arg_num_individuals)), nrow=arg_num_individuals),
       logs_sd_RE=runif(n = arg_d, min = 0, max = 2),
       cov_RE = runif(n = ((arg_d)*(arg_d)-(arg_d))/2, min = 0.1, max = 0.2))
}

remove_only_zero_rows <- function(m){
  m[rowSums(m) > 0,]
}


prepare_TMB_data_with_subset = function(subset_sigs, exposures=exposures, .keep_additional=T, add_x_vec=NULL, ilr_trans=T, ...){
  warning('If you get <Error in normalise_rw(exposures[, subset_sigs]) : 
  promise already under evaluation: recursive default argument reference or earlier problems?> error, add <exposures=exposures>')
  .exps = normalise_rw(exposures[,subset_sigs])
  if(ilr_trans){
    .irl_with_zeros = give_partial_irl(.exps, ...)
    .keep = (rowSums(.irl_with_zeros == 0) < (ncol(.exps) - 2) ) & .keep_additional
  }else{
    .irl_with_zeros <- .exps
    .keep <- rep(T, nrow(exposures)) & .keep_additional
  }
  .z_britroc = give_z_matrix_from_labels(patient.meta$PATIENT_ID[.keep])
  rownames(.z_britroc) <- patient.meta$SAMPLE_ID[.keep] ## 20220321
  colnames(.z_britroc) <- unique(patient.meta$PATIENT_ID[.keep]) ## 20220321
  .irl_with_zeros = .irl_with_zeros[.keep,]  ## if a sample has too many zeros (i.e. all but 1 or 2) remove
  
  if(!is.null(add_x_vec)){
    add_x = cbind(1, add_x_vec[.keep])
  }else{
    add_x= cbind(1, as.numeric(as.factor(patient.meta$group[.keep]))-1)
    }
  TMB_data_sim = list(Y = .irl_with_zeros,
                      num_individuals = ncol(.z_britroc),
                      d = length(subset_sigs)-1,
                      n = nrow(.irl_with_zeros),
                      x = add_x,
                      z = .z_britroc)
  return(TMB_data_sim)
}

give_summary_dimensions_TMBobj <- function(TMB_data){
  cat('Dimension of Y: ', dim(TMB_data$Y), '; Dimension of z: ',
      dim(TMB_data$z),
      '\nNumber of dimensions: ',  TMB_data$d, '; Number of samples: ',
      TMB_data$n,  '; Num of individuals: ', TMB_data$num_individuals)
}

compare_matrices <- function(mat1, mat2, remove_zeros=F, facets=F, groups=NULL, title=NULL, ncol=NULL){
  require(reshape2)
  require(ggplot2)
  .x <- melt(list(mat1,mat2))
  .x <- dcast(.x, Var1+Var2~L1, value.var = "value")
  if(!is.null(groups)){
    .x$group <- factor(groups[factor(.x$Var1, levels=unique(.x$Var1))])
  }else{
    .x$group <- '1'
  }
  if(remove_zeros){
    if(sum(.x$`1` == 0) > 0) .x[.x$`1` == 0,]$`1` <- NA
    if(sum(.x$`2` == 0) > 0) .x[.x$`2` == 0,]$`2` <- NA
  }
  
  if(facets){
    plt <- ggplot(.x, aes(x=`1`, y=`2`, shape=group, col=factor(Var2)))+geom_abline(slope = 1, intercept = 0, lty='dashed')+
      geom_point()+theme_bw()+labs(col='Category')
    if(is.null(ncol)){
      plt <- plt+facet_wrap(.~Var2)
    }else{
      plt <- plt+facet_wrap(.~Var2, ncol=ncol)
    }
  }else{
    plt <- ggplot(.x, aes(x=`1`, y=`2`, shape=group, col=factor(Var2)))+geom_abline(slope = 1, intercept = 0, lty='dashed')+
      geom_point()+theme_bw()+labs(col='Category')
  }
  if(!is.null(title)){
    plt <- plt+ggtitle(title)
  }
  
  plt
}

add_colnames <- function(a, colnames_arg){colnames(a) <- colnames_arg; a}
add_rownames <- function(a, rownames_arg){rownames(a) <- rownames_arg; a}

wrapper_run_TMB_use_nlminb <-  function(model, object=NULL, use_nlminb=F, initial_params=NULL, return_all=F, iter.max=150, basic_way_of_creating_groups=NULL){

  ## if the object of data and covariates is an argument
  data = object
  
  stopifnot(nrow(data$x) == nrow(data$Y)) ## 20220328
  # data$Y = matrix(data$Y, nrow=nrow(data$Y)) ## replaced by below 20220301
  data$Y = as(data$Y, 'matrix')
  data$x = (matrix(data$x, ncol=ncol(data$x)))
  print(c('Dim data$Y: ', dim(data$Y),
        'Dim data$x: ', dim(data$x)))
  
  d <- ncol(data$Y) ## number of signatures
  n <- nrow(data$z) ## number of samples
  n_individuals <- ncol(data$z)
  print(c('d: ', d,
        'n: ', n,
        'n_individuals:', n_individuals))
  if(is.null(n)){
    warning('No z matrix. Setting n=nrow(data$Y)')
    n=nrow(data$Y)
  }

  if(ncol(data$x) != 2){
    ## not the standard single-covariate analysis - make beta have as many rows as necessary
    beta_init = (matrix(rep(runif(1, min = -4, max = 4), ncol(data$x)*d),
                        nrow = ncol(data$x), byrow=TRUE))
  }else{
    beta_init = (matrix(rep(runif(1, min = -4, max = 4), 2*(d)),
                        nrow = 2, byrow=TRUE))
  }

  if(is.null(initial_params)){
    if(model %in% c("tmb_MVN_partial_ILR_FEb")){
      initial_params <- list(
        beta = beta_init,
        logsd=rep(1, d)
      )
    }else{
      warning('Change: it used to be       \n    << u_large = matrix(rep(1, (d)*n), nrow=n) >> , \nand now it is  
      << u_large = matrix(rep(1, (d)*n_individuals), nrow=n_individuals) >>,  \n (20220131)')
      initial_params <- list(
        beta = beta_init,
        u_large = matrix(rep(1, (d)*n_individuals), nrow=n_individuals),
        # u_large = matrix(rep(1, (d)*n), nrow=n),
        logs_sd_RE=rep(1, d),
        cov_RE = rep(1, ((d)*(d)-(d))/2)
      )
    }
  }

  if(model == "tmb_MVN_partial_ILR_noularge"){
    initial_params$u_large <- NULL
    dll_name <- "tmb_MVN_partial_ILR_noularge"
    rdm_vec <- NULL
  }else if(model == "tmb_MVN_partial_ILR"){
    stop('Use <model=partialILR>') ## 20220131
      dll_name <- "tmb_MVN_partial_ILR"
      rdm_vec <- "u_large"
  }else if(model == "partialILR"){
    dll_name <- "tmb_MVN_partial_ILR"
    rdm_vec <- "u_large"
    data$num_individuals <- n_individuals
  }else if(model == "partialILRnocor"){
    dll_name <- "tmb_MVN_partial_ILR_notcor"
    rdm_vec <- "u_large"
    data$num_individuals <- n_individuals
    initial_params$cov_RE <- NULL ## no correlations
  }else if(model == "partialILRnocorbetainterceptinRE"){
    stop('Go to  <partialILRnocoroutsidesd> instead')
    dll_name <- "tmb_MVN_partial_ILR_notcor_betainterceptinRE"
    # initial_params$u_large <- NULL
    # rdm_vec <- NULL
    # initial_params$logsd=rep(1, d)
    # initial_params$logs_sd_RE =  NULL
    initial_params$cov_RE =  NULL
    
    rdm_vec <- "u_large"
    # rdm_vec <-NULL
    data$num_individuals <- n_individuals
    #     # data$x = matrix(data$x[,2])
    #     # initial_params$betai = t(matrix(initial_params$beta[1,]))
    #     # initial_params$betas = t(matrix(initial_params$beta[2,]))
    #     # initial_params$beta <- NULL
    # initial_params$cov_RE <- NULL ## no correlations
    # print(initial_params)
    # print(data)
    # pdf("~/Desktop/z.pdf")
    # image(data$z)
    # dev.off()
  }else if(model == "partialILRnocoroutsidesd"){
    dll_name <- "tmb_MVN_partial_ILR_notcor_outsidesd"
    initial_params$cov_RE =  NULL
    rdm_vec <- "u_large"
    data$num_individuals <- n_individuals
  }else if(model == "partialILRoutsidesd"){
    dll_name <- "tmb_MVN_partial_ILR_outsidesd"
    rdm_vec <- "u_large"
    data$num_individuals <- n_individuals
  }else if(model == "partialILRnocoroverdisp"){
    dll_name <- "tmb_MVN_partial_ILR_notcor_overdisp"
    rdm_vec <- "u_large"
    data$num_individuals <- n_individuals
    initial_params$cov_RE <- NULL ## no correlations
    initial_params$lambda <- 1
  }else if(model == "tmb_MVN_partial_ILR_FEb"){
    initial_params$u_large <- NULL
    dll_name <- "tmb_MVN_partial_ILR_FEb"
    rdm_vec <- NULL
  }else if(model == "tmb_MVN_partial_ILR_FEc"){
    dll_name <- "tmb_MVN_partial_ILR_FEc"
    initial_params$u_large <- NULL
    rdm_vec <- NULL
  }else if(model == "tmb_MVN_partial_ILR_FEd"){
    dll_name <- "tmb_MVN_partial_ILR_FEd"
    rdm_vec <- NULL
    data$num_individuals <- NULL
  }else if(model == "tmb_MVN_partial_ILR_FEe"){
    dll_name <- "tmb_MVN_partial_ILR_FEe"
    rdm_vec <- NULL
    initial_params$Lcov = initial_params$cov_RE
    initial_params$u_large <- NULL
    initial_params$cov_RE <- NULL
    initial_params$logsd = initial_params$logs_sd_RE
    initial_params$logs_sd_RE <- NULL
    data$num_individuals <- NULL
  }else if(model == "tmb_MVN_partial_ILR_FEd_integrateularge"){
    dll_name <- "tmb_MVN_partial_ILR_FEd"
    rdm_vec <- "u_large"
    data$num_individuals <- NULL
  }else if(model == "mvnbetacor"){
    dll_name <- "tmb_MVN_with_mean_2"
    rdm_vec <- NULL
    data$num_individuals <- NULL
  }else if(model == "mvnbetacorONEGROUP"){
    dll_name <- "tmb_MVN_with_mean_2"
    rdm_vec <- NULL
    data$num_individuals <- NULL
    initial_params$beta <- t(matrix(initial_params$beta[1,]))
  }else if(model == "tmb_MVN_ILR"){
    stop('Use <tmbMVNILR> instead')
    # dll_name <- "tmb_MVN_ILR"
    # rdm_vec <- "u_large"
    data$num_individuals <- n_individuals
  }else if(model == "tmbMVNILR"){ ## repeated
    dll_name <- "tmb_MVN_ILR"
    rdm_vec <- "u_large"
    data$num_individuals <- n_individuals
  }else if(model == "tmb_MVN_ILRnocor"){
    stop("use <tmbMVNILRnocor> instead")
  }else if(model == "tmbMVNILRnocor"){
      dll_name <- "tmb_RE_20220222"
      rdm_vec <- "u_large"
      data$num_individuals <- n_individuals
      initial_params$cov_RE <- NULL ## no correlations}
  }else if(model == "tmbMVNILRnocorv2"){
    dll_name <- "tmb_RE_20220418"
    rdm_vec <- "u_large"
    data$num_individuals <- n_individuals
    initial_params$cov_RE <- NULL ## no correlations}
  }else if(model %in% c("mvn_beta_no_cor", "mvn_beta_no_corRESTRICTEDBETAS", "mvn_beta_no_cor_weightedbybaseline",
                        "mvn_beta_no_cor_nologsd", "mvn_beta_no_corRESTRICTEDBETAS_nologsd")){
    dll_name <- model
    rdm_vec <- NULL
    data$num_individuals <- NULL
    if(!(model %in% c("mvn_beta_no_cor_nologsd", "mvn_beta_no_corRESTRICTEDBETAS_nologsd"))){
      initial_params$logs_sd <- initial_params$logs_sd_RE
    }
    initial_params$logs_sd_RE <- NULL
    initial_params$u_large <- NULL
    initial_params$cov_RE <- NULL
    if(model %in% c('mvn_beta_no_corRESTRICTEDBETAS', 'mvn_beta_no_corRESTRICTEDBETAS_nologsd')){
      initial_params$beta <- initial_params$beta[1,]
    }
  }else if(model == 'tmb_MVN_with_mean_2_scaledsdpergroup'){
    dll_name <- model
    rdm_vec <- NULL
    data$num_individuals <- NULL
    data$lambda_accessory_mat = lambda_accessory_from_x(data$x)
    initial_params$log_lambda <- 1
  }else if(model == 'bernoulliFE'){
    dll_name <- "tmb_correlated_multinom_2_allFE_b"
    rdm_vec <- NULL
  }else if(model == 'bernoulliME'){
    dll_name <- "tmb_correlated_multinom_2"
    rdm_vec <- 'u_large'
  }else if(model == 'bernoulliMEnocor'){
    dll_name <- "tmb_correlated_multinom"
    rdm_vec <- "u_large"
    data$num_individuals <- n_individuals
  }else if(model == 'partialILRmatchedsubtract'){
    ## change the number of samples
    ## from z and x, get all the pairs of samples in each of the groups
    ## for each sample in the first group, get its equivalent/s in the second group
    data$z
    data$x[,2]
    # basic_way_of_creating_groups=T
    if(basic_way_of_creating_groups){
      cat('Creating groups using basic way of assuming samples are paired as 1, 2, 3, ..., 1, 2, 3, ...\n')
      data$group_1_idx <- ( 1:(n/2) - 1) ## zero-indexed
      data$group_2_idx <- ( (n/2+1):n ) - 1 ## zero-indexed
      data$x <- matrix(data$x[,1], ncol=1)
      data$num_individuals <- n/2
    }else{
      cat('Creating groups based on Z\n')
      indices_pairs <- do.call('rbind', lapply(which(data$x[,2] == 0), function(idx_first_group) rbind(which(data$z[,idx_first_group] == 1), idx_first_group)))
      indices_pairs <- t(indices_pairs[!(indices_pairs[,1] == indices_pairs[,2]),])
      print(indices_pairs)
      data$Y <- data$Y[c(indices_pairs[1,], indices_pairs[2,]),]
      data$group_1_idx <- indices_pairs[1,] - 1 ## zero-indexed
      data$group_2_idx <- indices_pairs[2,] - 1 ## zero-indexed
      data$x <- matrix(rep(1, nrow(data$Y)), ncol=1)
      data$num_individuals <- ncol(indices_pairs)
      print(dim(data$Y))
      print(dim(data$x))
      print(data$num_individuals)
      print(length(data$group_1_idx))
      print(length(data$group_2_idx))
    }
    dll_name <- "tmb_MVN_partial_ILR_matchedsubtractB"
    rdm_vec <- NULL
    data$z <- NULL
    initial_params$cov_RE <- NULL
    initial_params$beta <- matrix(beta_init[2,], nrow=1)
    print(dim(initial_params$beta))
    initial_params$u_large <- NULL
  }else{
    stop('Specify correct <model>\n')
  }
  
  if(!is.null(initial_params)){
    ## initial parameters are passed as arguments
    parameters <- initial_params
  }
  
  if(is.null(rdm_vec)){
    cat('Random vector is null: fixed-effects\n')
    ## fixed effects model
    obj <- MakeADFun(data, parameters, DLL=dll_name)
  }else{
    cat('Random vector is not null: random-effects\n')
    ## random effects model
    obj <- MakeADFun(data, parameters, DLL=dll_name, random = rdm_vec)
  }
  
  if(use_nlminb){
    opt = nlminb(start = obj$par, obj = obj$fn, gr = obj$gr, iter.max=iter.max)
  }else{
    obj$hessian <- TRUE
    opt <- do.call("optim", obj)
    opt
    opt$hessian ## <-- FD hessian from optim
  }

  return_report <- sdreport(obj)
  
  ## give_pairs_with_mvn_wrapper(matrix(return_report$par.random, ncol=object$d), common_lims = T) ## very correlated

  # if(sort_columns){
  #   ## return results in the original order
  #   order_cats
  #   return_report$par.fixed[grepl('beta', names(return_report$par.fixed))] = as.vector(matrix(python_like_select_name(return_report$par.fixed, 'beta'), nrow=2)[,order_cats])
  #   return_report$cov.fixed
  #   return_report
  # }
  
  if(return_all){
    return(list(sdreport=return_report, obj=obj, opt=opt))
  }else{
    return(return_report)
  }
}


imputate_TMB <- function(i, impute) {
  i$Y <- impute(i$Y, impute)
  i
}

give_ALR <- function(i) {
  i$Y <- as(compositions::alr(i$Y), 'matrix')
  i
}

lambda_accessory_from_x <- function(x){
  t(sapply(x[,2], function(i){.x <- rep(0, 2); .x[i+1] <- 1; .x}))
}

compare_TMB_fit_to_data <- function(rep_nlminb, TMB_data, ...){
  RE_res_nlminb <- matrix(rep_nlminb$par.random, ncol=TMB_data$d)
  fitted_res_nlminb <- TMB_data$x %*% matrix(python_like_select_name(rep_nlminb$par.fixed, 'beta'), nrow=2)+
    TMB_data$z %*% RE_res_nlminb
  compare_matrices(mat1 = TMB_data$Y,
                   mat2 = add_rownames(add_colnames(fitted_res_nlminb, colnames(TMB_data$Y)),
                                       rownames(TMB_data$Y)), ...) ## removing zero entries
}

compare_TMB_fits <- function(rep_nlminb, rep_nlminb2, TMB_data, ...){
  RE_res_nlminb <- matrix(rep_nlminb$par.random, ncol=TMB_data$d)
  RE_res_nlminb2 <- matrix(rep_nlminb2$par.random, ncol=TMB_data$d)
  fitted_res_nlminb <- TMB_data$x %*% matrix(python_like_select_name(rep_nlminb$par.fixed, 'beta'), nrow=2)+
    TMB_data$z %*% RE_res_nlminb
  fitted_res_nlminb2 <- TMB_data$x %*% matrix(python_like_select_name(rep_nlminb2$par.fixed, 'beta'), nrow=2)+
    TMB_data$z %*% RE_res_nlminb2
  
  compare_matrices(mat1 = add_rownames(add_colnames(fitted_res_nlminb2, colnames(TMB_data$Y)),
                                       rownames(TMB_data$Y)),
                   mat2 = add_rownames(add_colnames(fitted_res_nlminb, colnames(TMB_data$Y)),
                                       rownames(TMB_data$Y)), ...) ## removing zero entries
}

plot_id <- function(i,j,...){plot(i,j);abline(coef = c(0,1), lty='dashed', col='blue')}

rearrange_cols <- function(matrix, last_col){
  cbind(matrix[,-last_col], matrix[,last_col])
}


load_britroc_exposures <- function(path="../../../britroc-1/data/britroc_30kb_signature_data.rds"){
  ## avoid everything from this file from loading (it's actually a RData file, not RDS)
  load(path)
  return(t(sig_quants))
}

load_britroc_meta <- function(path="../../../britroc-1/data/britroc_30kb_signature_data.rds"){
  ## avoid everything from this file from loading (it's actually a RData file, not RDS)
  load(path)
  return(patient.meta)
}

give_rownames <- function(i, rw){
  rownames(i) <- rw; i}

give_colnames <- function(i, cl){
  colnames(i) <- cl; i}

give_rowcolnames <- function(i, n){
  colnames(i) <- n; rownames(i) <- n; i}

give_ternary_v3 <- function(TMB_data, exposures, add_par=T, opacity=0.2, col=NULL, pch='.', cex=5, legend_off=F, selected_sigs=c(3,4,5),...){
  TMB_data = give_subset_samples_TMBobj(TMB_data,
                                        samples_to_remove = names(which(rowSums(exposures[,selected_sigs]) == 0)))
  probs <- remove_all_NA(normalise_rw(exposures[,selected_sigs]))
  par(mar=c(0,0,0,0))
  TernaryPlot(atip = colnames(probs)[1], btip = colnames(probs)[2], ctip = colnames(probs)[3],
              grid.lines = 0, grid.col = NULL)
  # dens1 <- TernaryDensity(probs[TMB_data$x[,2] == 0,], resolution = 10L, col='blue')
  # dens2 <- TernaryDensity(probs[TMB_data$x[,2] == 1,], resolution = 10L, col='red')
  
  cls_legend = rbind(viridisLite::viridis(48L, alpha = 0.6),
                     seq(from = 0, to = 47, by=1))
  if(!legend_off){
    legend(x=-0.6,y=1.08,
           fill = cls_legend[1,][c(T,F,F,F,F)],
           legend = round(as.numeric(cls_legend[2,][c(T,F,F,F,F)])/sum(dens['z',]), 2), ncol=5,
           y.intersp=0.8,x.intersp=0.5,text.width=0.1, cex=0.9, bty = "n")
  }
  if(is.null(col)){
    col1 <- '#44d7a8'
    col2='purple'
  } else if(length(col) == 2){
    col1 = col[1]
    col2 = col[2]
  } else {
    stop("wrong cols specification")
  }
  TernaryPoints(probs, col = c(col1,
                               col2)[factor(TMB_data$x[,2])], pch = '.', cex=5)
  # TernaryDensityContour(probs[TMB_data$x[,2] == 0,], resolution = 30L, col=col1)
  # TernaryDensityContour(probs[TMB_data$x[,2] == 1,], resolution = 30L, col=col2)
  for(ncolz in 1:ncol(TMB_data$z)){
    cat(ncolz, '\n')
    TernaryLines(probs[which(TMB_data$z[,ncolz] == 1),], col=alpha('black', 0.2))
  }
}

python_like_select_name = function(vector, grep_substring){
  vector[grepl(pattern = grep_substring, x = names(vector))]
}

python_like_select_colnames = function(matrix, grep_substring){
  matrix[,grepl(pattern = grep_substring, x = colnames(matrix))]
}

python_like_select_rownames = function(matrix, grep_substring){
  matrix[grepl(pattern = grep_substring, x = rownames(matrix)),]
}

normalise_rw <- function(x){
  if(is.null(dim(x))){
    x/sum(x)
  }else{
    ## normalise row-wise
    sweep(x, 1, rowSums(x), '/')
  }
}

normalise_cl <- function(x){
  if(is.null(dim(x))){
    x/sum(x)
  }else{
    ## normalise col-wise
    t(sweep(x, 2, colSums(x), '/'))
  }
}

sim_LNM_FE <- function(TMB_res, covariates){
  ## without zeros
  logR <- (covariates) %*% t(matrix(python_like_select_rownames(summary(TMB_res), 'beta')[,1], ncol=2, byrow=T))
  res <- as(compositions::ilrInv(logR), 'matrix')
  colnames(res) <- paste0('s', 1:7)
  rownames(res) <- make.names(covariates[,2], unique = T)
  res
}

sim_LNM_FE_with_sd <- function(TMB_res, covariates){
  ## without zeros
  logR = matrix(NA, nrow = nrow(covariates), ncol=nrow(python_like_select_rownames(summary(TMB_res), 'beta'))/2)
  for(covariates_it in 1:nrow(covariates)){
    betas <- t(matrix(apply(python_like_select_rownames(summary(TMB_res), 'beta'), 1, 
                            function(i) rnorm(n = 1, mean = i[1], sd = i[2])), ncol=2, byrow=T))
    logR[covariates_it,] <- (covariates[covariates_it,]) %*% betas
  }
  res <- as(compositions::ilrInv(logR), 'matrix')
  colnames(res) <- paste0('s', 1:7)
  rownames(res) <- make.names(covariates[,2], unique = T)
  res
}