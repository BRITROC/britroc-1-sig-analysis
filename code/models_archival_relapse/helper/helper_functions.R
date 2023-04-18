#' Create the barplot given a matrix of exposures
createBarplot <- function(matrix_exposures, angle_rotation_axis = 0, order_labels=NULL,
                          remove_labels=FALSE, levels_signatures=NULL, includeMelt=NULL, Melt=NULL, verbose=TRUE){
  #' error due to it not being a matrix:  
  #'    No id variables; using all as measure variables
  #'    Rerun with Debug
  #'    Error in `[.data.frame`(.mat, , "Var1") : undefined columns selected \\
  #' (use tomatrix())
  
  
  if(is.null(colnames(matrix_exposures))) stop('columns must have names')
  require(reshape2)
  require(ggplot2)
  library(RColorBrewer)
  if(verbose){
    cat(paste0('Creating plot... it might take some time if the data are large. Number of samples: ', nrow(matrix_exposures), '\n'))
  }
  
  if( (!is.null(order_labels) & typeof(order_labels) == "logical")){if(!order_labels){cat('WARNING: Order labels is either a vector with desired order or NULL, not bool')}}
  
  if(!is.null(levels_signatures)){
    library(ggthemes)
    ggthemes_data$economist
  }else{
    levels_signatures <- colnames(matrix_exposures)
  }
  n <- 60
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
  
  myColors <- col_vector[1:length(levels_signatures)]
  names(myColors) <- levels_signatures
  if(is.null(order_labels)) order_labels = rownames(matrix_exposures)
  if(!is.null(includeMelt)){
    cat('For whatever reason sometimes the melt does not work. Here it is passed as argument.')
    .mat <- Melt
  }else{
    .mat <- melt(matrix_exposures)
  }
  .mat[,'Var1'] <- factor(.mat[,'Var1'], levels=order_labels)
  .mat[,'Var2'] <- factor(.mat[,'Var2'], levels=levels_signatures)
  ###rownames(.mat) <- rownames(matrix_exposures) ### new
  
  if(!remove_labels){
    if(!is.null(levels_signatures)){
      ggplot(.mat, aes(x=Var1, y=value, fill=factor(Var2, levels=levels_signatures)))+
        geom_bar(stat = 'identity')+
        theme(axis.text.x = element_text(angle = angle_rotation_axis, hjust = 1))+
        #theme(axis.title.x=element_blank(),
        #      axis.text.x=element_blank(),
        #      axis.ticks.x=element_blank())+
        guides(fill=guide_legend(title="Signature"))+
        scale_fill_manual(name = "grp",values = myColors)
    }else{
      ggplot(.mat, aes(x=Var1, y=value, fill=factor(Var2, levels=levels_signatures)))+
        geom_bar(stat = 'identity')+
        theme(axis.text.x = element_text(angle = angle_rotation_axis, hjust = 1))+
        #theme(axis.title.x=element_blank(),
        #      axis.text.x=element_blank(),
        #      axis.ticks.x=element_blank())+
        guides(fill=guide_legend(title="Signature"))
    }
  }else{
    if(!is.null(levels_signatures)){
      ggplot(.mat, aes(x=Var1, y=value, fill=factor(Var2, levels=levels_signatures)))+
        geom_bar(stat = 'identity')+
        theme(axis.text.x = element_text(angle = angle_rotation_axis, hjust = 1))+
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank())+
        guides(fill=guide_legend(title="Signature"))+
        scale_fill_manual(name = "grp",values = myColors)
    }else{
      ggplot(.mat, aes(x=Var1, y=value, fill=factor(Var2, levels=levels_signatures)))+
        geom_bar(stat = 'identity')+
        theme(axis.text.x = element_text(angle = angle_rotation_axis, hjust = 1))+
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank())+
        guides(fill=guide_legend(title="Signature"))
    }
  }
}


give_joint_histogram = function(some_list, bins=30, no_colour=FALSE){
  require(reshape2)
  require(ggplot2)
  .x = reshape2::melt(some_list)
  if(no_colour){
    ggplot(.x, aes(x=value, group=L1), alpha=0.2)+geom_histogram(bins = bins)+facet_wrap(.~L1, ncol=1, scales='free_y')
  }else{
    ggplot(.x, aes(x=value, fill=L1, group=L1), alpha=0.2)+geom_histogram(bins = bins)+facet_wrap(.~L1, ncol=1, scales='free_y')
  }
}


extractCopynumberFeatures = function(CN_data, cores = 1)
{
  #get chromosome lengths
  chrlen<-read.table(paste(this_path,"data/hg19.chrom.sizes.txt",sep="/"),sep="\t",stringsAsFactors = F)[1:24,]
  
  #get centromere locations
  gaps<-read.table(paste(this_path,"data/gap_hg19.txt",sep="/"),sep="\t",header=F,stringsAsFactors = F)
  centromeres<-gaps[gaps[,8]=="centromere",]
  
  if(cores > 1) {
    require(foreach)
    doMC::registerDoMC(cores)
    
    temp_list = foreach::foreach(i=1:6) %dopar% {
      if(i == 1){
        list(segsize = getSegsize(CN_data) )
      } else if (i == 2) {
        list(bp10MB = getBPnum(CN_data,chrlen) )
      } else if (i == 3) {
        list(osCN = getOscilation(CN_data,chrlen) )
      } else if (i == 4) {
        list(bpchrarm = getCentromereDistCounts(CN_data,centromeres,chrlen) )
      } else if (i == 5) {
        list(changepoint = getChangepointCN(CN_data) )
      } else {
        list(copynumber = getCN(CN_data) )
      }
      
    }
    unlist( temp_list, recursive = FALSE )
  } else {  
    
    segsize<-getSegsize(CN_data)
    bp10MB<-getBPnum(CN_data,chrlen)
    osCN<-getOscilation(CN_data,chrlen)
    bpchrarm<-getCentromereDistCounts(CN_data,centromeres,chrlen)
    changepoint<-getChangepointCN(CN_data)
    copynumber<-getCN(CN_data)
    
    list(segsize=segsize,bp10MB=bp10MB,osCN=osCN,bpchrarm=bpchrarm,changepoint=changepoint,copynumber=copynumber)
  }
  
}


getSampNames<-function(abs_profiles)
{
  if(class(abs_profiles)=="QDNAseqCopyNumbers")
  {
    samps<-colnames(abs_profiles)
  }
  else
  {
    samps<-names(abs_profiles)
  }
  samps
}

getBPnum2<-function(abs_profiles,chrlen){
  out<-c()
  cln = colnames(abs_profiles[[1]])
  abs_profiles = do.call('rbind', abs_profiles)
  colnames(abs_profiles) = cln
  # samps<-getSampNames(abs_profiles)
  # for(i in samps)
  # {
  segTab = abs_profiles#[[i]]
  chrs<-unique(segTab[,'chromosome'])
  allBPnum<-c()
  for(c in chrs)
  {
    currseg<-segTab[segTab[,'chromosome']==c,]
    intervals<-seq(1,chrlen[chrlen[,1]== c,2]+10000000,10000000)
    if(is.null(nrow(currseg))){
      res <- hist(as.numeric(currseg['end'][-1]),breaks=intervals,plot=FALSE)$counts
    }else{
      res <- hist(as.numeric(currseg[,'end'][-nrow(currseg)]),breaks=intervals,plot=FALSE)$counts
    }
    allBPnum<-c(allBPnum,res)
  }
  out<-rbind(out,cbind(ID=rep(i,length(allBPnum)),value=allBPnum))
  # out
  # }
  # rownames(out)<-NULL
  data.frame(out,stringsAsFactors = F)
}

wrapper_get_bp = function(df){
  for(i in 1:length(df)){
    df[[i]][,'chromosome'] = paste0('chr', df[[i]][,'chromosome'])
  }
  names(df) = paste0('df', 1:length(df))
  getBPnum2(df, chrlen = chrlen)
}

create_distrib_df = function(feature, name_value_col){
  df = list()
  df[['pcawg']] = pcawg_CN_features[[feature]] %>% group_by(ID) %>% mutate(ct1_num = as.numeric(get(name_value_col))) %>% summarize(mean(ct1_num), .groups = 'drop')
  df[['org']] = organoids_CN_features[[feature]] %>% group_by(ID) %>% mutate(ct1_num = as.numeric(get(name_value_col))) %>% summarize(mean(ct1_num), .groups = 'drop')
  df[['tcga']] = tcga_CN_features[[feature]] %>% group_by(ID) %>% mutate(ct1_num = as.numeric(get(name_value_col))) %>% summarize(mean(ct1_num), .groups = 'drop')
  df[['BriTROC']] = BriTROC_CN_features[[feature]] %>% group_by(ID) %>% mutate(ct1_num = as.numeric(get(name_value_col))) %>% summarize(mean(ct1_num), .groups = 'drop')
  return(df)
}

remove_infty = function(i){
  i = i[!is.infinite(i)]
  return(i)
}

impute = function(mat, inputation_value){
  mat[mat == 0] = inputation_value
  normalise_rw(mat)
}

give_dendrogram_generalised = function(df, modify_labels=T, keep_only_PDO, plot_dendro=T){
  if(keep_only_PDO & modify_labels){stop('Only one can be true: keep_only_PDO, modify_labels')}
  
  x = hclust(dist(df))
  if(modify_labels){
    x$labels[!grepl("PDO", x$labels)] = ''
    x$labels[grepl("PDO", x$labels)] = '*'
  }
  
  if(keep_only_PDO){
    x$labels[!grepl("PDO", x$labels)] = ''
  }
  
  y = x
  y$labels = gsub("Sample ", "PDO", x$labels)
  # plot(x)
  
  labelCol <- function(x) {
    if (is.leaf(x)) {
      ## fetch label
      label <- attr(x, "label")
      ## set label color to red for A and B, to blue otherwise
      attr(x, "nodePar") <- list(lab.col="blue")
    }
    return(x)
  }
  d <- dendrapply(as.dendrogram(x), labelCol)
  if(plot_dendro)  plot(d)
  return(y)
  
}

plot_rank = function(df_rank, nudge_scalar=1, size_labels=10){
  df_rank$value = as.numeric(df_rank$value)
  # df_rank = df_rank[(df_rank$Var1),]
  df_rank$L1= factor(df_rank$L1, levels = c('BriTROC', 'organoids', 'pcawg', 'tcga'))
  df_rank$labels = ifelse(as.character(df_rank$L1) == 'organoids', as.character(df_rank$Var1), NA)
  
  df_rank = cbind.data.frame(df_rank, organoids= grepl("organoids", df_rank$L1))
  ggplot(df_rank, aes(x=factor(Var1, as.character(Var1)[order(value)]),
                      y=value, label=labels))+
    geom_bar(stat = "identity", aes( fill= organoids), width=0.5)+
    geom_bar(data = df_rank[df_rank$L1 == "organoids",], stat = "identity",
             aes(x=factor(Var1, as.character(Var1)[order(value)]),
                 y=value, label=labels, width=5), width=100)+
    lims(y=c(-max(df_rank$value)*0.3, 2+max(df_rank$value)))+
    geom_label_repel(data=df_rank[c(T,F),], aes(x=factor(Var1, as.character(Var1)[order(value)]),
                                                y=value/2, label=labels), nudge_y=400*nudge_scalar, size=size_labels)+
    geom_label_repel(data=df_rank[c(F,T),], aes(x=factor(Var1, as.character(Var1)[order(value)]),
                                                y=value/2, label=labels), nudge_y=-400*nudge_scalar, size=size_labels)+
    scale_x_discrete(expand = c(.05, 0, .05, 0))+
    scale_fill_manual(values = c("#f2a5a5", "black"))+
    theme_cowplot()+
    theme(legend.position = "bottom", axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
}

give_distance_from_imputation <- function(impute_VALUE){
  dist(as(compositions::clr(impute(all_natgen[[which_natgen]], impute_VALUE)), 'matrix'))
}

give_dendrogram_from_imputation <- function(impute_VALUE, plot=T, exposures=NULL, return_grob=F,
                                            expand_vec=c(0.5, 0, 0.05, 0), ...){

  if(!plot & return_grob){
    stop('You cannot have plot=F and return_grob=T')
  }
  
  if(is.null(exposures)){
    .exposures <- all_natgen[[which_natgen]]
  }else{
    if(is.null(rownames(exposures))){
      stop('Row names of <exposures> should not be null\n')
    }
    .exposures <- exposures
  }
  
  dendroimputclr_all_lowerinput = give_dendrogram_generalised(as(compositions::clr(impute(.exposures, impute_VALUE)), 'matrix'), modify_labels=F, keep_only_PDO = F, ...)
  
  dend_data_inputclr0004 <- dendro_data(dendroimputclr_all_lowerinput, type = "rectangle")
  dend_data_inputclr0004$labels$label = as.character(dend_data_inputclr0004$labels$label)
  dend_data_inputclr0004$labels$label[!grepl('PDO', dend_data_inputclr0004$labels$label)] = ""
  
  if(plot){
    p_v2_0004 <- ggplot(dend_data_inputclr0004$segments) +
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
      geom_label_repel(data = dend_data_inputclr0004$labels, aes(x, y, label = gsub('Organoid ', '', label)),
                       hjust = 0, size = 3, vjust=0, nudge_y = -2)+
      ylim(-3, 15)+
      theme_bw()+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank())+
      scale_x_continuous(expand = c(extra_expand, extra_expand))+
      scale_y_continuous(expand = expand_vec)
    
    ## here there used  to be a problem in that I used "labels", but now using "label" the order needs to be changed
    heatmap_dendrogram_df_inputclr0004 = t(.exposures[rownames(.exposures)[match(gsub("Organoid ", "", label(dendroimputclr_all_lowerinput)[dendroimputclr_all_lowerinput$order]),rownames(.exposures))],])
    
    p2_inputclr_0004 = ggplot(reshape2::melt(heatmap_dendrogram_df_inputclr0004), aes(x=Var2, y=value, fill=Var1))+geom_bar(stat='identity')+theme_bw()+
      theme(axis.title.x=element_blank(),  legend.title=element_blank(),
            legend.text=element_text(size=10),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),  legend.position = "bottom",
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.border = element_blank())+
      scale_fill_brewer(palette="Dark2")+
      scale_x_discrete(expand = c(extra_expand_v2, extra_expand_v2))+
      guides(fill = guide_legend(nrow = 1))
    
    # grid.arrange(p_v2_0004, p2_inputclr_0004, heights=c(2,1), top=paste0('Imputation: ', impute_VALUE))
    
      if(return_grob){
        return_plt <- (arrangeGrob(p_v2_0004, p2_inputclr_0004, heights=c(2,1),
                                                         top=paste0('Imputation: ', impute_VALUE)))
        return(return_plt)
      }
    }else{
      return(list(plot_data=dend_data_inputclr0004,
                  dendogram_data=dendroimputclr_all_lowerinput))
    }
  
}

give_pca = function(data_matrix, center = T, title='', names, names_bool=T, give_loadings=F, print_labels=T, groups=NULL,groups_shape=NULL, nrow_legend=3, size_points=1, print_both_labels=FALSE, group_is_factor=T,
                    short_pdo_names=F, additional_df=NULL){
  prcomp_res = prcomp(data_matrix, scale. = TRUE, center = center)
  eigs <- prcomp_res$sdev^2
  if(!names_bool){
    names=NA
  }
  
  give_loadings_alt <- FALSE
  # if(give_loadings){
  if(give_loadings_alt){
    if(!is.null(groups_shape)){
      stop('Groups shape not yet implemented for loadings')
    }
    a = ggplot()+
      geom_point(data=cbind.data.frame(prcomp_res$x[,1:2], names=names), aes(x=PC1, y=PC2), size=size_points)+
      geom_segment(data=cbind.data.frame(prcomp_res$rotation[,1:2], names=paste0('n', 1:ncol(data_matrix))),
                   aes(x=0, y=0, xend=PC1*4, yend=PC2*4), col='red', arrow = arrow(length = unit(0.03, "npc")))+
      labs(x=paste0('PC1 (', round(100*eigs[1]/sum(eigs), 2), '%)' ),
           y=paste0('PC2 (', round(100*eigs[2]/sum(eigs), 2), '%)' ))+ggtitle(title)
    if(print_labels){
      a = a+geom_label_repel(data=cbind.data.frame(prcomp_res$rotation[,1:2], names=paste0('n', 1:ncol(data_matrix))),
                             aes(x=PC1*3, y=PC2*3, label=names), size=3, col='red')
      if(print_both_labels){
        a = a+geom_label_repel(data=cbind.data.frame(prcomp_res$x[,1:2], names=rownames(prcomp_res$x)),
                               aes(x=PC1, y=PC2, label=names), size=3, col='black')
      }
    }
    a
  }else{
    if(!is.null(groups)){
      if(group_is_factor){
        # df = cbind.data.frame(prcomp_res$x[,1:2], names=names, groups=as.factor(groups), additional_df)
        df = cbind.data.frame(prcomp_res$x[,1:2], names=names, groups=as.factor(groups))
      }else{
        # df = cbind.data.frame(prcomp_res$x[,1:2], names=names, groups=(groups), additional_df)
        df = cbind.data.frame(prcomp_res$x[,1:2], names=names, groups=(groups))
      }
    }else{
      # df = cbind.data.frame(prcomp_res$x[,1:2], names=names, additional_df)
      df = cbind.data.frame(prcomp_res$x[,1:2], names=names)
    }
    if(short_pdo_names){
      df$names <- gsub("PDO", "", df$names)
    }
    if(!is.null(groups_shape)){
      df = cbind(df, groups_shape=groups_shape)
    }
    a = ggplot(df,
               aes(x=PC1, y=PC2,label=gsub('Sample ', 'PDO', names)))+
      labs(x=paste0('PC1 (', round(100*eigs[1]/sum(eigs), 2), '%)' ),
           y=paste0('PC2 (', round(100*eigs[2]/sum(eigs), 2), '%)' ))+ggtitle(title)
    if(!is.null(groups)){
      if(!is.null(groups_shape)){
        a = a +geom_point(aes(col=groups, shape=groups_shape), size=size_points)
      }else{
        a = a +geom_point(aes(col=groups), size=size_points)
      }
      a = a + theme(legend.position = "bottom",
                    legend.key.size = unit(0.3, "cm"),
                    legend.key.width = unit(0.2,"cm"),
                    legend.title = element_blank())+
        guides(col=guide_legend(nrow=nrow_legend,byrow=TRUE))
    }else{
      if(!is.null(groups_shape)){
        a = a + geom_point(aes(shape=groups_shape), size=size_points)
      }else{
        a = a + geom_point(size=size_points)
      }
    }
    
    
    if(print_labels){
      a = a+geom_label_repel(size=3,)
    }

    if(give_loadings){
      a <- a + geom_segment(data=cbind.data.frame(prcomp_res$rotation[,1:2], names=paste0('n', 1:ncol(data_matrix))),
                   aes(x=0, y=0, xend=PC1*4, yend=PC2*4), col='red', arrow = arrow(length = unit(0.03, "npc")))+
        geom_label_repel(data=cbind.data.frame(prcomp_res$rotation[,1:2], names=paste0('n', 1:ncol(data_matrix))),
                         aes(x=PC1*3, y=PC2*3, label=names), size=3, col='red')
    }
    
    a
    
  }
}


give_heatmap_and_dendro <- function(dendrogram_arg, exposures_arg, extra_expand = .040, extra_expand_v2 = .040){
  dend_data_inputclr <- dendro_data(dendrogram_arg, type = "rectangle")
  dend_data_inputclr$labels$label = as.character(dend_data_inputclr$labels$label)
  dend_data_inputclr$labels$label[!grepl('PDO', dend_data_inputclr$labels$label)] = ""
  
  p_v2 <- ggplot(dend_data_inputclr$segments) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
    geom_label_repel(data = dend_data_inputclr$labels, aes(x, y, label = gsub('Organoid ', '', label)),
                     hjust = 0, size = 3, vjust=0, nudge_y = -2)+
    ylim(-3, 15)+
    theme_bw()+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())+
    scale_x_continuous(expand = c(extra_expand, extra_expand))+
    scale_y_continuous(expand = c(0.05, 0, 0.05, 0))
  # print(p)
  
  heatmap_dendrogram_df_inputclr = t(exposures_arg[rownames(exposures_arg)[match(labels(dendrogram_arg),rownames(exposures_arg))],])
  
  p2_inputclr = ggplot(melt(heatmap_dendrogram_df_inputclr), aes(x=Var2, y=value, fill=Var1))+geom_bar(stat='identity')+theme_bw()+
    theme(axis.title.x=element_blank(),  legend.title=element_blank(),
          legend.text=element_text(size=10),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),  legend.position = "bottom",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_blank())+
    scale_fill_brewer(palette="Dark2")+
    scale_x_discrete(expand = c(extra_expand_v2, extra_expand_v2))+
    guides(fill = guide_legend(nrow = 1))
  
  # p2_inputclr_with_ticks = ggplot(melt(heatmap_dendrogram_df_inputclr), aes(x=Var2, y=value, fill=Var1))+geom_bar(stat='identity')+theme_bw()+
  #   theme(axis.title.x=element_blank(),  legend.title=element_blank(),
  #         legend.text=element_text(size=10),
  #         axis.text.x=element_blank(),
  #         axis.ticks.x=element_blank(),
  #         legend.position = "bottom",
  #         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #         panel.border = element_blank())+
  #   scale_fill_brewer(palette="Dark2")+
  #   scale_x_discrete(expand = c(0, extra_expand_v2))+
  #   guides(fill = guide_legend(nrow = 1))
  
  grid.arrange(p_v2, p2_inputclr, heights=c(2,1))
}

wrapper_give_heatmap <- function(arg_exposures, imputation_value = 1e-2){
  .dendro = give_dendrogram_generalised(as(compositions::clr(impute(arg_exposures, imputation_value)), 'matrix'),
                                        modify_labels=F, keep_only_PDO = F)
  give_heatmap_and_dendro(dendrogram_arg = .dendro, exposures_arg = arg_exposures,
                          extra_expand = 0, extra_expand_v2 = 0)
}

give_amalgamation <- function(i, list_amalgamations){
  new_mat = sapply(list_amalgamations, function(j){
    grouped_exp <- i[,colnames(i) %in% j]
    if(!is.null(ncol(grouped_exp))){
      rowSums(grouped_exp)
    }else{
      grouped_exp
    }
  })
  new_mat
}


give_annotation_from_names <- function(i){
  sapply(i, function(j){
    if(grepl('^TCGA', j)){
      'TCGA'
    }else if(grepl('^IM', j)){
      'BriTROC-1'
    }else if(grepl('^JBLAB', j)){
      'BriTROC-1'
    }else if(grepl('^PDO', j)){
      'organoid'
    }else{
      'ICGC'
    }
  })
}
