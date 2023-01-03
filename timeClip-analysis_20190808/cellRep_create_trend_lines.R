#' Function to create PC1 lines from a list with PCs
#'   and their most important gene
#'   
#' @param pcs_and_their_genes a list of meta_classes that contains a list
#'   with pcs and list/matrix of genes
#' @param relative take the relative (to the max) of the PC value
#' @param title the title of the plot
#' @param palette RcolorBrewer palette name
#' @param discrete_palette is present overwrite palette
#'   
#' @return ggplot
#'
create_trend_lines_for_meta_classes <- function(pcs_and_their_genes, relative=TRUE,
                                                title="Trend", palette="Set1",
                                                discrete_palette=NULL) {
  library(ggplot2)
  multi_df <- list()
  for (meta_class in names(pcs_and_their_genes)){
    pcs <- pcs_and_their_genes[[meta_class]]$pcs
    day <- as.numeric(sub("d", "",colnames(pcs)))
    PC1 = pcs[1,]
    ylab_value = "PC1"
    if (relative) {
      PC1=PC1/max(abs(PC1))
      ylab_value = paste("Relative", ylab_value)
    }
    
    multi_df[[meta_class]] <- data.frame(PC1=PC1,day=day, meta_class=meta_class)
  }
  data <- do.call(rbind, multi_df)
  p <- ggplot(data) +
    geom_line(aes(x=day, y=PC1, color=meta_class), lwd=2) +
    theme_bw() +
    ylab(ylab_value) +
    scale_x_continuous(breaks = data$day, labels = data$day) +
    ggtitle(title) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  if (!is.null(discrete_palette)){
    p = p + scale_colour_manual(values=discrete_palette)
  } else {
    p = p + scale_color_brewer(palette=palette)
  }
  
  p
}

create_trend_lines_for_genes <- function(pcs_and_their_genes, expr_data, relative=TRUE,
                                         scale=T,
                                         title="Trend", palette="Set1",
                                         discrete_palette=NULL) {
  library(ggplot2)
  source("extract_PC_average_profiles.R")
  
  multi_df <- list()
  for (meta_class in names(pcs_and_their_genes)){
    pc1_genes <- pcs_and_their_genes[[meta_class]]$pcs2genes[,1]
    pc1_avg <- extract_average_gene_profiles(pc1_genes, expr_data,scale = scale)
    day <- as.numeric(sub("d", "",names(pc1_avg)))
    ylab_value = "PC1 avg"
    if (relative) {
      pc1_avg=pc1_avg/max(abs(pc1_avg))
      ylab_value = paste("Relative", ylab_value)
    }
    
    multi_df[[meta_class]] <- data.frame(PC1=pc1_avg,day=day, meta_class=meta_class)
  }
  data <- do.call(rbind, multi_df)
  p <- ggplot(data) +
    geom_line(aes(x=day, y=PC1, color=meta_class), lwd=2) +
    theme_bw() +
    ylab(ylab_value) +
    scale_x_continuous(breaks = data$day, labels = data$day) +
    ggtitle(title) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  if (!is.null(discrete_palette)){
    p = p + scale_colour_manual(values=discrete_palette)
  } else {
    p = p + scale_color_brewer(palette=palette)
  }
  
  p
}


create_trend_lines_comparison <- function(pcs_and_their_genes, meta_class,
                                          expr_data, expr_data_cmp, relative=TRUE, 
                                          scale=TRUE,
                                          title="Trend", palette="Set1",
                                          discrete_palette=NULL) {
  library(ggplot2)

  pc1_genes <- pcs_and_their_genes[[meta_class]]$pcs2genes[,1]
  pc1_avg <- extract_average_gene_profiles(pc1_genes, expr_data, scale=scale)
  
  if (!length(pc1_avg))
    stop("At least one gene must be present in the first comparison term.")
  
  day <- as.numeric(sub("d", "",names(pc1_avg)))
  ylab_value = "PC1 avg"
  
  if (relative) {
    pc1_avg=pc1_avg/max(abs(pc1_avg))
    ylab_value = paste("Relative", ylab_value)
  }
  
  data <- data.frame(PC1=pc1_avg,day=day, condition="wt")
  
  pc1_avg_cmp <- extract_average_gene_profiles(pc1_genes, expr_data_cmp, scale = scale)
  if (!is.null(pc1_avg_cmp)){
    day_cmp <- as.numeric(sub("d", "",names(pc1_avg_cmp)))
    if (relative)
      pc1_avg_cmp=pc1_avg_cmp/max(abs(pc1_avg_cmp))
    data <- rbind(data, data.frame(PC1=pc1_avg_cmp,day=day_cmp, condition="ko"))
  }
  
  p <- ggplot(data) +
    geom_line(aes(x=day, y=PC1, color=condition), lwd=2) +
    theme_bw() +
    ylab(ylab_value) +
    scale_x_continuous(breaks = data$day, labels = data$day) +
    ggtitle(title) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  if (!is.null(discrete_palette)){
    p = p + scale_colour_manual(values=discrete_palette)
  } else {
    p = p + scale_color_brewer(palette=palette)
  }
  
  p
}
