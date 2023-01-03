loadCellRData <- function(RData, objname) {
  load(RData)
  get(objname)
}

annotation_matrix <- function(names) {
  classes <- sapply(strsplit(names, "_"), function(x) paste0("d",x[2]))
  m <- sapply(unique(classes), function(class) {
    classes==class
  })
  row.names(m) <- names
  m
}

create_meta_days <- function(days, breaks = c(0,1,4,10,15)-1, labels=c("steady", "early", "middle", "late","inf")) {
  names(breaks) <- labels
  rep <- cut(days, breaks=breaks)
  levels(rep) <- names(breaks[-length(breaks)])
  variables <- as.character(rep)
  var_levels <- levels(rep)
  factor(variables, levels=var_levels)
}

get_pcs <- function(data_obj) {
  # require(EDASeq)
  eset <- data_obj$esetRUV$setEmpirical
  norm <- EDASeq::normCounts(eset)
  # degs <- data_obj$DEGS
  scaleNorm <- scale(t(log2(norm+1)), scale=FALSE, center = TRUE)
  s <- svd(scaleNorm)
  pcs <- data.frame(s$u)
  colnames(pcs) <- paste0("PC", seq_along(row.names(scaleNorm)))
  row.names(pcs) <- row.names(scaleNorm)
  pcs
}

additional_clustering_index <- function(pcs, clusters) {
  days <- sapply(strsplit(row.names(pcs), "_"), function(x) as.numeric(x[2]))
  annotation <- tapply(row.names(pcs), paste0("d",days), as.character)
  ann_mat <- annotation_matrix(row.names(pcs))
  meta_days <- create_meta_days(days)
  
  intern <- clValid(pcs[,1:2], nClust = 2:8, 
                    clMethods = c("hierarchical","kmeans","pam"),
                    validation = c("internal"))
  
  biological_stability_index <- BSI(clusters, meta_days, annotation = ann_mat)
  biological_homogeneity_index = BHI(clusters, annotation = ann_mat)
  list(cl_valid = intern, biological_homogeneity_index=biological_homogeneity_index, biological_stability_index=biological_stability_index)
}

format_axis <- function() {
  return(list(theme(axis.text=element_text(size=12))))
}

evaluate_BHI <- function(pcs, clusters, ylim=NULL) {
  require(ggplot2)
  days <- sapply(strsplit(row.names(pcs), "_"), function(x) as.numeric(x[2]))
  annotation <- tapply(row.names(pcs), paste0("d",days), as.character)
  ann_mat <- annotation_matrix(row.names(pcs))
  meta_days <- create_meta_days(days)
  
  bhis <- sapply(clusters, function(k) {
    set.seed(1)
    km_clusters <- kmeans(pcs[,1:2], k, 100)
    BHI(factor(km_clusters$cluster), ann_mat)
  })
  df <- data.frame(clusters=factor(as.character(clusters), levels = as.character(clusters)),
                   BHI=bhis)

  plot <- ggplot(df, aes(x=clusters, y=BHI, group=1)) + 
    geom_point() + geom_line() +
    # xlab("bio homogeneity index") + 
    theme(legend.position = "none") +
    # theme(axis.title.y = element_blank()) +
    ylab("bio homogeneity") 
  
  if (!is.null(ylim))
    plot <- plot + ylim(ylim)
  list(plot=plot, df=df, best=clusters[which.max(bhis)])
}

silhuette_plot <- function(nbclust_obj, ylim=NULL) {
  df_silh <- data.frame(clusters=factor(names(nbclust_obj$All.index), levels=names(nbclust_obj$All.index)), 
                        silh=nbclust_obj$All.index)
  silh_plot <- ggplot(df_silh, aes(x=clusters, y=silh, group=1)) +
    geom_point() + geom_line() + 
    # xlab("silhuette") +
    theme(legend.position = "none") + 
    # theme(axis.title.y = element_blank()) +
    ylab("silhuette")
  
  if (!is.null(ylim))
    silh_plot <- silh_plot + ylim(ylim)
  return(list(plot=silh_plot, df=df_silh))
}

plot_2index_clusters <- function(pcs_data, silh_obj, bhi_obj, title="", choose=NULL) {
  library(ggplot2)
  library(cowplot)
  library(ggrepel)
  if (!all(c("cluster_silh", "cluster_bhi") %in% colnames(pcs_data) ))
    stop("cluster information is missing in pcs_data")
  k_silh = length(levels(pcs_data$cluster_silh))
  k_bhi = length(levels(pcs_data$cluster_bhi))
  
  days <- sapply(strsplit(row.names(pcs_data), "_"), function(x) as.numeric(x[2]))
  
  title <- ggdraw() + 
    draw_label(title, fontface = 'bold', size=8, x = 0, hjust = 0) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7)
    )
  
  r1 <- cowplot::plot_grid(silh_obj$plot, bhi_obj$plot, ncol=2, nrow = 1)
  r2 <- ggplot(pcs_data, aes(x=PC1,y=PC2, color=cluster_silh)) + geom_point() +
    geom_label_repel(label=days) +
    theme_classic() +
    theme(legend.position = "none") +
    # ggtitle(paste0("K = ", k_silh)) +
    format_axis()
  
  r3 <- ggplot(pcs_data, aes(x=PC1,y=PC2, color=cluster_bhi)) + geom_point() + 
    geom_label_repel(label=days) +
    theme_classic() + 
    theme(legend.position = "none") +
    # ggtitle(paste0("K = ", k_bhi)) +
    format_axis()
  
  if (is.null(title)) {
    if (is.null(choose)){
      cowplot::plot_grid(r1, r2, r3, nrow=3, ncol = 1, rel_heights = c(0.17,0.415,0.415))
    } else if (choose=="silhuette") {
      cowplot::plot_grid(r1, r2, nrow=2, ncol = 1, rel_heights = c(0.34,0.66))
    } else if (choose=="bhi") {
      cowplot::plot_grid(r1, r3, nrow=2, ncol = 1, rel_heights = c(0.34,0.66))
    } else {
      stop("Unrecognized method choose between silhuette or bhi")
    }
  } else {
    if (is.null(choose)){
      cowplot::plot_grid(title, r1, r2, r3, nrow=4, ncol = 1, rel_heights = c(0.03, 0.17,0.4,0.4))
    } else if (choose=="silhuette") {
      cowplot::plot_grid(title, r1, r2, nrow=3, ncol = 1, rel_heights = c(0.06, 0.34,0.6))
    } else if (choose=="bhi") {
      cowplot::plot_grid(title, r1, r3, nrow=3, ncol = 1, rel_heights = c(0.06, 0.34,0.6))
    } else {
      stop("Unrecognized method choose between silhuette or bhi")
    }
  }
  
}
