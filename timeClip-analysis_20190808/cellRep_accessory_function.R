
check_significance_of_pathways <- function(timeClipRes, pathways, cell_cond="ec-wt", pathDB="kegg", watchCliques=TRUE) {
  timeRes <- timeClipRes[[cell_cond]][[pathDB]][pathways]
  
  if (watchCliques){
    sig <-sapply(timeRes, function(k) {
      any(k$alphas <= 0.05)
    })
  } else {
    sig <-sapply(timeRes, function(k) {
      any(k$alphaPCs <= 0.05)
    })
  }
  names(which(sig))
}

extrac_pathway_sig_cliques_pca_master_profile <- function(tc_results, sig_list, cell_cond) {
  # cell_cond <- "ec-wt"
  idx <- which(names(tc_results[[cell_cond]]$reactome) %in% sig_list[[cell_cond]])
  path_idx = idx[1]
  
  multi_mat <- lapply(idx, function(path_idx) {
    path <- tc_results[[cell_cond]]$reactome[[path_idx]]
    sig_idx <- which(path$alphas <= 0.05)
    sigPCS <- lapply(seq_along(sig_idx), function(id) {
      path$cliquesPCs[[id]][, path$bestPcs[[id]]]
    })
    
    matrix_pcs <- do.call(rbind, sigPCS)
    row.names(matrix_pcs) <- paste0("clique_", sig_idx)
    # pheatmap::pheatmap(matrix_pcs, cluster_cols = F, cluster_rows = F, scale = "row")
    matrix_pcs
  })
  names(multi_mat) <- names(tc_results[[cell_cond]]$reactome)[idx]
  multi_mat

}


extract_shared_time_point <- function(x, y) {
  xd <- sapply(strsplit(colnames(x), "_"), function(e) e[1])
  yd <- sapply(strsplit(colnames(y), "_"), function(e) e[1])
  
  list(x_names=xd, y_names=yd, shared=intersect(xd, yd))
}

extrac_clique_best_pca_master_profile <- function(tc_results, pathway, clique, cell_cond, best=T) {
  # cell_cond <- "ec-wt"
  path <- tc_results[[cell_cond]]$reactome[[pathway]]
  
  if (!best)
    return(t(path$cliquesPCs[[clique]]))
  best <- unlist(path$bestPcs[clique])
  t(path$cliquesPCs[[clique]][, best, drop=F])
}

extrac_pathway_pca <- function(tc_results, pathway, cell_cond) {
  # cell_cond <- "ec-wt"
  tc_results[[cell_cond]]$reactome[[pathway]]$cliquePCs
}

extrac_pathway_loadings <- function(tc_results, pathway, cell_cond) {
  # cell_cond <- "ec-wt"
  tc_results[[cell_cond]]$reactome[[pathway]]$cliqueLoadings
}
