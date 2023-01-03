extract_pc_profile_and_related_genes <- function(cell_data) {
  
  lapply(cell_data, function(meta_class) {
  pcs_values <- t(meta_class$pcs)
  lds <- meta_class$cliqueLoadings
  # revert sign of PC accotding to its best load gene
  best_idx_per_pc <- apply(abs(lds),2, which.max) # choose the best gene
  best_per_pc <- sapply(seq_along(best_idx_per_pc), function(i) {
    lds[best_idx_per_pc[i],i]
  }) #extract the loading value for that gene
  names(best_per_pc) <- row.names(lds)[best_idx_per_pc]
  pcs <- pcs_values * round(best_per_pc) # invert sign when needed
  
  important_genes <- apply(lds != 0, 2, function(x) {
    idx <- which(x)
    row.names(lds)[idx]
  }) #extract the genes with load != 0
  return(list(pcs=pcs, pcs2genes=important_genes))
  })
}
