extract_genes_of_the_8GO_map <- function(table, cluster_map) {
  
  # table=sigTable; cluster_map=reduced_8GO_Clusters_map
  go_categories <- as.character(table$Description)
  selector <- go_categories %in% cluster_map[,1]
  
  good_categories <- table$Description[selector]
  good_geneID <- table$geneID[selector]
  
  categories <- cluster_map[cluster_map[,1] %in% good_categories, , drop=F]
  
  tapply(categories[,1], categories[,2], function(gos) {
    genes <- good_geneID[go_categories%in% gos]
    unique(unlist(strsplit(genes, "/")))
  })
}

invert_map <- function(map) {
  if (is.null(names(map)))
    stop("No names provided. Impossible to invert the map")
  
  if (any(duplicated(map)))
    stop("Duplicated future names in vector map. Impossible to invert")
  
  imap <- names(map)
  names(imap) <- map
  imap
}
