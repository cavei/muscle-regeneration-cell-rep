# Day by Day Cell Population - GO BP Analysis with enricher

library(org.Mm.eg.db)
library(clusterProfiler)
library(ggplot2)
library(cellCB)

# We load the list of terms we manually selected for the importance in muscle regeneration.

manually_curated_GO_list <- read.table("../manual-selection-of-go-terms/keep_me.txt", sep="\t",
                      quote="\"", stringsAsFactors = F, check.names = F)
manually_curated_GO_list <- manually_curated_GO_list$V1

# We created the database of all the BP terms. The database is created using a wrapper of clusterProfiler function. This speed up the analysis. If the file has alreasy been computer, the database is simply loaded.

if (!(file.exists("go_bp_full_data.RData"))){
  go_bp_full_data <- cellCB:::create.GO_DB()
  save(go_bp_full_data, file="go_bp_full_data.RData")
} else {
  load("go_bp_full_data.RData")
}

# We set up the remote dir that contains the results and one where we are going to save our tables.

resultsDir <- "../preprocessing-data/rpkm-RDatas-padj0.01-lfc2/"

tableDir = "GO_enrichment_tables"
dir.create(tableDir, showWarnings = F)

# We start from the total muscle wt. We load the genes that are active at each day and the normalized expression matrix.

activeAtDaysTot <- cellCB:::getActiveAtDaysFromRData(paste0(resultsDir, "/tot-wt.RData"), "tot_wt")

expressionNormTot <- cellCB::getNormExpressionFromRData(paste0(resultsDir, "/", "tot", "-wt.RData"), paste0("tot", "_wt"))
cellExpressedTot <- row.names(expressionNormTot)

# We are going to perform the enrichment. We perform a regular enrichment using compareCluster (clusterProfiler R package). The results are returned in a data frame. This data frame is filtered to contains only the GO category that we manually curated.

if (file.exists("clusterCompare_total_muscle.RData")) {
  load("clusterCompare_total_muscle.RData")
} else {
  set.seed(1234)
  clusterCompare_total_muscle.df <- my.compareCluster(activeAtDaysTot, go_bp_full_data,
                                      minGSSize=10,
                                      universe = cellExpressedTot, 
                                      pvalueCutoff = 0.1, qvalueCutoff = 0.05)
  manual_selection <- clusterCompare_total_muscle.df$df$Description %in% manually_curated_GO_list
  clusterCompare_total_muscle <- clusterCompare_total_muscle.df$df[manual_selection, , drop=F]
  
  clusterCompare_total_muscle <- new("compareClusterResult",
                                     compareClusterResult = clusterCompare_total_muscle,
                                     geneClusters = clusterCompare_total_muscle.df$df,
                                     fun = clusterCompare_total_muscle.df$fun,
                                     .call = clusterCompare_total_muscle.df$.call)
  
  save(clusterCompare_total_muscle, file="clusterCompare_total_muscle.RData")

# Using the dotplot function we can plot the enriched categories. Using the parameter _showCategory = 10_ we are telling the function that we want to show, for each day, the top 10 enriched GO BP category. To cope with the overwelming quantity of data, we selected the GO BP that fit the most with the case study.
# We created a list of the category that were keeped for each cell population.

keeped_go <- list()

# We are going to collapse similar GO categoy into meta classes (meta category) to simplify the visualization. We manually created the GO category clusters (meta classes). We store them and load the for the analysis.

library(clusterProfiler)

create_brand_new_filter_annotation <- FALSE

cell="tot"
sigTable <- as.data.frame(clusterCompare_total_muscle)

if (create_brand_new_filter_annotation) {
  terms <- sort(unique(sigTable$Description))
  # terms <- cellCB:::exclude_go_terms_by_counts(sigTable, lowerCut = 2)
  pure_terms <- cellCB:::remove_words(terms)
  store = "manual-filtered-and-annotated-terms"
  dir.create(store, showWarnings = F)
  output <- cbind(names(pure_terms), pure_terms)
  colnames(output) <- c("full", "pure_terms")
  write.table(output, file=paste0(store, "/my_terms_filtered_190726.txt"), row.names=F, sep="\t", quote=F)
}

store = "manual-filtered-and-annotated-terms"

if (!file.exists(paste0(store,"/total_GO_clusters.RData"))) {
  manual_cluster <- read.table(paste0(store,"/190726-terms_filtered_new_classes_EG.txt"), header=T, sep="\t", quote="\"", stringsAsFactors = F, check.names = F)
  
  if (length(setdiff(sigTable$Description, manual_cluster$full)) != 0)
    warning("Consider re reunning the manual selection")
  
  total_GO_clusters <- manual_cluster$EG_NEW
  names(total_GO_clusters) <- manual_cluster$full
  
  total_GO_Cluster_map <- do.call(rbind, lapply(seq_along(total_GO_clusters), function(i) {
    orid <- names(total_GO_clusters)[i]
    meta_cat <- unlist(strsplit(total_GO_clusters[i], "/"))
    data.frame(orid, meta_cat, stringsAsFactors = F)
  }))
  row.names(total_GO_Cluster_map) <- NULL
  
  save(total_GO_clusters, total_GO_Cluster_map, file=paste0(store, "/","total_GO_clusters.RData"))
} else {
  load(paste0(store,"/total_GO_clusters.RData"))
}

# "total_GO_clusters" is a named vector: names are the category, values are the meta classes. Using our function guess_best_sets we are going to identify the most representative GO category that will represent the meta class for the total muscle. We save the result in the keeped_go list.

keep_go <- unname(cellCB:::guess_best_sets_map(total_GO_Cluster_map, sigTable))
keeped_go[[cell]] <- list(keep_go=keep_go, grp=cellCB:::groupByCluster_map(total_GO_Cluster_map))

# We can create a plot clustering the terms in manually curated superclasses for the total muscle.

cellCB:::personal.dotplot.compareClusterResult(clusterCompare_total_muscle, showCategory=10,
                                               title=paste0(cell, "-wt"),
                                               filterSets = keep_go,
                                               relabel_description = keeped_go[[cell]]$grp) + 
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10))

# We saved the table that created the plot in "tableDir" a.k.a. "GO_enrichment_tables".

write.table(sigTable, file=paste0(tableDir, "/full-", cell, "-wt-enrichGO.txt"), sep="\t", quote=F)

# The GO terms that were found in the total muscle are going to be reused in the cell population analysis.
# We saved the list in a new object called valid_go_clusters.

valid_go_clusters <- total_GO_Cluster_map$orid

# Now it is time to analyze the cell populations. We are going to load expression and the lists of the genes active at each day.

cells <- c("fap", "ec", "mp", "inf", "per")

expressionNormCell <- lapply(cells, function(cell) {
  getNormExpressionFromRData(paste0(resultsDir, "/", cell, "-wt.RData"), paste0(cell, "_wt"))
})
names(expressionNormCell) <- cells

cellExpressed <- lapply(expressionNormCell, row.names)

activeAtDaysCell <- lapply(cells, function(cell) {
  getActiveAtDaysFromRData(paste0(resultsDir, "/", cell, "-wt.RData"), paste0(cell, "_wt"))
})
names(activeAtDaysCell) <- cells

load("../associate-gene-to-cell/gene2cellAssociation.RData")

genesBehavioursCell <- lapply(cells, function(cell) {
  extractGenesBehavious(activeAtDaysCell[[cell]], activeAtDaysTot, gene2cellAssociation[[cell]])
})
names(genesBehavioursCell) <- cells

movingGenesCell <- lapply(genesBehavioursCell, function(cell){
  moving <- lapply(names(cell$activeAtDays), findAWay::defineDailyExpressed, obj=cell)
  names(moving) <- names(cell$activeAtDays)
  moving
})


# As you noticed we loaded also the list of genes constitutively active in the different cells.
# With this lists we can perform the enrichment analysis using clusterProfiler wrapper.

if (file.exists("cell_wise_GO_enrichment_dfs.RData")) {
  load("cell_wise_GO_enrichment_dfs.RData")
} else {
  cell_wise_GO_enrichment_dfs <- lapply(names(activeAtDaysCell), function(cell) {
    my.compareCluster(activeAtDaysCell[[cell]], go_bp_full_data,
                                      minGSSize=10,
                                      universe = cellExpressed[[cell]], 
                                      pvalueCutoff = 0.1, qvalueCutoff = 0.05)
  })
  names(cell_wise_GO_enrichment_dfs) <- names(activeAtDaysCell)
  save(cell_wise_GO_enrichment_dfs, file="cell_wise_GO_enrichment_dfs.RData")
}

# As previously seen for the total muscle, we created the enrichment data frames. Now we need to filter the results by the GO used in the total tissue that have been stored in "valid_go_clusters".

cell_wise_GO_enrichment <- lapply(names(cell_wise_GO_enrichment_dfs), function(cell) {
  df <- cell_wise_GO_enrichment_dfs[[cell]]$df
  geneCls <- activeAtDaysCell[[cell]]
  small.df <- df[as.character(df$Description) %in% valid_go_clusters, ,drop=F]
  new("compareClusterResult", compareClusterResult = small.df, geneClusters = geneCls, fun = "my.enrich", .call = match.call())
})
names(cell_wise_GO_enrichment) <- names(cell_wise_GO_enrichment_dfs)

# After collecting the enriched object already filtered

cell="fap"
cell_sigTable <- as.data.frame(cell_wise_GO_enrichment[[cell]])

keeped_go[[cell]] = sigTable_filter_by_clusters_map(cell_sigTable, total_GO_Cluster_map)

cellCB:::personal.dotplot.compareClusterResult(cell_wise_GO_enrichment[[cell]], showCategory=10,
                                               title=paste0(cell, "-wt"),
                                               filterSets = keeped_go[[cell]]$keep_go,
                                               relabel_description = keeped_go[[cell]]$grp) +
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10))

cell="ec"
cell_sigTable <- as.data.frame(cell_wise_GO_enrichment[[cell]])
keeped_go[[cell]] = sigTable_filter_by_clusters_map(cell_sigTable, total_GO_Cluster_map)

cellCB:::personal.dotplot.compareClusterResult(cell_wise_GO_enrichment[[cell]], showCategory=10,
                                               title=paste0(cell, "-wt"),
                                               filterSets = keeped_go[[cell]]$keep_go,
                                               relabel_description = keeped_go[[cell]]$grp) +
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10))

cell="mp"
cell_sigTable <- as.data.frame(cell_wise_GO_enrichment[[cell]])
keeped_go[[cell]] = sigTable_filter_by_clusters_map(cell_sigTable, total_GO_Cluster_map)

cellCB:::personal.dotplot.compareClusterResult(cell_wise_GO_enrichment[[cell]], showCategory=10, 
                                               title=paste0(cell, "-wt"),
                                               filterSets = keeped_go[[cell]]$keep_go,
                                               relabel_description = keeped_go[[cell]]$grp) +
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10))

cell="inf"
cell_sigTable <- as.data.frame(cell_wise_GO_enrichment[[cell]])

keeped_go[[cell]] = sigTable_filter_by_clusters_map(cell_sigTable, total_GO_Cluster_map)

cellCB:::personal.dotplot.compareClusterResult(cell_wise_GO_enrichment[[cell]], showCategory=10, 
                                               title=paste0(cell, "-wt"),
                                               filterSets = keeped_go[[cell]]$keep_go,
                                               relabel_description = keeped_go[[cell]]$grp) +
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10))

cell="per"
cell_sigTable <- as.data.frame(cell_wise_GO_enrichment[[cell]])
keeped_go[[cell]] = sigTable_filter_by_clusters_map(cell_sigTable, total_GO_Cluster_map)

cellCB:::personal.dotplot.compareClusterResult(cell_wise_GO_enrichment[[cell]], showCategory=10, 
                                               title=paste0(cell, "-wt"),
                                               filterSets = keeped_go[[cell]]$keep_go,
                                               relabel_description = keeped_go[[cell]]$grp) +
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10))

# We are going to save the tables of the enrichment. You will find them in the files in dir "GO_enrichment_tables".

for (cell in names(cell_wise_GO_enrichment)){
  table <- as.data.frame(cell_wise_GO_enrichment[[cell]])
  write.table(table, file=paste0(tableDir, "/full-", cell, "-wt-enrichGO.txt"), sep="\t", quote=F)
}

# Now we analyze the KO cell population.
# We need to load the data for the KOs.

cells <- c("fap", "ec", "mp")
expressionNormCell.ko <- lapply(cells, function(cell) {
  getNormExpressionFromRData(paste0(resultsDir, "/", cell, "-ko.RData"), paste0(cell, "_ko"))
})
names(expressionNormCell.ko) <- cells
cellExpressed.ko <- lapply(expressionNormCell.ko, row.names)

activeAtDaysCell.ko <- lapply(cells, function(cell) {
  getActiveAtDaysFromRData(paste0(resultsDir, "/", cell, "-ko.RData"), paste0(cell, "_ko"))
})
names(activeAtDaysCell.ko) <- cells

# In this case we don't have the total muscle ko dataset.
# We can enrich. The same setting used for total muscle and WT were also used for KO. We started by creating the series of enriched data frames.

if (file.exists("cell_wise_GO_enrichment_ko_dfs.RData")) {
  load("cell_wise_GO_enrichment_ko_dfs.RData")
} else {
  cell_wise_GO_enrichment_ko_dfs <- lapply(names(activeAtDaysCell.ko), function(cell) {
    my.compareCluster(activeAtDaysCell.ko[[cell]], go_bp_full_data,
                                      minGSSize=10,
                                      universe = cellExpressed.ko[[cell]], 
                                      pvalueCutoff = 0.1, qvalueCutoff = 0.05)
  })
  names(cell_wise_GO_enrichment_ko_dfs) <- names(activeAtDaysCell.ko)
  save(cell_wise_GO_enrichment_ko_dfs, file="cell_wise_GO_enrichment_ko_dfs.RData")
}

# Next, we created the object clusterCompare filtering by the classes stored in the valid clusters.

cell_wise_GO_enrichment_ko <- lapply(names(cell_wise_GO_enrichment_ko_dfs), function(cell) {
  df <- cell_wise_GO_enrichment_ko_dfs[[cell]]$df
  geneCls <- activeAtDaysCell.ko[[cell]]
  small.df <- df[df$Description %in% valid_go_clusters, ,drop=F]
  new("compareClusterResult", compareClusterResult = small.df, geneClusters = geneCls, fun = "my.enrich", .call = match.call())
})
names(cell_wise_GO_enrichment_ko) <- names(cell_wise_GO_enrichment_ko_dfs)

## Create the data frame.
# We explore the GO BP for FAP KO dataset performing the same plots. Once again with _showCategory=10_ parameter. In this case we used the category keeped for FAP wt.

cellCB:::personal.dotplot.compareClusterResult(cell_wise_GO_enrichment_ko[[cell]], showCategory=10,
                                               title=paste0(cell, "-ko"), 
                                               filterSets = keeped_go[[cell]]$keep_go,
                                               relabel_description = keeped_go[[cell]]$grp) + 
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10))

# EC KO

cell="ec"
cellCB:::personal.dotplot.compareClusterResult(cell_wise_GO_enrichment_ko[[cell]], showCategory=10,
                                               title=paste0(cell, "-ko"), 
                                               filterSets = keeped_go[[cell]]$keep_go,
                                               relabel_description = keeped_go[[cell]]$grp) + 
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10))


# MP KO

cell="mp"
cellCB:::personal.dotplot.compareClusterResult(cell_wise_GO_enrichment_ko[[cell]], showCategory=10,
                                               title=paste0(cell, "-ko"),
                                               filterSets = keeped_go[[cell]]$keep_go,
                                               relabel_description = keeped_go[[cell]]$grp) + 
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10))

# We are going to save the tables of the enrichment also for the KO. You will find them in the GO_enrichment_tables dir.

for (cell in names(cell_wise_GO_enrichment_ko)){
  table <- as.data.frame(cell_wise_GO_enrichment_ko[[cell]])
  write.table(table, file=paste0(tableDir, "/", cell, "-ko-enrichGO.txt"), sep="\t", quote=F)
}

## Compare the KO and WT

remove_uneven_days <- TRUE
cellDEG_df <- lapply(names(activeAtDaysCell.ko), function(cell) {
  cellDEG_df.ko <- do.call(rbind, 
                           lapply(names(activeAtDaysCell.ko[[cell]]), function(d){
                             data.frame(symbol=activeAtDaysCell.ko[[cell]][[d]], day=d,
                                        cell=cell, condition="ko", stringsAsFactors = F)
                             }))
  
  cellDEG_df.wt <- do.call(rbind,
                           lapply(names(activeAtDaysCell[[cell]]), function(d){
                             data.frame(symbol=activeAtDaysCell[[cell]][[d]], day=d,
                                        cell=cell, condition="wt", stringsAsFactors = F)
                             }))
  
  cellDEG_df <- rbind(cellDEG_df.wt,cellDEG_df.ko)
  
  if (remove_uneven_days) {
    remove_this_days = names(which(table(c(names(activeAtDaysCell.ko[[cell]]), names(activeAtDaysCell[[cell]])))==1))
    if (!is.null(remove_this_days))
      cellDEG_df <- cellDEG_df[!cellDEG_df$day %in% remove_this_days, , drop=F]
  }
  
  day_label <- as.numeric(gsub("d","",cellDEG_df$day))
  cellDEG_df$day[day_label<=9] <- paste0("d0",day_label[day_label<=9])
  cellDEG_df
})
names(cellDEG_df) <- names(activeAtDaysCell.ko)

if (file.exists("cell_wise_GO_enrichment_wt_vs_ko_dfs.RData")) {
  load("cell_wise_GO_enrichment_wt_vs_ko_dfs.RData")
} else {
  cell_wise_GO_enrichment_wt_vs_ko_dfs <- lapply(names(cellDEG_df), function(cell){
    my.compareCluster(symbol~day+condition, go_bp_full_data, data=cellDEG_df[[cell]],
                      minGSSize=10,
                      universe = unique(cellExpressed[[cell]], cellExpressed.ko[[cell]]),
                      pvalueCutoff = 0.1, qvalueCutoff = 0.05)
  })
  names(cell_wise_GO_enrichment_wt_vs_ko_dfs) <- names(cellDEG_df)
  save(cell_wise_GO_enrichment_wt_vs_ko_dfs, file="cell_wise_GO_enrichment_wt_vs_ko_dfs.RData")
}

cell_wise_GO_enrichment_wt_vs_ko <- lapply(cell_wise_GO_enrichment_wt_vs_ko_dfs, function(clsRes) {
  df <- clsRes$df
  geneCls <- clsRes$geneClusters
  small.df <- df[df$Description %in% valid_go_clusters, ,drop=F]
  new("compareClusterResult", compareClusterResult = small.df, geneClusters = geneCls,
      fun = clsRes$fun, .call = clsRes$.call)
})

cell="ec"
cellCB:::personal.dotplot.compareClusterResult(cell_wise_GO_enrichment_wt_vs_ko[[cell]], x=~day, showCategory=20, title=paste0(cell), filterSets = keeped_go[[cell]]$keep_go, relabel_description = keeped_go[[cell]]$grp) +
  ggplot2::facet_grid(~condition) +
  theme(axis.text.x = element_text(size=10, angle=80, hjust=1), axis.text.y = element_text(size=8))

cell="fap"
cellCB:::personal.dotplot.compareClusterResult(cell_wise_GO_enrichment_wt_vs_ko[[cell]], x=~day, showCategory=10, title=paste0(cell), filterSets = keeped_go[[cell]]$keep_go, relabel_description = keeped_go[[cell]]$grp) +
  ggplot2::facet_grid(~condition) +
  theme(axis.text.x = element_text(size=10, angle=80, hjust=1), axis.text.y = element_text(size=8))

cell="mp"
cellCB:::personal.dotplot.compareClusterResult(cell_wise_GO_enrichment_wt_vs_ko[[cell]], x=~day, showCategory=10, title=paste0(cell), filterSets = keeped_go[[cell]]$keep_go, relabel_description = keeped_go[[cell]]$grp) +
  ggplot2::facet_grid(~condition) +
  theme(axis.text.x = element_text(size=10, angle=80, hjust=1), axis.text.y = element_text(size=8))

# ec-wt-vs-ko-enrichGO
for (cell in names(cell_wise_GO_enrichment_wt_vs_ko)){
  table <- as.data.frame(cell_wise_GO_enrichment_wt_vs_ko[[cell]])
  write.table(table, file=paste0(tableDir, "/", cell, "-wt-vs-ko-auto-paracrine-enrich-go.txt"), sep="\t", quote=F)
}

