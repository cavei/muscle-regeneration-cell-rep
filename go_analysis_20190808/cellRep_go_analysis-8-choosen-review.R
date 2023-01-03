
# Day by Day Cell Population - GO BP Analysis with enricher


library(org.Mm.eg.db)
library(clusterProfiler)
library(ggplot2)
library(cellCB)
source("cellRep_extract_genes_of_the_8GO_map.R")

# We set up the remote dir that contains the results.

resultsDir <- "../preprocessing-data/rpkm-RDatas-padj0.01-lfc2/"
tableDir = "GO_enrichment_tables_8_class_review"
dir.create(tableDir, showWarnings = F)


# Let's start from the whole muscle. We load directly the GO analysis.

if (file.exists("clusterCompare_total_muscle.RData")) {
  load("clusterCompare_total_muscle.RData")
} else {
  warning("Comparison must be there. Check GO analysis")
}

store = "manual-filtered-and-annotated-terms"
keeped_go <- list()
selected_8GO_to_genes <- list()

library(clusterProfiler)

cell="tot"
if (!file.exists(paste0(store, "/","reduced_8GO_clusters.RData"))) {
  
  if (!file.exists(paste0(store,"/total_GO_clusters.RData"))) {
    warning("File need to be created from go analysis")
  } else {
    load(paste0(store,"/total_GO_clusters.RData"))
  }
  
  chosen_classes <- c("apoptosis",
  "phagocytosis",
  "immune response",
  "chemotaxis",
  "angiogenesis",
  "muscle tissue morphogenesis",
  "ECM remodeling",
  "proliferation")
  
  reduced_8GO_Clusters_map <- total_GO_Cluster_map[total_GO_Cluster_map[,2] %in% chosen_classes, , drop=F]
  
  save(reduced_8GO_Clusters_map, file=paste0(store, "/","reduced_8GO_Clusters_map.RData"))
} else {
  load(paste0(store, "/","reduced_8GO_Clusters_map.RData"))
}

reduced_8GO_Clusters_map[reduced_8GO_Clusters_map$meta_cat=="proliferation",]

valid <- reduced_8GO_Clusters_map[,1]
proliferation <- reduced_8GO_Clusters_map[reduced_8GO_Clusters_map[,2]=="proliferation", 1]
sigTable <- as.data.frame(clusterCompare_total_muscle)
sigTable <- sigTable[as.character(sigTable$Description) %in% valid, ,drop=F]

# sigTable[as.character(sigTable$Description) %in% proliferation,1:6]

keep_go <- cellCB:::guess_best_sets_map(reduced_8GO_Clusters_map, sigTable)
keeped_go[[cell]] <- list(keep_go=keep_go, grp=cellCB:::groupByCluster_map(reduced_8GO_Clusters_map))

# Total muscle WT - GO BP

selected_8GO_to_genes[[cell]] <- extract_genes_of_the_8GO_map(sigTable, reduced_8GO_Clusters_map)
p <- cellCB:::personal.dotplot.compareClusterResult(clusterCompare_total_muscle, showCategory=10,
                                               title=paste0(cell, "-wt"),
                                               filterSets = keep_go,
                                               relabel_description = NULL,
                                               relabel_names = invert_map(keep_go)) + 
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10))

res_table <- cellCB:::personal.dotplot.compareClusterResult(clusterCompare_total_muscle,
                                                            showCategory=10,output_df = TRUE,
                                               title=paste0(cell, "-wt"),
                                               filterSets = keep_go,
                                               relabel_description = NULL,
                                               relabel_names = invert_map(keep_go),
                                               )
names_ordered = rev(levels(res_table$Description))
ggsave(p, file="fig-1g-total-muscle-GO-8-category-analysis_review.pdf", height=5, width=7)
write.table(sigTable, file=paste0(tableDir, "/full-", cell, "-wt-enrichGO.txt"), sep="\t", quote=F)

# The GO terms that were founded in the total muscle to be reused in the cell population analysis.
# We need to load the original cluster definition for each cell

cells = c("fap", "ec", "mp", "inf", "per")
activeAtDaysCell <- lapply(cells, function(cell) {
  getActiveAtDaysFromRData(paste0(resultsDir, "/", cell, "-wt.RData"), paste0(cell, "_wt"))
})
names(activeAtDaysCell) <- cells

if (file.exists("cell_wise_GO_enrichment_dfs.RData")) {
  load("cell_wise_GO_enrichment_dfs.RData")
} else {
  warning("File should be there. Check fulla analysis GO")
}

# Now we filterResults by the GO Used in the total tissue.

## We filtered the results using the valid object

cell_wise_GO_enrichment <- lapply(names(cell_wise_GO_enrichment_dfs), function(cell) {
  df <- cell_wise_GO_enrichment_dfs[[cell]]$df
  geneCls <- activeAtDaysCell[[cell]]
  small.df <- df[as.character(df$Description) %in% valid, ,drop=F]
  new("compareClusterResult", compareClusterResult = small.df, geneClusters = geneCls, fun = "my.enrich", .call = match.call())
})
names(cell_wise_GO_enrichment) <- names(cell_wise_GO_enrichment_dfs)

# We collected all the data frame.

cell="fap"

plotList <- list()
cell_sigTable <- as.data.frame(cell_wise_GO_enrichment[[cell]])

keeped_go[[cell]] = sigTable_filter_by_clusters_map(cell_sigTable, reduced_8GO_Clusters_map)
selected_8GO_to_genes[[cell]] <- extract_genes_of_the_8GO_map(cell_sigTable, reduced_8GO_Clusters_map)

p <- cellCB:::personal.dotplot.compareClusterResult(cell_wise_GO_enrichment[[cell]], showCategory=10,
                                               title=paste0(cell, "-wt"),
                                               filterSets = keeped_go[[cell]]$keep_go,
                                               relabel_names = invert_map(keeped_go[[cell]]$keep_go),
                                               name_orders = names_ordered) +
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10))
plotList[[cell]] <- p

cell="ec"
cell_sigTable <- as.data.frame(cell_wise_GO_enrichment[[cell]])

keeped_go[[cell]] = sigTable_filter_by_clusters_map(cell_sigTable, reduced_8GO_Clusters_map)

selected_8GO_to_genes[[cell]] <- extract_genes_of_the_8GO_map(cell_sigTable, reduced_8GO_Clusters_map)

p <- cellCB:::personal.dotplot.compareClusterResult(cell_wise_GO_enrichment[[cell]], showCategory=10,
                                               title=paste0(cell, "-wt"),
                                               filterSets = keeped_go[[cell]]$keep_go,
                                               relabel_names = invert_map(keeped_go[[cell]]$keep_go),
                                               name_orders = names_ordered) +
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10))
plotList[[cell]] <- p

cell="mp"
cell_sigTable <- as.data.frame(cell_wise_GO_enrichment[[cell]])

keeped_go[[cell]] = sigTable_filter_by_clusters_map(cell_sigTable, reduced_8GO_Clusters_map)

selected_8GO_to_genes[[cell]] <- extract_genes_of_the_8GO_map(cell_sigTable, reduced_8GO_Clusters_map)

p <- cellCB:::personal.dotplot.compareClusterResult(cell_wise_GO_enrichment[[cell]], showCategory=10,
                                               title=paste0(cell, "-wt"),
                                               filterSets = keeped_go[[cell]]$keep_go,
                                               relabel_names = invert_map(keeped_go[[cell]]$keep_go),
                                               name_orders = names_ordered) +
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10))
plotList[[cell]] <- p

cell="inf"
cell_sigTable <- as.data.frame(cell_wise_GO_enrichment[[cell]])

keeped_go[[cell]] = sigTable_filter_by_clusters_map(cell_sigTable, reduced_8GO_Clusters_map)

selected_8GO_to_genes[[cell]] <- extract_genes_of_the_8GO_map(cell_sigTable, reduced_8GO_Clusters_map)

p <- cellCB:::personal.dotplot.compareClusterResult(cell_wise_GO_enrichment[[cell]], showCategory=10,
                                               title=paste0(cell, "-wt"),
                                               filterSets = keeped_go[[cell]]$keep_go,
                                               relabel_names = invert_map(keeped_go[[cell]]$keep_go),
                                               name_orders = names_ordered) +
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10))

plotList[[cell]] <- p

cell="per"
cell_sigTable <- as.data.frame(cell_wise_GO_enrichment[[cell]])

keeped_go[[cell]] = sigTable_filter_by_clusters_map(cell_sigTable, reduced_8GO_Clusters_map)

selected_8GO_to_genes[[cell]] <- extract_genes_of_the_8GO_map(cell_sigTable, reduced_8GO_Clusters_map)

p <- cellCB:::personal.dotplot.compareClusterResult(cell_wise_GO_enrichment[[cell]], showCategory=10,
                                               title=paste0(cell, "-wt"),
                                               filterSets = keeped_go[[cell]]$keep_go,
                                               relabel_names = invert_map(keeped_go[[cell]]$keep_go),
                                               name_orders = names_ordered) +
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10))
plotList[[cell]] <- p

# We are going to save the tables of the enrichment. You will find them in the excel file.

for (cell in names(cell_wise_GO_enrichment)){
  table <- as.data.frame(cell_wise_GO_enrichment[[cell]])
  write.table(table, file=paste0(tableDir, "/full-", cell, "-wt-enrichGO.txt"), sep="\t", quote=F)
}

outPlot <- lapply(plotList, function(sp) {
  sp +
    guides(col=guide_legend(title.theme = element_text(size=8),
                            label.theme = element_text(size=8)),
           size=guide_legend(title.theme = element_text(size=8),
                            label.theme = element_text(size=8), order=1))
})

library(cowplot)
plot_grid(plotlist = outPlot, ncol=1)
ggsave(file="supple-all-cells-GO-8-category-analysis-review.pdf", height=20, width=7)

# Now we can move to analyze the KO cell population.
# We need to load the data for the KOs.

cells <- c("fap", "ec", "mp")
activeAtDaysCell.ko <- lapply(cells, function(cell) {
  getActiveAtDaysFromRData(paste0(resultsDir, "/", cell, "-ko.RData"), paste0(cell, "_ko"))
})
names(activeAtDaysCell.ko) <- cells
selected_8GO_to_genes_ko <- list()

# We can enrich. The same setting used for total muscle and WT were used also for KO.

if (file.exists("cell_wise_GO_enrichment_ko_dfs.RData")) {
  load("cell_wise_GO_enrichment_ko_dfs.RData")
} else {
  warning("File should be there. Check fulla GO analysis")
}

cell_wise_GO_enrichment_ko <- lapply(names(cell_wise_GO_enrichment_ko_dfs), function(cell) {
  df <- cell_wise_GO_enrichment_ko_dfs[[cell]]$df
  geneCls <- activeAtDaysCell.ko[[cell]]
  small.df <- df[df$Description %in% valid, ,drop=F]
  new("compareClusterResult", compareClusterResult = small.df, geneClusters = geneCls, fun = "my.enrich", .call = match.call())
})
names(cell_wise_GO_enrichment_ko) <- names(cell_wise_GO_enrichment_ko_dfs)

## Create the data frame.
# We explore the GO BP for FAP KO dataset performing the same plots. Once again with _showCategory=10_ parameter.

cell="fap"

cell_sigTable <- as.data.frame(cell_wise_GO_enrichment_ko[[cell]])
selected_8GO_to_genes_ko[[cell]] <- extract_genes_of_the_8GO_map(cell_sigTable, reduced_8GO_Clusters_map)

cellCB:::personal.dotplot.compareClusterResult(cell_wise_GO_enrichment_ko[[cell]], showCategory=10,
                                               title=paste0(cell, "-ko"),
                                               filterSets = keeped_go[[cell]]$keep_go,
                                               relabel_names = invert_map(keeped_go[[cell]]$keep_go),
                                               name_orders = names_ordered) + 
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10))

# EC KO

cell="ec"

cell_sigTable <- as.data.frame(cell_wise_GO_enrichment_ko[[cell]])
selected_8GO_to_genes_ko[[cell]] <- extract_genes_of_the_8GO_map(cell_sigTable, reduced_8GO_Clusters_map)

cellCB:::personal.dotplot.compareClusterResult(cell_wise_GO_enrichment_ko[[cell]], showCategory=10,
                                               title=paste0(cell, "-ko"),
                                               filterSets = keeped_go[[cell]]$keep_go,
                                               relabel_names = invert_map(keeped_go[[cell]]$keep_go),
                                               name_orders = names_ordered) + 
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10))


# MP KO

cell="mp"

cell_sigTable <- as.data.frame(cell_wise_GO_enrichment_ko[[cell]])
selected_8GO_to_genes_ko[[cell]] <- extract_genes_of_the_8GO_map(cell_sigTable, reduced_8GO_Clusters_map)

cellCB:::personal.dotplot.compareClusterResult(cell_wise_GO_enrichment_ko[[cell]], showCategory=10,
                                               title=paste0(cell, "-ko"),
                                               filterSets = keeped_go[[cell]]$keep_go,
                                               relabel_names = invert_map(keeped_go[[cell]]$keep_go),
                                               name_orders = names_ordered) + 
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10))

# We are going to save the tables of the enrichment also for the KO. You will find them in the excel file.

for (cell in names(cell_wise_GO_enrichment_ko)){
  table <- as.data.frame(cell_wise_GO_enrichment_ko[[cell]])
  write.table(table, file=paste0(tableDir, "/", cell, "-ko-enrichGO.txt"), sep="\t", quote=F)
}

save(selected_8GO_to_genes_ko, selected_8GO_to_genes, file="selected_8GO_to_genes_review.RData")

## Compare the KO and WT
if (file.exists("cell_wise_GO_enrichment_wt_vs_ko_dfs.RData")) {
  load("cell_wise_GO_enrichment_wt_vs_ko_dfs.RData")
} else {
  warning("File should be there. Check GO analysis")
}

cell_wise_GO_enrichment_wt_vs_ko <- lapply(cell_wise_GO_enrichment_wt_vs_ko_dfs, function(clsRes) {
  df <- clsRes$df
  geneCls <- clsRes$geneClusters
  small.df <- df[df$Description %in% valid, ,drop=F]
  new("compareClusterResult", compareClusterResult = small.df, geneClusters = geneCls,
      fun = clsRes$fun, .call = clsRes$.call)
})

cell="ec"
cellCB:::personal.dotplot.compareClusterResult(cell_wise_GO_enrichment_wt_vs_ko[[cell]], x=~day,
                                               title=paste0(cell), 
                                               filterSets = keeped_go[[cell]]$keep_go, 
                                               relabel_names = invert_map(keeped_go[[cell]]$keep_go)) +
  ggplot2::facet_grid(~condition) +
  theme(axis.text.x = element_text(size=10, angle=80, hjust=1), axis.text.y = element_text(size=8))

cell="fap"
cellCB:::personal.dotplot.compareClusterResult(cell_wise_GO_enrichment_wt_vs_ko[[cell]], x=~day,
                                               title=paste0(cell), 
                                               filterSets = keeped_go[[cell]]$keep_go, 
                                               relabel_names = invert_map(keeped_go[[cell]]$keep_go)) +
  ggplot2::facet_grid(~condition) +
  theme(axis.text.x = element_text(size=10, angle=80, hjust=1), axis.text.y = element_text(size=8))

cell="mp"
cellCB:::personal.dotplot.compareClusterResult(cell_wise_GO_enrichment_wt_vs_ko[[cell]], x=~day,
                                               title=paste0(cell), 
                                               filterSets = keeped_go[[cell]]$keep_go, 
                                               relabel_names = invert_map(keeped_go[[cell]]$keep_go)) +
  ggplot2::facet_grid(~condition) +
  theme(axis.text.x = element_text(size=10, angle=80, hjust=1), axis.text.y = element_text(size=8))

# ec-wt-vs-ko-enrichGO

for (cell in names(cell_wise_GO_enrichment_wt_vs_ko)){
  table <- as.data.frame(cell_wise_GO_enrichment_wt_vs_ko[[cell]])
  write.table(table, file=paste0(tableDir, "/", cell, "-wt-vs-ko-auto-paracrine-enrich-go.txt"), sep="\t", quote=F)
}
