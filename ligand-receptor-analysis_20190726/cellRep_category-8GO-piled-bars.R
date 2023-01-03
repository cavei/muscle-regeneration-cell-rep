# Summarize Autocrine / Paracrine Contributions

library(org.Mm.eg.db)
library(clusterProfiler)
library(ggplot2)
library(cellCB)

all_enrichments <- read.table("all-cell-enrichment.net.txt", header=T, sep="\t", check.names = F, stringsAsFactors = F)
all_enrichments$lid <- NULL


load("../go_analysis_20190726/manual-filtered-and-annotated-terms/total_GO_clusters.RData")
selected_all <- total_GO_Cluster_map[total_GO_Cluster_map[,2] == "immune response", 1]

load("../go_analysis_20190726//manual-filtered-and-annotated-terms/reduced_8GO_Clusters_map.RData")
selected_all_8 <- reduced_8GO_Clusters_map[reduced_8GO_Clusters_map[,2]=="immune response", 1]
identical(selected_all, selected_all_8)

all_enrichments_8GO <- all_enrichments[all_enrichments$Description %in% reduced_8GO_Clusters_map[,1], , drop=F]

extract_metaGO_class_related_gene <- function(all_enrichments, map) {
  rmap <- map[map[,1] %in% all_enrichments$Description, ,drop=F]
  cls <- cellCB:::groupByCluster_map(rmap)
  lapply(cls, function(gos) {
    class_specific <- all_enrichments[all_enrichments$Description %in% gos, , drop=F]
    class_specific_genes <- unique(unlist(strsplit(class_specific$geneID, "/")))
  })
}

metaGO_class_related_gene_8GO <- extract_metaGO_class_related_gene(all_enrichments_8GO, reduced_8GO_Clusters_map)
save(metaGO_class_related_gene_8GO, file="metaGO_class_related_gene_8GO.RData")

all_nets <- read.table("all-cell-network-wt.net.txt", header=T,
                       sep="\t", check.names = F, stringsAsFactors = F)
colnames(all_nets) <- c("ligand","source","cell_l","day","receptor","cell_r")

create_ligand_spec_db <- function(net) {
  feed <- net$feed
  id <- net$id
  dict <- tapply(feed, id, function(x) paste(sort(unique(x)), collapse=";"))
  dict
}

metaGO_class_8GO_related_net <- lapply(metaGO_class_related_gene_8GO, function(class_related) {
  selection <- all_nets$ligand %in% class_related & all_nets$receptor %in% class_related
  class_net <- all_nets[selection, , drop=F]
  class_net$feed <- ifelse(class_net$cell_l == class_net$cell_r, "autocrine", "paracrine")
  class_net$id <- paste(class_net$ligand, class_net$receptor, class_net$day, class_net$cell_r, sep="_")
  
  class_dict <- create_ligand_spec_db(class_net)
  
  ord <- paste(class_net$ligand,class_net$receptor,class_net$day, class_net$cell_r, sep="_")
  class_net$feed <- class_dict[ord]
  class_net
})


cells <- c("ec", "fap", "mp", "inf", "per")
auto_para_dataframes <- lapply(cells, function(cell) {
  cell_class_net <- lapply(names(metaGO_class_8GO_related_net), function(name) {
    class_net<- metaGO_class_8GO_related_net[[name]]
    cell_receptor_net<-class_net[class_net$cell_r==cell, ,drop=F]
    
    pairs <- unique(apply(as.matrix(cell_receptor_net[, c("ligand", "receptor", "feed"), drop=F]), 1, paste, collapse="_"))
    if (length(pairs) == 0)
      return(NULL)
    feed_system <- sapply(strsplit(pairs, "_"), function(x) x[3])
    out <- rep(0,3)
    names(out) <- c( "autocrine", "autocrine;paracrine", "paracrine")
    out_t <- table(feed_system)
    out[names(out_t)] <-out_t 
    out
  })
  names(cell_class_net) <- names(metaGO_class_8GO_related_net)
  
  df_bars <- do.call(rbind, cell_class_net)
  
  keep_class <- row.names(df_bars)[keep <- apply(df_bars, 1, sum) >= 0] # Do not filter
  df_bars <- data.frame(class = row.names(df_bars),
                        autocrine=df_bars[,1], shared=df_bars[,2], paracrine=df_bars[,3],
               stringsAsFactors = F)
  library(tidyr)
  df_bars=tidyr::gather(df_bars, "feed", "counts", -c("class"))
  
  list(df=df_bars, meaningfull_class=keep_class)
})
names(auto_para_dataframes) <- cells

keep_this_classes <- unique(unlist(lapply(auto_para_dataframes, function(x) {
       x$meaningfull_class
})))
all_bars <- do.call(rbind, lapply(auto_para_dataframes, '[[', 1))
meaningfull_ylim <- c(min(all_bars$counts), max(all_bars$counts))

plots <- lapply(names(auto_para_dataframes), function(cell){
  x <- auto_para_dataframes[[cell]]
  df_bars <- x$df
  # select <- df_bars$class %in% keep_this_classes
  missing <- keep_this_classes[! keep_this_classes %in% df_bars$class]
  missing_df=NULL
  if (length(missing)){
    missing_df <- lapply(missing, function(class) {
      data.frame(class=class, feed=c("autocrine", "paracrine"), counts=0, stringsAsFactors = F)
    })
    missing_df <- do.call(rbind, missing_df)
  }
  
  df_bars <- rbind(df_bars, missing_df)
  df_bars$class <- sapply(strsplit(df_bars$class, " "), paste, collapse="\n")
  p <- ggplot(data=df_bars, aes(x=class, y=counts, fill=feed)) +
      geom_bar(stat="identity") +
      theme_minimal() + scale_fill_manual(values=c('#ffb142',"#33d9b2", "#218c74")) +
      ggtitle(paste0(cell, " - wt")) + 
      theme(axis.text.x = element_text(size=8, angle = 90, vjust = 0.5),
            axis.text.y = element_text(size=8), axis.title.x = element_blank()) +
    guides(fill=guide_legend(title.theme = element_blank(),
                                 label.theme = element_text(size=8)))
    
  p
})

names(plots) <- cells
dir.create("8GO_categories_plot", showWarnings = F)
for (cell in names(plots)) {
  ggsave(filename = paste0("8GO_categories_plot/", cell, "-lig-rec-8GO-piled.pdf"), plots[[cell]], height=5)
}

# Create small plots

extrafont::loadfonts()
splots <- lapply(names(auto_para_dataframes), function(cell){
  x <- auto_para_dataframes[[cell]]
  df_bars <- x$df
  # select <- df_bars$class %in% keep_this_classes
  missing <- keep_this_classes[! keep_this_classes %in% df_bars$class]
  missing_df=NULL
  if (length(missing)){
    missing_df <- lapply(missing, function(class) {
      data.frame(class=class, feed=c("autocrine", "paracrine"), counts=0, stringsAsFactors = F)
    })
    missing_df <- do.call(rbind, missing_df)
  }
  
  df_bars <- rbind(df_bars, missing_df)
  df_bars$class <- sapply(strsplit(df_bars$class, " "), paste, collapse="\n")
  
  p <- ggplot(data=df_bars, aes(x=class, y=counts, fill=feed)) +
      geom_bar(stat="identity") +
      theme_minimal() + scale_fill_manual(values=c('#ffb142',"#33d9b2", "#218c74")) +
      ggtitle(paste0(cell, " - wt")) + 
      theme(axis.text.x = element_text(size=8, angle = 90, vjust = 0.5),
            axis.text.y = element_text(size=8)) +
      theme(axis.text = element_text(family = "Arial"),
            axis.title = element_text(size=8),
            plot.title = element_text(size=8),
            legend.position = "none") +
      theme(axis.title = element_blank()) +
      theme(axis.line.x = element_blank(), axis.text.x = element_blank(),
            axis.ticks.length = unit(0,"pt")) 
      
  p
})
splots[[1]] 
splots[[3]] <- splots[[3]] + ylim(0,100)
splots[[5]] <- splots[[5]] + ylim(0,100)

df_bars <- auto_para_dataframes[["ec"]]$df
df_bars$counts <- 0
  # select <- df_bars$class %in% keep_this_classes
missing <- keep_this_classes[! keep_this_classes %in% df_bars$class]
missing_df=NULL
if (length(missing)){
  missing_df <- lapply(missing, function(class) {
    data.frame(class=class, feed=c("autocrine", "paracrine"), counts=0, stringsAsFactors = F)
  })
  missing_df <- do.call(rbind, missing_df)
}

df_bars <- rbind(df_bars, missing_df)
df_bars$class <- sapply(strsplit(df_bars$class, " "), paste, collapse="\n")
  
p <- ggplot(data=df_bars, aes(x=class, y=counts, fill=feed)) +
    geom_bar(stat="identity") +
    theme_minimal() + scale_fill_manual(values=c('#ffb142',"#33d9b2", "#218c74")) +
    ggtitle(paste0(cell, " - wt")) + 
    theme(axis.text.x = element_text(size=8, angle = 90, vjust = 0.5),
          axis.text.y = element_text(size=8)) +
    theme(axis.text = element_text(family = "Arial"),
          axis.title = element_text(size=8),
          plot.title = element_blank(),
          legend.position = "none") +
    theme(axis.title = element_blank()) +
    theme(axis.line.x = element_blank(),
          legend.position = "bottom") + 
  ylim(c(0,100)) + guides(fill=guide_legend(title.theme = element_blank(),
                                 label.theme = element_text(size=6),
                                 ncol=3))
      
p
ggsave(p, filename = paste0("8GO_categories_plot/", "x-axsis-lig-rec-8GO-piled.pdf"), height=3, width = 4)

library(cowplot)
plot_grid(plotlist = splots[c(1,3:5)], nrow=4, ncol=1)

ggsave(filename = paste0("8GO_categories_plot/", "cells-lig-rec-8GO-piled.pdf"), height=6, width = 4)
