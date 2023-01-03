# Load the libraries

library(NbClust)
library(factoextra)
library(clValid)
library(ggplot2)

source("cellRep_pca-clustering-functions.R")


# Setting up the remote dir.

rdatadir="../preprocessing-data/rpkm-RDatas-padj0.01-lfc2/"
save_dir = "clt_eval"
dir.create(save_dir,showWarnings = F)

## Load tot
cell="tot"
condition="wt"
tot_wt <- loadCellRData(paste0(rdatadir, cell,"-", condition, ".RData"), paste0(cell, "_", condition))
tot_wt_pcs <- get_pcs(tot_wt)

# We can use Nbclust Silhouette or Biological Homogeneity Index (clClust)

clust_num_tot_wt <- NbClust(data=tot_wt_pcs[,1:2], method="kmeans", index = "silhouette", max.nc = 10)

tot_wt_silh <- silhuette_plot(clust_num_tot_wt, ylim=c(0,0.8))
tot_wt_bhi <- evaluate_BHI(tot_wt_pcs, 2:10, ylim=c(0,0.8))

tot_wt_silh$plot <- tot_wt_silh$plot + format_axis()
tot_wt_bhi$plot <- tot_wt_bhi$plot + format_axis()

set.seed(1)
tot_wt_km <- kmeans(tot_wt_pcs[,1:2], tot_wt_bhi$best, 100)

tot_wt_pcs_data <- data.frame(tot_wt_pcs, cluster_silh=factor(clust_num_tot_wt$Best.partition))
tot_wt_pcs_data$cluster_bhi <- factor(tot_wt_km$cluster)

plot_2index_clusters(tot_wt_pcs_data, tot_wt_silh, tot_wt_bhi, title=NULL, choose = "bhi")
ggsave(file=paste0(save_dir, "/fig-S1C-tot_wt_clustering_evaluation.pdf"), height = 5, width = 5)
