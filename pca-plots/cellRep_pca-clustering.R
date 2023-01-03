# Load the libraries

library(NbClust)
library(factoextra)
library(clValid)
library(ggplot2)

source("cellRep_pca-clustering-functions.R")

rdatadir="../preprocessing-data/rpkm-RDatas-padj0.01-lfc2/"
save_dir = "clt_eval"
dir.create(save_dir,showWarnings = F)

# FAP
cell="fap"
condition="wt"
fap_wt <- loadCellRData(paste0(rdatadir, cell,"-", condition, ".RData"), paste0(cell, "_", condition))
fap_wt_pcs <- get_pcs(fap_wt)

# We use Nbclust Silhouette or Biological Homogeneity Index (clClust)

clust_num_fap_wt <- NbClust(data=fap_wt_pcs[,1:2], method="kmeans", index = "silhouette", max.nc = 10)
fap_wt_silh <- silhuette_plot(clust_num_fap_wt, ylim=c(0,0.8))
fap_wt_bhi <- evaluate_BHI(fap_wt_pcs, 2:10, ylim=c(0,0.8))
set.seed(1)
fap_wt_km <- kmeans(fap_wt_pcs[,1:2], fap_wt_bhi$best, 100)

fap_wt_pcs_data <- data.frame(fap_wt_pcs, cluster_silh=factor(clust_num_fap_wt$Best.partition))
fap_wt_pcs_data$cluster_bhi <- factor(fap_wt_km$cluster)

fap_wt_silh$plot <- fap_wt_silh$plot + format_axis()
fap_wt_bhi$plot <- fap_wt_bhi$plot + format_axis()

plot_2index_clusters(fap_wt_pcs_data, fap_wt_silh, fap_wt_bhi, title=NULL, choose = "silhuette")
ggsave(file=paste0(save_dir, "/fap_wt_clustering_evaluation.pdf"), height = 5, width = 5)


## FAP ko
cell="fap"
condition="ko"
fap_ko <- loadCellRData(paste0(rdatadir, cell,"-", condition, ".RData"), paste0(cell, "_", condition))
fap_ko_pcs <- get_pcs(fap_ko)

clust_num_fap_ko <- NbClust(data=fap_ko_pcs[,1:2], method="kmeans", index = "silhouette", max.nc = 10)

fap_ko_silh <- silhuette_plot(clust_num_fap_ko, ylim=c(0,0.8))
fap_ko_bhi <- evaluate_BHI(fap_ko_pcs, 2:10, ylim=c(0,0.8))

set.seed(1)
fap_ko_km <- kmeans(fap_ko_pcs[,1:2], fap_ko_bhi$best, 100)

fap_ko_pcs_data <- data.frame(fap_ko_pcs, cluster_silh=factor(clust_num_fap_ko$Best.partition))
fap_ko_pcs_data$cluster_bhi <- factor(fap_ko_km$cluster)

fap_ko_silh$plot <- fap_ko_silh$plot + format_axis()
fap_ko_bhi$plot <- fap_ko_bhi$plot + format_axis()

# Use K=4 as for WT 
fap_ko_km <- kmeans(fap_ko_pcs[,1:2], 4, 100)
fap_ko_pcs_data$cluster_bhi <- factor(fap_ko_km$cluster)

plot_2index_clusters(fap_ko_pcs_data, fap_ko_silh, fap_ko_bhi, title=NULL, choose = "bhi")
ggsave(file=paste0(save_dir, "/fap_ko_clustering_evaluation.pdf"), height = 5, width = 5)

# EC
cell="ec"
condition="wt"
ec_wt <- loadCellRData(paste0(rdatadir,
                               cell,"-", condition, ".RData"), paste0(cell, "_", condition))
ec_wt_pcs <- get_pcs(ec_wt)

clust_num_ec_wt <- NbClust(data=ec_wt_pcs[,1:2], method="kmeans", index = "silhouette", max.nc = 10)

ec_wt_silh <- silhuette_plot(clust_num_ec_wt,  ylim=c(0,0.8))
ec_wt_bhi <- evaluate_BHI(ec_wt_pcs, 2:10,  ylim=c(0,0.8))

set.seed(1)
ec_wt_km <- kmeans(ec_wt_pcs[,1:2], ec_wt_bhi$best, 100)

ec_wt_silh$plot <- ec_wt_silh$plot + format_axis()
ec_wt_bhi$plot <- ec_wt_bhi$plot + format_axis()

ec_wt_pcs_data <- data.frame(ec_wt_pcs, cluster_silh=factor(clust_num_ec_wt$Best.partition))
ec_wt_pcs_data$cluster_bhi <- factor(ec_wt_km$cluster)

plot_2index_clusters(ec_wt_pcs_data, ec_wt_silh, ec_wt_bhi, title = NULL, choose = "bhi")
ggsave(file=paste0(save_dir, "/fig-1d-ec_wt_clustering_evaluation.pdf"), height = 5, width = 5)

#EC KO
cell="ec"
condition="ko"
ec_ko <- loadCellRData(paste0(rdatadir, cell,"-", condition, ".RData"), paste0(cell, "_", condition))
ec_ko_pcs <- get_pcs(ec_ko)
clust_num_ec_ko <- NbClust(data=ec_ko_pcs[,1:2], method="kmeans", index = "silhouette", max.nc = 10)

ec_ko_silh <- silhuette_plot(clust_num_ec_ko, ylim=c(0,0.8))
ec_ko_bhi <- evaluate_BHI(ec_ko_pcs, 2:10,  ylim=c(0,0.8))

ec_ko_silh$plot <- ec_ko_silh$plot + format_axis()
ec_ko_bhi$plot <- ec_ko_bhi$plot + format_axis()

set.seed(1)
ec_ko_km <- kmeans(ec_ko_pcs[,1:2], ec_ko_bhi$best, 100)

ec_ko_pcs_data <- data.frame(ec_ko_pcs, cluster_silh=factor(clust_num_ec_ko$Best.partition))
ec_ko_pcs_data$cluster_bhi <- factor(ec_ko_km$cluster)

plot_2index_clusters(ec_ko_pcs_data, ec_ko_silh, ec_ko_bhi, title = NULL, choose = "bhi")
ggsave(file=paste0(save_dir, "/fig-1e-ec_ko_clustering_evaluation.pdf"), height = 5, width = 5)


# Inflammatory 

cell="inf"
condition="wt"
inf_wt <- loadCellRData(paste0(rdatadir,
                               cell,"-", condition, ".RData"), paste0(cell, "_", condition))
inf_wt_pcs <- get_pcs(inf_wt)

clust_num_inf_wt <- NbClust(data=inf_wt_pcs[,1:2], method="kmeans", index = "silhouette", max.nc = 10)

inf_wt_silh <- silhuette_plot(clust_num_inf_wt,  ylim=c(0,0.8))
inf_wt_bhi <- evaluate_BHI(inf_wt_pcs, 2:10, ylim=c(0,0.8))

set.seed(1)
inf_wt_km <- kmeans(inf_wt_pcs[,1:2], inf_wt_bhi$best, 100)

inf_wt_pcs_data <- data.frame(inf_wt_pcs, cluster_silh=factor(clust_num_inf_wt$Best.partition))
inf_wt_pcs_data$cluster_bhi <- factor(inf_wt_km$cluster)

plot_2index_clusters(inf_wt_pcs_data, inf_wt_silh, inf_wt_bhi, title = "INF WT")
plot_2index_clusters(inf_wt_pcs_data, inf_wt_silh, inf_wt_bhi, title = "INF WT", choose = "silhuette")
ggsave(file="inf_wt_clustering_evaluation.pdf", height = 7, width = 9)

# MP
cell="mp"
condition="wt"
mp_wt <- loadCellRData(paste0(rdatadir,
                               cell,"-", condition, ".RData"), paste0(cell, "_", condition))
mp_wt_pcs <- get_pcs(mp_wt)


clust_num_mp_wt <- NbClust(data=mp_wt_pcs[,1:2], method="kmeans", index = "silhouette", max.nc = 10)

mp_wt_silh <- silhuette_plot(clust_num_mp_wt, ylim=c(0,0.8))
mp_wt_bhi <- evaluate_BHI(mp_wt_pcs, 2:10, ylim=c(0,0.8))

mp_wt_silh$plot <- mp_wt_silh$plot + format_axis()
mp_wt_bhi$plot <- mp_wt_bhi$plot + format_axis()

set.seed(1)
mp_wt_km <- kmeans(mp_wt_pcs[,1:2], mp_wt_bhi$best, 100)

mp_wt_pcs_data <- data.frame(mp_wt_pcs, cluster_silh=factor(clust_num_mp_wt$Best.partition))
mp_wt_pcs_data$cluster_bhi <- factor(mp_wt_km$cluster)

plot_2index_clusters(mp_wt_pcs_data, mp_wt_silh, mp_wt_bhi, title = NULL)

ggsave(file=paste0(save_dir, "/mp_wt_clustering_evaluation.pdf"), height = 5, width = 5)

# MP KO
cell="mp"
condition="ko"
mp_ko <- loadCellRData(paste0(rdatadir,
                               cell,"-", condition, ".RData"), paste0(cell, "_", condition))
mp_ko_pcs <- get_pcs(mp_ko)

clust_num_mp_ko <- NbClust(data=mp_ko_pcs[,1:2], method="kmeans", index = "silhouette", max.nc = 6)

mp_ko_silh <- silhuette_plot(clust_num_mp_ko, ylim=c(0,0.8))
mp_ko_bhi <- evaluate_BHI(mp_ko_pcs, 2:6, ylim=c(0,0.8))

set.seed(1)
mp_ko_km <- kmeans(mp_ko_pcs[,1:2], mp_ko_bhi$best, 100)

mp_ko_pcs_data <- data.frame(mp_ko_pcs, cluster_silh=factor(clust_num_mp_ko$Best.partition))
mp_ko_pcs_data$cluster_bhi <- factor(mp_ko_km$cluster)

plot_2index_clusters(mp_ko_pcs_data, mp_ko_silh, mp_ko_bhi, title = "mp ko")
ggsave(file="mp_ko_clustering_evaluation.pdf", height = 9, width = 7)
