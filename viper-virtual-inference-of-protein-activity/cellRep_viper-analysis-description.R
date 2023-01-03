# Viper analysis

library(org.Mm.eg.db)

## Methods

### Protein activity inference with Viper
# To infer protein activity, we applied Viper (Viper R package) analysis to our data. Viper assumes that the expression of transcriptional targets of a protein (regulon) can be a reporter of the protein activity itself.
# Following the approach described in [Ref 1] we created the regulatory networks from the cell population dataset (WT and KO) and the total wt tissue. Both datasets were normalized separately as describe earlier (see Filtering and Normalization Step).
# We set up a candidate regulator list. A regulator that share two categories was assigned following this priority: ligands (Rezza et al.), receptors (Rezza et al.), transcription factors (TF; GO:0003700, GO:0004677, GO:0030528, GO:0004677), transcriptional cofactors (COF; GO:0003712) and signaling pathway related genes (SP: GO:0007165, GO:0005622, GO:0005886). 
# To create the regulatory networks of the regulator list, we run ARACNE with 100 bootstrap and mutual information (MI) p-values threshold 10-8. Cell and total tissue regulatory networks were merged and regulators with a regulon size smaller than 25 co-expressed genes were excluded from the analysis.
# Viper function were used to create the protein activity matrix of the regulators and the cell population dataset.

### TimeClip analysis
# Time clip analysis were performed using the protein activity matrix. For each day – cell combination replicates were averaged. TimeClip was performd on a subset of Reactome pathways that contain the gene Kdr.
# 
# _Ref 1: Functional characterization of somatic mutations in cancer using network-based inference of rotein activity; Alvarez 2016_

# Analysis

## Load data for cells
# First we loaded the eset with all the cells.

load("all-cell-eset.RData")

## Run viper
# To run viper we need the regulon network that we computed with ARACNe. After reading the network we run viper.

library(viper)
if (!file.exists("all_cell_activity.RData")) {
  regulon <- aracne2regulon("ARACNe-AP/tot-all-merged-network-gt24.adj", eset, verbose = FALSE)
  all_cell_activity <- viper(eset, regulon, verbose = FALSE)
  save(all_cell_activity, file="all_cell_activity.RData")
} else {
  load("all_cell_activity.RData")
}

# We can create the protein activity matrix.

activity_inferred <- exprs(all_cell_activity)
activity_inferred_a <- cbind(gene=row.names(activity_inferred), activity_inferred)
write.table(activity_inferred_a, file="all_cell_activity_inferred.txt", sep="\t", quote=F, row.names=F)

# We manually curated a set of receptor / markers of interest. We used the list to extract the corresponding protein activity matrix.

cellFrom <- list(fap=c("Tgfbr1", "Tgfbr2", "Pdgfrb", "Pdgfra", "Tnfrsf14", "Tnfrsf12a", "EphB2", "Itga5"),
                 ec=c("Kdr", "Tek", "Nrp1", "Nrp2", "Notch4", "Flt1", "EphB4", "Il1r1", "Il1r2"),
                 mp=c("Itga7"),
                 inf=c("Ccr5", "Ccr2", "Csf1r", "Ifngr1"),
                 per=c("Itga7", "Notch3", "Pdgfrb"),
                 additional=c("Fgfr3", "Fgfr1", "Fgfr2", "Fgfr4"))

# ligands=c("Vegfa","Vegfc", "Angpt1", "Angpt2","Dll4", "Dll1", "Tgfb2", "Tnf", "Pdgfa")

interestingReceptors <- unlist(cellFrom)
interestingReceptors <- intersect(interestingReceptors, row.names(activity_inferred))
interesting_receptor_activity <- activity_inferred[interestingReceptors, ]

# Now we load a custom made function to average replicates. We averaged both whole matrix and receptor matrix.

source("average_replicates.R")
activity_inferred_avg <- average_replicates(activity_inferred)
activity_inferred_avg <- cbind(gene=row.names(activity_inferred_avg), activity_inferred_avg)
write.table(activity_inferred_avg, file="all_cell_activity_inferred_avergae.txt", sep="\t", quote=F, row.names=F)
activity_inferred_avg <- average_replicates(activity_inferred)

avg_interesting_rec_activity <- average_replicates(interesting_receptor_activity)
ann_columns <- do.call(rbind, strsplit(colnames(avg_interesting_rec_activity), "_"))
ann_columns <- data.frame(row.names=colnames(avg_interesting_rec_activity), cell=ann_columns[,2])

order <- unique(as.character(ann_columns[,1]))
numbers <- table(as.character(ann_columns[,1]))[order]

# Now we set up the cell specific colors that will be used in the annotation.

library(RColorBrewer)
colors <- brewer.pal(5, "Set1")
names(colors) <- c("ec-wt","fap-wt","inf-wt","mp-wt","per-wt")
colors.ko <- brewer.pal(5, "Set3")[c(4,5,3)]
names(colors.ko) <- c("ec-ko","fap-ko","mp-ko")
cls <- c(colors, colors.ko)

# Now we can perform some heatmap.

new_colnames <- sapply(strsplit(colnames(avg_interesting_rec_activity), "_"), function(x) x[1])

pheatmap::pheatmap(avg_interesting_rec_activity, cluster_cols = F, scale="row",
                   annotation_col = ann_columns, 
                   gaps_col = cumsum(numbers),
                   labels_col = new_colnames,
                   annotation_colors = list(cell=cls),
                   fontsize = 10)

pheatmap::pheatmap(avg_interesting_rec_activity, cluster_cols = F, scale="row",
                   annotation_col = ann_columns, 
                   gaps_col = cumsum(numbers),
                   labels_col = new_colnames,
                   annotation_colors = list(cell=cls),
                   fontsize = 8,
                   fontsize_col = 6,
                   filename = "fig-s3c-receptor-activity-scaled-row.pdf", width=8, height=4)

# We now need to zoom in into particular population. First MP.

watch_these <- c("Kdr", "Tek", "Flt1")
intersect(row.names(activity_inferred_avg), watch_these)
mpIdx <-grep("mp-", colnames(avg_interesting_rec_activity))

mp_ann <- droplevels(ann_columns[mpIdx, ,drop=F])
mp_cls <- cls[c("mp-wt", "mp-ko")]

pheatmap::pheatmap(avg_interesting_rec_activity[watch_these, mpIdx], cluster_cols = F, scale="row",
                   annotation_col = mp_ann, 
                   gaps_col = c(9),
                   # gaps_col = cumsum(numbers), 
                   annotation_colors = list(cell=mp_cls))

# What about a scatter plot.

source("functions.R")
library(ggplot2)

scaled_avg_interesting_rec_activity <- t(scale(t(avg_interesting_rec_activity[, mpIdx]), center=T, scale=T))

# choose between scaled or not
# exp <- scaled_avg_interesting_rec_activity["Kdr", ]
exp <- avg_interesting_rec_activity["Kdr", mpIdx]
days <- extract_days(names(exp))
cell <- extract_cell(names(exp))

data <- data.frame(Kdr=exp, day=days, cell=factor(cell, levels = names(mp_cls)), row.names=names(exp))

kdr_mp <- ggplot(data = data, aes(x = day, y = Kdr)) +
  geom_point(pch = 1, color= mp_cls[data$cell]) +
  geom_smooth(aes(colour=cell), formula=y ~ x, se=F) +
  scale_color_manual(values=mp_cls) + 
  ylab("Activity score") + 
  theme_classic() + 
  scale_x_continuous(breaks = c(unique(days),14)) +
  theme(axis.text = element_text(size=8),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "top") +
  guides(col=guide_legend(title = element_blank(), label.theme = element_text(size=6)))
  
# Now its time to watch closer EC populations.

watch_these <- c("Kdr", "Tek", "Flt1")
ecIdx <-grep("ec-", colnames(avg_interesting_rec_activity))

ec_ann <- droplevels(ann_columns[ecIdx, ,drop=F])
ec_cls <- cls[c("ec-wt", "ec-ko")]

pheatmap::pheatmap(avg_interesting_rec_activity[watch_these, ecIdx], cluster_cols = F, scale="row",
                   annotation_col = ec_ann, 
                   gaps_col = c(10),
                   # gaps_col = cumsum(numbers), 
                   annotation_colors = list(cell=ec_cls))

source("functions.R")
library(ggplot2)
scaled_avg_interesting_rec_activity <- t(scale(t(avg_interesting_rec_activity[, ecIdx]), center=T, scale=T))

# exp <- scaled_avg_interesting_rec_activity["Kdr", ]
exp <- avg_interesting_rec_activity["Kdr", ecIdx]
days <- extract_days(names(exp))
cell <- extract_cell(names(exp))

data <- data.frame(Kdr=exp, day=days, cell=factor(cell, levels = names(ec_cls)), row.names=names(exp))

kdr_ec <- ggplot(data = data, aes(x = day, y = Kdr)) +
  geom_point(pch = 1, color= ec_cls[data$cell]) +
  geom_smooth(aes(colour=cell), formula=y ~ x, se=F) +
  scale_color_manual(values=ec_cls) + 
  ylab("Kdr scaled activity") +
  theme_classic() + 
  scale_x_continuous(breaks = unique(days))+
  theme(axis.text = element_text(size=8),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "top") +
  guides(col=guide_legend(title = element_blank(), label.theme = element_text(size=6)))

library(cowplot)
library(ggplot2)
library(gridExtra)
library(grid)
y.grob <- textGrob("Activity score", gp=gpar(fontface="bold", fontsize=8), rot=90)
x.grob <- textGrob("day", gp=gpar(fontface="bold", fontsize=8))

fig.s4a <- plot_grid(kdr_ec, kdr_mp, rel_widths = c(14,10))

grid.arrange(arrangeGrob(fig.s4a, left = y.grob, bottom = x.grob))

ggsave(filename = "fig-s4a-kdr-profile-in-ec-and-mp.pdf", width=8, height=4)

library(graphite)
load("../timeClip-analysis_20190726/reactome.RData")

genes <- gsub(pattern = "SYMBOL:", "", nodes(reactome[["Signaling by VEGF"]]))
genes <- intersect(genes, row.names(activity_inferred_avg))

ecIdx <-grep("ec-", colnames(activity_inferred_avg))

ec_ann <- droplevels(ann_columns[ecIdx, ,drop=F])
ec_cls <- cls[c("ec-wt", "ec-ko")]


pheatmap::pheatmap(activity_inferred_avg[genes, ], cluster_cols = F, scale="row",
                   annotation_col = ann_columns, 
                   # gaps_col = cumsum(numbers), 
                   annotation_colors = list(cell=cls))

pheatmap::pheatmap(activity_inferred_avg[genes, ecIdx], cluster_cols = F,
                   scale = "row",
                   annotation_col = ec_ann, 
                   gaps_col = 10,
                   annotation_colors = list(cell=ec_cls))
