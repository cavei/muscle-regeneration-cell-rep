# Evaluate LR specificity

# We created the cell specific color palette.

library(RColorBrewer)
colors <- brewer.pal(5, "Set1")
names(colors) <- c("ec","fap","inf","mp","per")
colors.ko <- brewer.pal(5, "Set3")[c(4,5,3)]
names(colors.ko) <- c("ec","fap","mp")

# We load the table were all the valid ligand receptor interaction has been collapsed for the cell that express the receptors. In this case we lost the paracrine autocrine division and also the temporal patterns.
# We can see the specificity of each couple of ligand receptor founded in at least one time point between 2 cell.

receptor_table <- read.table("all-valid-ligand-receptor-interaction-with-cell-wt-collapsed.txt", sep="\t",
           quote="\"", stringsAsFactors = F, check.names = F)

## Pie chart for all wt lig-receptors: receptor side

source("cellRep_pie-chart-template.R")
all_plots <- list()

stats <- table(receptor_table$V3)/nrow(receptor_table)*100
g_idx <- grep(";", names(stats))
gross_labels <- names(stats)
gross_labels[g_idx] <- "shared"
df <- data.frame(label=names(stats), gross_labels=gross_labels, perc=as.numeric(stats), stringsAsFactors = F)

all_plots[[1]] <- my_piePlot(df) + scale_fill_manual(values=c(colors, shared="#b8bbc2")) +
  theme(legend.position = "bottom") + 
  guides(col=guide_legend(title.theme = element_blank(),
                                 label.theme = element_text(size=6),
                                 ncol=3))

## Pie chart for all wt lig-receptors: ligand side

tableLig <- read.table("all-valid-ligand-receptor-interaction-with-cell-wt-ligand-side-collapsed.txt", sep="\t", quote="\"", stringsAsFactors = F, check.names = F)

stats <- table(tableLig$V3)/nrow(tableLig)*100
g_idx <- grep(";", names(stats))
gross_labels <- names(stats)
gross_labels[g_idx] <- "shared"
df <- data.frame(label=names(stats), gross_labels=gross_labels, perc=as.numeric(stats), stringsAsFactors = F)
head(df)
all_plots[[2]] <- my_piePlot(df) + scale_fill_manual(values=c(colors, shared="#b8bbc2")) +
  theme(legend.position = "bottom") + 
  guides(col=guide_legend(title.theme = element_blank(),
                                 label.theme = element_text(size=6),
                                 ncol=3))
all_plots[[2]]

## Pie chart for all ko lig-receptors: receptor side
receptor_table.ko <- read.table("all-valid-ligand-receptor-interaction-with-cell-ko-collapsed.txt", sep="\t",
           quote="\"", stringsAsFactors = F, check.names = F)

source("pie-chart-template.R")
stats <- table(receptor_table.ko$V3)/nrow(receptor_table.ko)*100
g_idx <- grep(";", names(stats))
gross_labels <- names(stats)
gross_labels[g_idx] <- "shared"
df <- data.frame(label=names(stats), gross_labels=gross_labels, perc=as.numeric(stats), stringsAsFactors = F)
head(df)
my_piePlot(df) + scale_fill_manual(values=c(colors.ko, shared="#b8bbc2"))

## Pie chart for all ko lig-receptors: ligand side

tableLig.ko <- read.table("all-valid-ligand-receptor-interaction-with-cell-ko-ligand-side-collapsed.txt", sep="\t",
           quote="\"", stringsAsFactors = F, check.names = F)

stats <- table(tableLig.ko$V3)/nrow(tableLig.ko)*100
g_idx <- grep(";", names(stats))
gross_labels <- names(stats)
gross_labels[g_idx] <- "shared"
df <- data.frame(label=names(stats), gross_labels=gross_labels, perc=as.numeric(stats), stringsAsFactors = F)
head(df)
my_piePlot(df) + scale_fill_manual(values=c(colors.ko, shared="#b8bbc2")) 

## Pie chart for angiogenesis related wt lig-receptors: receptor side

load("metaGO_class_related_gene.RData")
source("pie-chart-template.R")
angio_receptor_table <- receptor_table[receptor_table$V2 %in% metaGO_class_related_gene$angiogenesis, ]
stats <- table(angio_receptor_table$V3)/nrow(angio_receptor_table)*100

g_idx <- grep(";", names(stats))
gross_labels <- names(stats)
gross_labels[g_idx] <- "shared"
df <- data.frame(label=names(stats), gross_labels=gross_labels, perc=as.numeric(stats), stringsAsFactors = F)

all_plots[[3]] <- my_piePlot(df) + scale_fill_manual(values=c(colors, shared="#b8bbc2"))+
  theme(legend.position = "bottom") + 
  guides(col=guide_legend(title.theme = element_blank(),
                                 label.theme = element_text(size=6),
                                 ncol=3))

## Pie chart for angiogenesis related wt lig-receptors: ligand side

angio_ligand_table <- tableLig[tableLig$V1 %in% metaGO_class_related_gene$angiogenesis, ]
stats <- table(angio_ligand_table$V3)/nrow(angio_ligand_table)*100

g_idx <- grep(";", names(stats))
gross_labels <- names(stats)
gross_labels[g_idx] <- "shared"
df <- data.frame(label=names(stats), gross_labels=gross_labels, perc=as.numeric(stats), stringsAsFactors = F)
all_plots[[4]] <- my_piePlot(df) + scale_fill_manual(values=c(colors, shared="#b8bbc2"))+
  theme(legend.position = "bottom") + 
  guides(col=guide_legend(title.theme = element_blank(),
                                 label.theme = element_text(size=6),
                                 ncol=3))
## Paper Figure
library(cowplot)
legend <- get_legend(all_plots[[1]])
all_p <- lapply(all_plots, function(plot) {
  plot + theme(legend.position="none")
})
main <- plot_grid(plotlist = all_p, ncol=2)
plot_grid(main, legend, ncol = 1, rel_heights = c(1,0.1))
ggsave(file="fig2b-pie-chart-lig-reg-specificity.pdf", height = 4.2, width = 4)

## Pie chart for immune response related wt lig-receptors: receptor side

load("metaGO_class_related_gene.RData")
source("pie-chart-template.R")
immune_receptor_table <- receptor_table[receptor_table$V2 %in% metaGO_class_related_gene$`immune response`, ]
stats <- table(immune_receptor_table$V3)/nrow(immune_receptor_table)*100

g_idx <- grep(";", names(stats))
gross_labels <- names(stats)
gross_labels[g_idx] <- "shared"
df <- data.frame(label=names(stats), gross_labels=gross_labels, perc=as.numeric(stats), stringsAsFactors = F)
my_piePlot(df) + scale_fill_manual(values=c(colors, shared="#b8bbc2"))

## Pie chart for immune response related wt lig-receptors: ligand side

immune_ligand_table <- tableLig[tableLig$V1 %in% metaGO_class_related_gene$`immune response`, ]
stats <- table(immune_ligand_table$V3)/nrow(immune_ligand_table)*100

g_idx <- grep(";", names(stats))
gross_labels <- names(stats)
gross_labels[g_idx] <- "shared"
df <- data.frame(label=names(stats), gross_labels=gross_labels, perc=as.numeric(stats), stringsAsFactors = F)
my_piePlot(df) + scale_fill_manual(values=c(colors, shared="#b8bbc2"))
