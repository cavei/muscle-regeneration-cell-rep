# timeClip analysis
library(timeClip)

# We have store all our data in the dirctory "viper-virtual-inference-of-protein-activity".
dataDirname <- "../viper-virtual-inference-of-protein-activity/"

file <- "all_cell_activity_inferred_avergae.txt"
activity <- read.table(file=paste0(dataDirname,file), sep="\t", header=T, row.names=1, quote="\"", check.names = F, stringsAsFactors = F)
pheno_mat <- do.call(rbind, strsplit(colnames(activity), "_"))

pheno_mat <- data.frame(row.names=colnames(activity),
                        dday= pheno_mat[,1],
                        day=as.numeric(gsub("d", "", pheno_mat[,1])),
                        cell_cond = pheno_mat[,2])

# Now we need the set of pathways from graphite. We used both KEGG and REACTOME databases. We keep only pathways with more than 10 genes.

library(graphite)
if (file.exists("reactome.RData")) {
  load("reactome.RData")
} else {
  reactome <- pathways(species = "mmusculus", database = "reactome")
  reactome <- convertIdentifiers(reactome, "symbol")
  save(reactome, file="reactome.RData")
}
if (file.exists("kegg.RData")) {
  load("kegg.RData")
} else {
  kegg <- pathways(species = "mmusculus", database = "kegg")
  kegg <- convertIdentifiers(kegg, "symbol")
  save(kegg, file="kegg.RData")
}

nodeLength <- sapply(reactome, function(g) length(graphite::nodes(g)))
reactomeSmall <- reactome[names(reactome)[nodeLength>=10]]

nodeLength <- sapply(kegg, function(g) length(graphite::nodes(g)))
keggSmall <- kegg[names(kegg)[nodeLength>=10]]

kegg_pathways_with_Kdr <- names(which(sapply(lapply(keggSmall, nodes), function(nds) "SYMBOL:Kdr" %in% nds)))
reactome_pathways_with_Kdr <- names(which(sapply(lapply(reactomeSmall, nodes), function(nds) "SYMBOL:Kdr" %in% nds)))

k <- grep("cancer", kegg_pathways_with_Kdr, invert = T)
kegg_pathways_with_Kdr <- kegg_pathways_with_Kdr[k]

reactome_pathways_with_Kdr <- reactome_pathways_with_Kdr[!reactome_pathways_with_Kdr %in% c("Developmental Biology", "Signaling Pathways", "Axon guidance")]

reactome_Kdr_nodes <- lapply(reactomeSmall[reactome_pathways_with_Kdr], nodes)
lapply(reactome_Kdr_nodes, function(nds) {"SYMBOL:Flt1" %in% nds})

# We need to get the sequence of equally spaced time points for each cell population.

daysCell <- tapply(pheno_mat$day, pheno_mat$cell_cond, unique)
equiTimeCell <- lapply(daysCell, function(x) {
  1:8
})
equiTimeCell$`ec-ko` <- 1:7
equiTimeCell$`mp-ko` <- 3:6
equiTimeCell$`per-wt` <- 1:4

# Now we can create the chunk to perform the analysis over all the cell population.
# Let's cicle over all cells.
# We start from the whole pathway analysis.

if (file.exists("activity_allCellTimeClipPathway.RData")) {
  load("activity_allCellTimeClipPathway.RData")
} else {
  cell_cond <- unique(as.character(pheno_mat$cell_cond))
  
  allCellTimeClipPathway <- lapply(cell_cond, function(cell) {
    cellExp <- activity[,as.character(pheno_mat$cell_cond) %in% cell, drop=F]
    row.names(cellExp) <- paste0("SYMBOL:",row.names(cellExp))
    times <- as.numeric(pheno_mat$day[as.character(pheno_mat$cell_cond) %in% cell])
  
    pathwayTimeKegg <- lapply(keggSmall, function(g) {
      g <- graphite::pathwayGraph(g)
      set.seed(1234)
      tryCatch(pathwayTimeCourse(cellExp, g, times, pcsMethod = "topological", 
                                 eqids = equiTimeCell[[cell]], maxPCs = 3),
               error = function(e) list(bestPc=e, alphas=NULL))
      })
    
    pathwayTimeReactome <- lapply(reactomeSmall, function(g) {
      g <- graphite::pathwayGraph(g)
      set.seed(1234)
      tryCatch(pathwayTimeCourse(cellExp, g, times, pcsMethod = "topological",
                                 eqids = equiTimeCell[[cell]], maxPCs = 3),
               error = function(e) list(bestPc=e, alphas=NULL))
      })
  
    list(kegg=pathwayTimeKegg, reactome=pathwayTimeReactome)
  })
  names(allCellTimeClipPathway) <- cell_cond
  save(allCellTimeClipPathway, file="activity_allCellTimeClipPathway.RData")
}
```
We run the clique analysis as well.
```{r}
if (file.exists("activity_allCellTimeClipModules.RData")) {
  load("activity_allCellTimeClipModules.RData")
} else {
  allCellTimeClipModules <- lapply(cell_cond, function(cell) {
    cellExp <- activity[,as.character(pheno_mat$cell_cond) %in% cell, drop=F]
    row.names(cellExp) <- paste0("SYMBOL:",row.names(cellExp))
    times <- as.numeric(pheno_mat$day[as.character(pheno_mat$cell_cond) %in% cell])
  
    cliqueTimeKegg <- lapply(keggSmall, function(g) {
      g <- graphite::pathwayGraph(g)
      set.seed(1234)
      tryCatch(cliqueTimeCourseTest(cellExp, g, times, pcsMethod = "sparse",
                                    eqids = equiTimeCell[[cell]],
                                    maxPCs = 3), error = function(e) list(bestPc=e, alphas=NULL))
    })
  
    cliqueTimeReactome <- lapply(reactomeSmall, function(g) {
      g <- graphite::pathwayGraph(g)
      set.seed(1234)
      tryCatch(cliqueTimeCourseTest(cellExp, g, times, pcsMethod = "sparse",
                                    eqids = equiTimeCell[[cell]], maxPCs = 3), error = function(e) list(bestPc=e, alphas=NULL))
    })
  
    list(kegg=cliqueTimeKegg, reactome=cliqueTimeReactome)
  })
  names(allCellTimeClipModules) <- cell_cond
  save(allCellTimeClipModules, file="activity_allCellTimeClipModules.RData")
}

# We would like to compare the pathways containing Kdr with at least one module significant.

kegg_sig_pathways_list <- list()
reactome_sig_pathways_list <- list()

source("accessory_function.R")
cell_cond = "ec-wt"
kegg_sig_pathways_list[[cell_cond]] <- check_significance_of_pathways(allCellTimeClipModules, pathways = kegg_pathways_with_Kdr, cell_cond = cell_cond, pathDB = "kegg", watchCliques = TRUE)

reactome_sig_pathways_list[[cell_cond]] <- check_significance_of_pathways(allCellTimeClipModules, pathways = reactome_pathways_with_Kdr, cell_cond = cell_cond, pathDB = "reactome", watchCliques = TRUE)

cell_cond = "fap-wt"
kegg_sig_pathways_list[[cell_cond]] <- check_significance_of_pathways(allCellTimeClipModules, pathways = kegg_pathways_with_Kdr, cell_cond = cell_cond, pathDB = "kegg", watchCliques = TRUE)

reactome_sig_pathways_list[[cell_cond]] <- check_significance_of_pathways(allCellTimeClipModules, pathways = reactome_pathways_with_Kdr, cell_cond = cell_cond, pathDB = "reactome", watchCliques = TRUE)

cell_cond = "mp-wt"
kegg_sig_pathways_list[[cell_cond]] <- check_significance_of_pathways(allCellTimeClipModules, pathways = kegg_pathways_with_Kdr, cell_cond = cell_cond, pathDB = "kegg", watchCliques = TRUE)

reactome_sig_pathways_list[[cell_cond]] <- check_significance_of_pathways(allCellTimeClipModules, pathways = reactome_pathways_with_Kdr, cell_cond = cell_cond, pathDB = "reactome", watchCliques = TRUE)

cell_cond = "inf-wt"
kegg_sig_pathways_list[[cell_cond]] <- check_significance_of_pathways(allCellTimeClipModules, pathways = kegg_pathways_with_Kdr, cell_cond = cell_cond, pathDB = "kegg", watchCliques = TRUE)

reactome_sig_pathways_list[[cell_cond]] <- check_significance_of_pathways(allCellTimeClipModules, pathways = reactome_pathways_with_Kdr, cell_cond = cell_cond, pathDB = "reactome", watchCliques = TRUE)


cell_cond = "per-wt"
kegg_sig_pathways_list[[cell_cond]] <- check_significance_of_pathways(allCellTimeClipModules, pathways = kegg_pathways_with_Kdr, cell_cond = cell_cond, pathDB = "kegg", watchCliques = TRUE)

reactome_sig_pathways_list[[cell_cond]] <- check_significance_of_pathways(allCellTimeClipModules, pathways = reactome_pathways_with_Kdr, cell_cond = cell_cond, pathDB = "reactome", watchCliques = TRUE)

## KO style

cell_cond = "ec-ko"
kegg_sig_pathways_list[[cell_cond]] <- check_significance_of_pathways(allCellTimeClipModules, pathways = kegg_pathways_with_Kdr, cell_cond = cell_cond, pathDB = "kegg", watchCliques = TRUE)

reactome_sig_pathways_list[[cell_cond]] <- check_significance_of_pathways(allCellTimeClipModules, pathways = reactome_pathways_with_Kdr, cell_cond = cell_cond, pathDB = "reactome", watchCliques = TRUE)


cell_cond = "fap-ko"
kegg_sig_pathways_list[[cell_cond]] <- check_significance_of_pathways(allCellTimeClipModules, pathways = kegg_pathways_with_Kdr, cell_cond = cell_cond, pathDB = "kegg", watchCliques = TRUE)

reactome_sig_pathways_list[[cell_cond]] <- check_significance_of_pathways(allCellTimeClipModules, pathways = reactome_pathways_with_Kdr, cell_cond = cell_cond, pathDB = "reactome", watchCliques = TRUE)


cell_cond = "mp-ko"
kegg_sig_pathways_list[[cell_cond]] <- check_significance_of_pathways(allCellTimeClipModules, pathways = kegg_pathways_with_Kdr, cell_cond = cell_cond, pathDB = "kegg", watchCliques = TRUE)

reactome_sig_pathways_list[[cell_cond]] <- check_significance_of_pathways(allCellTimeClipModules, pathways = reactome_pathways_with_Kdr, cell_cond = cell_cond, pathDB = "reactome", watchCliques = TRUE)

# Now we can collapse and create a dot plot.

library(ggplot2)
wts <- c("ec-wt", "fap-wt", "mp-wt", "inf-wt", "per-wt")
df <- cellCB:::meltListOfVectors(kegg_sig_pathways_list[wts])

ggplot(df) +
  geom_point(aes(x = nm, y = content))


library(ggplot2)

library(RColorBrewer)
colors <- brewer.pal(5, "Set1")
names(colors) <- c("ec-wt","fap-wt","inf-wt","mp-wt","per-wt")
colors.ko <- brewer.pal(5, "Set3")[c(4,5,3)]
names(colors.ko) <- c("ec-ko","fap-ko","mp-ko")
cls <- c(colors, colors.ko)

wts <- c("ec-wt", "fap-wt", "mp-wt", "inf-wt", "per-wt")
df <- cellCB:::meltListOfVectors(reactome_sig_pathways_list[wts])

extrafont::loadfonts()
ggplot(df) +
  geom_point(aes(x = nm, y = content, color=nm), size=5) +
  theme_bw() +
  theme(axis.title=element_blank()) +
  scale_color_manual(values = colors) + 
  theme(legend.position = "none") +
  theme(text = element_text(size=8, family = "Arial"),
        axis.text = element_text(size=8, family = "Arial"),
        plot.title = element_text(size=8))

ggsave(file="reactome-activation-timeClip.pdf")
# ggsave(file="fig-2d-reactome-activation-timeClip.pdf", height=4, width = 5.5)
ggsave(file="fig-2c-reactome-activation-timeClip.pdf", height=4, width = 5.5)


library(ggplot2)
kos <- c("ec-ko", "fap-ko", "mp-ko")
df <- cellCB:::meltListOfVectors(kegg_sig_pathways_list[kos])

ggplot(df) +
  geom_point(aes(x = nm, y = content))

library(ggplot2)
kos <- c("ec-ko", "fap-ko", "mp-ko")
df <- cellCB:::meltListOfVectors(reactome_sig_pathways_list[kos])

ggplot(df) +
  geom_point(aes(x = nm, y = content))


# Cdh5 profile.

ec_wt <- allCellTimeClipModules$`ec-wt`$reactome[["Signaling by VEGF"]]$cliquesExpr[[8]]
# allCellTimeClipModules$`ec-wt`$reactome[["Signaling by VEGF"]]$cliquesLoadings[[8]]

ec_ko <- allCellTimeClipModules$`ec-ko`$reactome[["Signaling by VEGF"]]$cliquesExpr[[8]]
ec <- cbind(ec_wt, ec_ko)
row.names(ec) <- gsub(pattern = "SYMBOL:", "", row.names(ec))
# pheatmap::pheatmap(ec, cluster_cols = F, scale = "row", gaps_col = 10)
# pheatmap::pheatmap(ec, cluster_cols = F, gaps_col = 10)
pheatmap::pheatmap(ec[c("Cdh5", "Ctnnb1"),], cluster_cols = F, scale = "row", gaps_col = 10)


ec_wt_cli8 <- allCellTimeClipModules$`ec-wt`$reactome[["Signaling by VEGF"]]$cliquesExpr[[8]]
# allCellTimeClipModules$`ec-wt`$reactome[["Signaling by VEGF"]]$cliquesLoadings[[8]]

ec_wt_kdr <- allCellTimeClipModules$`ec-wt`$reactome[["Signaling by VEGF"]]$cliquesExpr[[1]]["SYMBOL:Kdr",,drop=F]

ec_ko_cli8 <- allCellTimeClipModules$`ec-ko`$reactome[["Signaling by VEGF"]]$cliquesExpr[[8]]

ec_ko_kdr <- allCellTimeClipModules$`ec-ko`$reactome[["Signaling by VEGF"]]$cliquesExpr[[1]]["SYMBOL:Kdr",,drop=F]

ec_cli8 <- cbind(ec_wt_cli8, ec_ko_cli8)
row.names(ec_cli8) <- gsub(pattern = "SYMBOL:", "", row.names(ec_cli8))
ec_kdr <- cbind(ec_wt_kdr, ec_ko_kdr)
row.names(ec_kdr) <- gsub(pattern = "SYMBOL:", "", row.names(ec_kdr))

pheno_wt <- do.call(rbind, strsplit(colnames(ec_wt_kdr), "_"))
pheno_ko <- do.call(rbind, strsplit(colnames(ec_ko_kdr), "_"))

pheno <- rbind(pheno_wt, pheno_ko)

df <- data.frame(day=pheno[,1], cell=pheno[,2], kdr=t(ec_kdr), t(ec_cli8[c("Cdh5", "Ctnnb1"), , drop=F]), stringsAsFactors = F)
df$day <- factor(df$day, levels=unique(df$day))
df_plot <- tidyr::gather(df, "gene", "value", -c("day","cell"))

# We create a numeric version of the days

df$day_num <- as.numeric(gsub("d", "", df$day))
df_plot <- tidyr::gather(df, "gene", "value", -c("day","cell", "day_num"))
```

```{r}
library(ggplot2)
gene <- "Kdr"
ggplot(df_plot[df_plot$gene==gene, ], aes(x=day, color=cell, y=value, group=cell)) +
  geom_point() +
  geom_path() +
  ylab("activity score") +
  scale_color_manual(values=c(colors.ko["ec-ko"], colors["ec-wt"]))+
  theme_classic() +
  ggtitle(gene)
```

```{r}
library(ggplot2)
gene <- "Kdr"

ggplot(df_plot[df_plot$gene==gene, ], aes(x=day_num, color=cell, y=value, group=cell)) +
  geom_point() +
  stat_smooth(se=F) +
  ylab("activity score") +
  xlab("days") +
  scale_color_manual(values=c(colors.ko["ec-ko"], colors["ec-wt"])) +
  scale_x_continuous(breaks = unique(df_plot$day_num))+
  theme_classic() +
  ggtitle(gene)

ggsave(file=paste0(gene,"-ec-activity_profile.pdf"))

```


```{r}
library(ggplot2)
gene="Cdh5"
ggplot(df_plot[df_plot$gene==gene, ], aes(x=day, color=cell, y=value, group=cell)) +
  geom_point() +
  geom_path() +
  ylab("activity score") +
  scale_color_manual(values=c(colors.ko["ec-ko"], colors["ec-wt"]))+
  theme_classic() +
  ggtitle(gene)
```

```{r}
library(ggplot2)
gene <- "Cdh5"
ggplot(df_plot[df_plot$gene==gene, ], aes(x=day_num, color=cell, y=value, group=cell)) +
  geom_point() +
  stat_smooth(se=F) +
  ylab("activity score") +
  xlab("days") +
  scale_color_manual(values=c(colors.ko["ec-ko"], colors["ec-wt"]))+
  scale_x_continuous(breaks = unique(df_plot$day_num))+
  theme_classic() +
  ggtitle(gene)

ggsave(file=paste0(gene,"-ec-activity_profile.pdf"))

```

```{r}
library(ggplot2)
gene="Ctnnb1"
ggplot(df_plot[df_plot$gene==gene, ], aes(x=day, color=cell, y=value, group=cell)) +
  geom_point() +
  geom_path() +
  ylab("activity score") +
  scale_color_manual(values=c(colors.ko["ec-ko"], colors["ec-wt"]))+
  theme_classic() +
  ggtitle(gene)

# ec <- cbind(ec_wt, ec_ko)
# row.names(ec) <- gsub(pattern = "SYMBOL:", "", row.names(ec))
# pheatmap::pheatmap(ec, cluster_cols = F, scale = "row", gaps_col = 10)
# pheatmap::pheatmap(ec, cluster_cols = F, gaps_col = 10)
# pheatmap::pheatmap(ec[c("Cdh5", "Ctnnb1"),], cluster_cols = F, scale = "row", gaps_col = 10)

```

```{r}
library(ggplot2)
gene <- "Ctnnb1"
ggplot(df_plot[df_plot$gene==gene, ], aes(x=day_num, color=cell, y=value, group=cell)) +
  geom_point() +
  stat_smooth(se=F) +
  ylab("activity score") +
  scale_color_manual(values=c(colors.ko["ec-ko"], colors["ec-wt"]))+
  scale_x_continuous(breaks = unique(df_plot$day_num))+
  theme_classic() +
  ggtitle(gene)

ggsave(file=paste0(gene,"-ec-activity_profile.pdf"))

```

```{r}
which(allCellTimeClipModules$`mp-wt`$reactome[["Signaling by VEGF"]]$alphas<=0.05)

genes_ld_cl1 <- names(which(abs(allCellTimeClipModules$`mp-wt`$reactome[["Signaling by VEGF"]]$cliquesLoadings[[1]][,"PC2"]) >0.1))

genes_ld_cl4 <- names(which(abs(allCellTimeClipModules$`mp-wt`$reactome[["Signaling by VEGF"]]$cliquesLoadings[[4]][,"PC2"]) >0.1))

mp_wt <- rbind(allCellTimeClipModules$`mp-wt`$reactome[["Signaling by VEGF"]]$cliquesExpr[[4]][genes_ld_cl4,],
               allCellTimeClipModules$`mp-wt`$reactome[["Signaling by VEGF"]]$cliquesExpr[[1]][genes_ld_cl1,])
row.names(mp_wt) <- gsub("SYMBOL:", "", row.names(mp_wt))

mp_ko <- rbind(allCellTimeClipModules$`mp-ko`$reactome[["Signaling by VEGF"]]$cliquesExpr[[4]][genes_ld_cl4,],
               allCellTimeClipModules$`mp-ko`$reactome[["Signaling by VEGF"]]$cliquesExpr[[1]][genes_ld_cl1,])
row.names(mp_ko) <- gsub("SYMBOL:", "", row.names(mp_ko))

mp <- cbind(mp_wt, mp_ko)

pheatmap::pheatmap(mp, cluster_cols = F, gaps_col = 9)

pheno <- do.call(rbind, strsplit(colnames(mp), "_"))
df <- data.frame(day=pheno[,1], cell=pheno[,2], t(mp), stringsAsFactors = F)

df$day <- factor(df$day, levels=unique(df$day))
df_plot <- tidyr::gather(df, "gene", "value", -c("day","cell"))

# We create a numeric version of the days

df$day_num <- as.numeric(gsub("d", "", df$day))
df_plot <- tidyr::gather(df, "gene", "value", -c("day","cell", "day_num"))

library(ggplot2)
gene="Kdr"
ggplot(df_plot[df_plot$gene==gene, ], aes(x=day, color=cell, y=value, group=cell)) +
  geom_point() +
  geom_path() +
  ylab("activity score") +
  scale_color_manual(values=c(colors.ko["mp-ko"], colors["mp-wt"]))+
  theme_classic()+
  ggtitle(gene)

library(ggplot2)
gene="Kdr"
ggplot(df_plot[df_plot$gene==gene, ], aes(x=day_num, color=cell, y=value, group=cell)) +
  geom_point() +
  stat_smooth(se=F) +
  ylab("activity score") +
  xlab("days") + 
  scale_color_manual(values=c(colors.ko["mp-ko"], colors["mp-wt"]))+
  scale_x_continuous(breaks = unique(df_plot$day_num))+
  theme_classic()+
  ggtitle(gene)

ggsave(file=paste0(gene,"-mp-activity_profile.pdf"))

row.names(mp)
# "Itgav" "Itgb3" "Flt1"  "Kdr"   "Ptk2"  "Abi1"  "Pak1"
library(ggplot2)
gene="Pak1"
ggplot(df_plot[df_plot$gene==gene, ], aes(x=day, color=cell, y=value, group=cell)) +
  geom_point() +
  geom_path() +
  ylab("activity score") +
  scale_color_manual(values=c(colors.ko["mp-ko"], colors["mp-wt"]))+
  theme_classic()+
  ggtitle(gene)


gene="Pak1"
ggplot(df_plot[df_plot$gene==gene, ], aes(x=day_num, color=cell, y=value, group=cell)) +
  geom_point() +
  stat_smooth(se=F) +
  ylab("activity score") +
  xlab("days") + 
  scale_color_manual(values=c(colors.ko["mp-ko"], colors["mp-wt"]))+
  scale_x_continuous(breaks = unique(df_plot$day_num))+
  theme_classic()+
  ggtitle(gene)
ggsave(file=paste0(gene,"-mp-activity_profile.pdf"))

gene="Abi1"
ggplot(df_plot[df_plot$gene==gene, ], aes(x=day, color=cell, y=value, group=cell)) +
  geom_point() +
  geom_path() +
  ylab("activity score") +
  scale_color_manual(values=c(colors.ko["mp-ko"], colors["mp-wt"]))+
  theme_classic()+
  ggtitle(gene)


gene="Abi1"
ggplot(df_plot[df_plot$gene==gene, ], aes(x=day_num, color=cell, y=value, group=cell)) +
  geom_point() +
  stat_smooth(se=F) +
  ylab("activity score") +
  xlab("days") + 
  scale_color_manual(values=c(colors.ko["mp-ko"], colors["mp-wt"]))+
  scale_x_continuous(breaks = unique(df_plot$day_num))+
  theme_classic()+
  ggtitle(gene)

ggsave(file=paste0(gene,"-mp-activity_profile.pdf"))

gene="Itgav"
ggplot(df_plot[df_plot$gene==gene, ], aes(x=day, color=cell, y=value, group=cell)) +
  geom_point() +
  geom_path() +
  ylab("activity score") +
  scale_color_manual(values=c(colors.ko["mp-ko"], colors["mp-wt"]))+
  theme_classic()+
  ggtitle(gene)

gene="Itgav"
ggplot(df_plot[df_plot$gene==gene, ], aes(x=day_num, color=cell, y=value, group=cell)) +
  geom_point() +
  stat_smooth(se=F) +
  ylab("activity score") +
  xlab("days") + 
  scale_color_manual(values=c(colors.ko["mp-ko"], colors["mp-wt"]))+
  scale_x_continuous(breaks = unique(df_plot$day_num))+
  theme_classic()+
  ggtitle(gene)

ggsave(file=paste0(gene,"-mp-activity_profile.pdf"))
