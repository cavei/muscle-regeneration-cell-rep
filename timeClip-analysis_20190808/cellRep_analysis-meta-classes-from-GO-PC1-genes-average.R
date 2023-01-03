# Meta class analysis with PC1 gene average

library(timeClip)

dataDirname <- "../preprocessing-data/rpkm-RDatas-padj0.01-lfc2/"
resultsDir <- "../preprocessing-data/rpkm-RDatas-padj0.01-lfc2/"

# Now let's load all the data. We create the normPsedo Counts for all the cel population. After that we average the same time points to have one replicate per time point.

library(cellCB)
expressionTot <- getNormExpressionFromRData(paste0(resultsDir, "/", "tot", "-wt.RData"), paste0("tot", "_wt"))
averageRuvExprTot <- averageReplicate(expressionTot)

expressionMatCell <- lapply(c("fap", "ec", "mp", "inf", "per"), function(cell) {
  getNormExpressionFromRData(paste0(resultsDir, "/", cell, "-wt.RData"), paste0(cell, "_wt"))
})
names(expressionMatCell) <- c("fap", "ec", "mp", "inf", "per")
averageRuvExprCell <- lapply(expressionMatCell, averageReplicate)

expressionMatCell.ko <- lapply(c("fap", "ec", "mp"), function(cell) {
  getNormExpressionFromRData(paste0(resultsDir, "/", cell, "-ko.RData"), paste0(cell, "_ko"))
})
names(expressionMatCell.ko) <- c("fap", "ec", "mp")
averageRuvExprCell.ko <- lapply(expressionMatCell.ko, averageReplicate)

# Two list with cells and their psudoCounts, plus the total muscle.
# Now we need the set of genes from GO analysis. We load an object that will release two objects selected_8GO_to_genes_ko and selected_8GO_to_genes.

# load dataset from GO
load("../go_analysis_20190808/selected_8GO_to_genes.RData")

# We create a list with equally spaced time points for each cell population.

timeTot <- as.numeric(sub("d", "",colnames(averageRuvExprTot)))
equiTimeTot <- 1:8
  
timesCell <- lapply(averageRuvExprCell, function(cellExp) {
  times <- as.numeric(sub("d", "",colnames(cellExp)))
})
equiTimeCell <- lapply(averageRuvExprCell, function(cellExp) {
  1:8
})
equiTimeCell$per <- 1:4

timesCell.ko <- lapply(averageRuvExprCell.ko, function(cellExp) {
  times <- as.numeric(sub("d", "",colnames(cellExp)))
})

equiTimeCell.ko <- lapply(averageRuvExprCell.ko, function(cellExp) {
  1:8
})
equiTimeCell.ko$ec <- 1:7
equiTimeCell.ko$mp <- 3:6


# Now we can create the chunk to perform the analysis over all the cell populations.
# Let's cicle over all cells.
# We are going to produce a timeClip result to extract PC information.

if (file.exists("totTimeClipGO_Sets.RData")) {
  load("totTimeClipGO_Sets.RData")
} else {
  normPseudoCount <- log2(averageRuvExprTot+1)
  times <- as.numeric(sub("d", "",colnames(normPseudoCount)))
  
  totTimeClipGO_Sets <- lapply(selected_8GO_to_genes[["tot"]], function(g) {
    set.seed(1234)
    tryCatch(singleCliqueTest(g, normPseudoCount, times, method = "sparse",
                               eqids = equiTimeTot, maxPCs = 8),
             error = function(e) list(bestPc=e, alphas=NULL))
  })
  save(totTimeClipGO_Sets, file="totTimeClipGO_Sets.RData")
}

source("cellRep_compute_pcs_only.R")
if (file.exists("allCellTimeClipGO_Sets.RData")) {
  load("allCellTimeClipGO_Sets.RData")
} else {
  allCellTimeClipGO_sets <- lapply(names(averageRuvExprCell), function(cell) {
    
    cellExp <- averageRuvExprCell[[cell]]
    normPseudoCount <- log2(cellExp+1)
    # normPseudoCount <- t(scale(t(normPseudoCount), center = T, scale = T))
    times <- as.numeric(sub("d", "",colnames(normPseudoCount)))
    
    go_gene_set_timeClip <- lapply(selected_8GO_to_genes[[cell]], function(g) {
      set.seed(1234)
      tryCatch(singleCliqueTest(g, normPseudoCount, times, method = "sparse",
                                 eqids = equiTimeCell[[cell]], maxPCs = 8),
               error = function(e) compute_pcs_only_no_test(g, normPseudoCount, method = "sparse", maxPCs = 8))
      })
  })
  names(allCellTimeClipGO_sets) <- names(averageRuvExprCell)
  save(allCellTimeClipGO_sets, file="allCellTimeClipGO_Sets.RData")
}

if (file.exists("allCellKoTimeClipGO_Sets.RData")) {
  load("allCellKoTimeClipGO_Sets.RData")
} else {
  allCellKoTimeClipGO_sets <- lapply(names(averageRuvExprCell.ko), function(cell) {
    
    cellExp <- averageRuvExprCell.ko[[cell]]
    normPseudoCount <- log2(cellExp+1)
    # normPseudoCount <- t(scale(t(normPseudoCount), center = T, scale = T))
    times <- as.numeric(sub("d", "",colnames(normPseudoCount)))
    
    go_gene_set_timeClip <- lapply(selected_8GO_to_genes[[cell]], function(g) {
      set.seed(1234)
      g <- intersect(g, row.names(normPseudoCount))
      if (length(g)==0)
        return(list(bestPc=NULL, alphas=NULL))
      tryCatch(singleCliqueTest(g, normPseudoCount, times, method = "sparse",
                                 eqids = equiTimeCell.ko[[cell]], maxPCs = 8),
               error = function(e) compute_pcs_only_no_test(g, normPseudoCount, method = "sparse", maxPCs = 8))
      })
  })
  names(allCellKoTimeClipGO_sets) <- names(averageRuvExprCell.ko)
  save(allCellKoTimeClipGO_sets, file="allCellKoTimeClipGO_sets.RData")
}

# We need to inspect the data for WT. First, we are going to plot the PC1 for each category.

cell_data <- totTimeClipGO_Sets
source("cellRep_extract_pc_profile_and_related_genes.R")

tot_pcs_and_their_genes <- extract_pc_profile_and_related_genes(cell_data)

cell_pcs_and_their_genes <- lapply(allCellTimeClipGO_sets, function(cell_data) {
  pcs_and_their_genes <- extract_pc_profile_and_related_genes(cell_data)
})

cell_ko_pcs_and_their_genes <- lapply(allCellKoTimeClipGO_sets, function(cell_data) {
  pcs_and_their_genes <- extract_pc_profile_and_related_genes(cell_data)
})

# Now we create a palette for the meta_classes.

cols <- RColorBrewer::brewer.pal(8, "Dark2")
cols <- cols[c(1:4, 8, 5, 6, 7)]

names(cols) <- names(tot_pcs_and_their_genes)

# Plot the PC1 lines for Total Muscle

source("cellRep_create_trend_lines.R")
save_to_dir = "cell_population_8category_PC1_average"
dir.create(save_to_dir, showWarnings = F)

p <- create_trend_lines_for_genes(tot_pcs_and_their_genes, averageRuvExprTot, title="Total Muscle",
                                  discrete_palette = cols,
                                  relative = F, scale=F)
ggsave(p, filename = paste0(save_to_dir, "/", "total_muscle_8_categories_trend.pdf"), height=7, width = 10)

### CELLS pop
# We can go on analyzyng the cell populations. First FAP wt.

cell="fap"
p <- create_trend_lines_for_genes(cell_pcs_and_their_genes[[cell]], averageRuvExprCell[[cell]],
                                  title="FAP wt", discrete_palette = cols,
                                  relative = F, scale=F)
ggsave(p, filename = paste0(save_to_dir, "/", cell, "_wt_8_categories_trend.pdf"),  height=7, width = 10)
p

cell="ec"
p <- create_trend_lines_for_genes(cell_pcs_and_their_genes[[cell]], averageRuvExprCell[[cell]],
                                  title="EC wt", discrete_palette = cols,
                                  relative = F, scale=F)
ggsave(p, filename = paste0(save_to_dir, "/", cell, "_wt_8_categories_trend.pdf"),  height=7, width = 10)

cell="mp"
p <- create_trend_lines_for_genes(cell_pcs_and_their_genes[[cell]], averageRuvExprCell[[cell]],
                                  title="MP wt", discrete_palette = cols,
                                  relative = F, scale=F)
ggsave(p, filename = paste0(save_to_dir, "/", cell, "_wt_8_categories_trend.pdf"),  height=7, width = 10)

cell="inf"
p <- create_trend_lines_for_genes(cell_pcs_and_their_genes[[cell]], averageRuvExprCell[[cell]],
                                  title="INF wt", discrete_palette = cols,
                                  relative = F, scale=F)
ggsave(p, filename = paste0(save_to_dir, "/", cell, "_wt_8_categories_trend.pdf"),  height=7, width = 10)

cell="per"
p <- create_trend_lines_for_genes(cell_pcs_and_their_genes[[cell]], averageRuvExprCell[[cell]],
                                  title="PER wt", discrete_palette = cols,
                                  relative = F, scale=F)
ggsave(p, filename = paste0(save_to_dir, "/", cell, "_wt_8_categories_trend.pdf"),  height=7, width = 10)

cell="fap"
p <- create_trend_lines_for_genes(cell_ko_pcs_and_their_genes[[cell]], averageRuvExprCell.ko[[cell]],
                                  title="FAP ko", discrete_palette = cols,
                                  relative = F, scale=F)
ggsave(p, filename = paste0(save_to_dir, "/", cell, "_ko_8_categories_trend.pdf"),  height=7, width = 10)
p

cell="ec"
p <- create_trend_lines_for_genes(cell_ko_pcs_and_their_genes[[cell]], averageRuvExprCell.ko[[cell]],
                                  title="EC ko", discrete_palette = cols,
                                  relative = F, scale=F)
ggsave(p, filename = paste0(save_to_dir, "/", cell, "_ko_8_categories_trend.pdf"),  height=7, width = 10)

cell="mp"
p <- create_trend_lines_for_genes(cell_ko_pcs_and_their_genes[[cell]], averageRuvExprCell.ko[[cell]],
                                  title="MP ko", discrete_palette = cols,
                                  relative = F, scale=F)
ggsave(p, filename = paste0(save_to_dir, "/", cell, "_ko_8_categories_trend.pdf"),  height=7, width = 10)

cell="fap"
plots_p <- lapply(names(cell_pcs_and_their_genes[[cell]]), function(meta_class) {
  create_trend_lines_comparison(pcs_and_their_genes=cell_pcs_and_their_genes[[cell]],
                                meta_class=meta_class, 
                                expr_data=averageRuvExprCell[[cell]],
                                expr_data_cmp=averageRuvExprCell.ko[[cell]],
                                palette="Set2", title = paste(cell, meta_class),
                                relative = F, scale=F)
})
library(gridExtra)
mp <- marrangeGrob(grob=plots_p, ncol=3, nrow=3)
ggsave(mp, filename = paste0(save_to_dir, "/", cell, "_wt_ko_8_categories_trend.pdf"),  height=7, width = 10)


cell="ec"
plots_p <- lapply(names(cell_pcs_and_their_genes[[cell]]), function(meta_class) {
  create_trend_lines_comparison(pcs_and_their_genes=cell_pcs_and_their_genes[[cell]],
                                meta_class=meta_class, 
                                expr_data=averageRuvExprCell[[cell]],
                                expr_data_cmp=averageRuvExprCell.ko[[cell]],
                                palette="Set2", title = paste(cell, meta_class),
                                relative = F, scale=F)
})
library(gridExtra)
mp <- marrangeGrob(grob=plots_p, ncol=3, nrow=3)
ggsave(mp, filename = paste0(save_to_dir, "/", cell, "_wt_ko_8_categories_trend.pdf"),  height=7, width = 10)


cell="mp"
plots_p <- lapply(names(cell_pcs_and_their_genes[[cell]]), function(meta_class) {
  create_trend_lines_comparison(pcs_and_their_genes=cell_pcs_and_their_genes[[cell]],
                                meta_class=meta_class, 
                                expr_data=averageRuvExprCell[[cell]],
                                expr_data_cmp=averageRuvExprCell.ko[[cell]],
                                palette="Set2", title = paste(cell, meta_class),
                                relative = F, scale=F)
})
library(gridExtra)
mp <- marrangeGrob(grob=plots_p, ncol=3, nrow=3)
ggsave(mp, filename = paste0(save_to_dir, "/", cell, "_wt_ko_8_categories_trend.pdf"),  height=7, width = 10)

