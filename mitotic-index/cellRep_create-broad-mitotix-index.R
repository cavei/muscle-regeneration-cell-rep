# Mitotic index collapsing times

library(pheatmap)
library(ggplot2)
library(ggpubr)

dataDirname <- "../preprocessing-data/rpkm-RDatas-padj0.01-lfc2/"
resultsDir="../preprocessing-data/rpkm-RDatas-padj0.01-lfc2/"

# Now let's load all the data.
# The computation of the Mitotic Index is based on the following papers:
# https://academic.oup.com/nar/article/46/14/7022/5035173
# and its original paper
# "https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1064-3#Sec11".
# 
# __Caveats:__ RuvSeq Data are collected without filter. RuvSeq data are collected with filters.
# ```{r load the data, echo=FALSE}

library(cellCB)

expressionMatCell <- lapply(c("fap", "ec", "mp", "inf", "per"), function(cell) {
  getNormExpressionFromRData(paste0(resultsDir, "/", cell, "-wt.RData"), paste0(cell, "_wt"))
})
names(expressionMatCell) <- c("fap", "ec", "mp", "inf", "per")
averageRuvExprCell <- lapply(expressionMatCell, averageReplicate)

expressionRpkmCell <- lapply(c("fap", "ec", "mp", "inf", "per"), function(cell) {
  getRpkmExpressionFromRData(paste0(resultsDir, "/", cell, "-wt.RData"), paste0(cell, "_wt"))
})

names(expressionRpkmCell) <- c("fap", "ec", "mp", "inf", "per")
averageRpkmExprCell <- lapply(expressionRpkmCell, averageReplicate)

# We load the set of genes for the mitotic index. I manually translated human Symbols in Mouse orthologs.

mitoticList <- read.table("mRNA-for-mitotic-risk.txt", sep="\t", header=T, quote="\"", stringsAsFactors = F, check.names = F)

# Now we are ready to calculate.

source("cellRep_plot_index.R")
wideMitoticActivityRuv <- lapply(expressionRpkmCell, computeMitoticIndex, migenes=mitoticList$mouse)

ko.expressionMatCell <- lapply(c("fap", "ec", "mp"), function(cell) {
  getNormExpressionFromRData(paste0(resultsDir, "/", cell, "-ko.RData"), paste0(cell, "_ko"))
})
names(ko.expressionMatCell) <- c("fap", "ec", "mp")
ko.averageRuvExprCell <- lapply(ko.expressionMatCell, averageReplicate)

ko.expressionRpkmCell <- lapply(c("fap", "ec", "mp"), function(cell) {
  getRpkmExpressionFromRData(paste0(resultsDir, "/", cell, "-ko.RData"), paste0(cell, "_ko"))
})
names(ko.expressionRpkmCell) <- c("fap", "ec", "mp")

ko.averageRpkmExprCell <- lapply(ko.expressionRpkmCell, averageReplicate)

wide.ko.mitoticActivityRuv <- lapply(ko.expressionRpkmCell, computeMitoticIndex, migenes=mitoticList$mouse)

# It's time for plots. WE cicle and create all the Mitotic Indexes.

library(RColorBrewer)
colors <- brewer.pal(5, "Set1")
names(colors) <- c("ec","fap","inf","mp","per")
colors.ko <- brewer.pal(5, "Set3")[c(4,5,3)]
names(colors.ko) <- c("ec","fap","mp")

breaks <- c(0,1,4,10,15)-1
names(breaks) <- c("steady", "early", "middle", "late","inf")

dfs <- lapply(names(wide.ko.mitoticActivityRuv), create_mitotic_index_ggpubr_data, wt=wideMitoticActivityRuv, ko=wide.ko.mitoticActivityRuv, force_cats=breaks)
names(dfs) <- names(wide.ko.mitoticActivityRuv)
source("ggpubr_barplot.R")

barplot_ggpubr(dfs$fap, colors=unname(c(colors["fap"], colors.ko["fap"])), l_height = 56)
barplot_ggpubr(dfs$ec, colors=unname(c(colors["ec"], colors.ko["ec"])), l_height = 50)
barplot_ggpubr(dfs$mp, colors=unname(c(colors["mp"], colors.ko["mp"])), l_height = 90)

lonely_dfs <- lapply(c("inf", "per"), create_mitotic_index_ggpubr_data, wt=wideMitoticActivityRuv, ko=wide.ko.mitoticActivityRuv, removeCondition="ko", force_cats=breaks)
names(lonely_dfs) <- c("inf", "per")

# Barplot Paper

p <- barplot_ggpubr(dfs$fap, colors=c(rgb(15,0,255, max=255)	, rgb(30, 239, 228, max=255)), l_height = 56)
ggsave(p, filename = "fap-mitotec.index.activity-broad-plot.pdf", width=5, height=3.2)

p <- barplot_ggpubr(dfs$ec, colors=c(rgb(251, 0, 5	, max=255)	, rgb(236, 139, 126	, max=255)), l_height = 50)
ggsave(p, filename = "ec-mitotec.index.activity-broad-plot.pdf", width=5, height=3.2)

p <- barplot_ggpubr(dfs$mp, colors=c(rgb(87, 0, 255	, max=255)	, rgb(196, 153, 229	, max=255)), l_height = 90)
ggsave(p, filename = "mp-mitotec.index.activity-broad-plot.pdf", width=5, height=3.2)

