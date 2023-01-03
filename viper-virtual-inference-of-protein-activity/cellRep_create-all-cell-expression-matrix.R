library(cellCB)
resultsDir <- paste0("../preprocessing-data/","/rpkm-Rdatas-padj0.01-lfc2/")

all_cells <- getNormExpressionFromRData(paste0(resultsDir, "/", "all-wt_ko.RData"),"all_wt_ko")

all_cell_expr <- data.frame(symbol=row.names(all_cells), all_cells, check.names = F)

write.table(all_cell_expr, file="all_cell_expr-norm-counts.txt", sep="\t", quote=F, row.names=F)

all_cell_log <- log2(all_cells+1)
all_cell_expr_log <- data.frame(symbol=row.names(all_cells), all_cell_log, check.names = F)
write.table(all_cell_expr_log, file="all_cell_expr-log-norm-counts.txt", sep="\t", quote=F, row.names=F)


collapsed_pheno <- gsub("-", "_", colnames(all_cell_log))
pheno <- data.frame(do.call(rbind, strsplit(collapsed_pheno, "_")),
                    row.names=colnames(all_cell_log), stringsAsFactors = F)
colnames(pheno) <- c("rep", "day", "cell", "condition")

library(Biobase)
eset <- ExpressionSet(assayData = as.matrix(all_cell_log), phenoData = AnnotatedDataFrame(pheno))

save(eset, file="all-cell-eset.RData")
