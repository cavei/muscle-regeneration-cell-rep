
# Assign gene to cell

library(cellCB)

result_dir <- "../preprocessing-data/rpkm-RDatas-padj0.01-lfc2/"
rdata_storage <- "rdata_storage"
dir.create(rdata_storage, showWarnings = F)

# We compare the wt populations.

if (!(file.exists("gene2cellAssociation.RData"))) {
  cellExpressionTable <- getCellExpressionSignature(RData = paste0(result_dir, "/all-wt.RData"))
  gene2cellAssociationTables <- cellCB:::associateGeneToCells(cellExpressionTable)
  gene2cellAssociation <- cellCB:::tables2maps(gene2cellAssociationTables)
  file=paste0(rdata_storage, "/gene2cellAssociation-padj0.01-lfc2-", as.character(Sys.Date()), ".RData")
  link="gene2cellAssociation.RData"
  save(gene2cellAssociation, file=file)
  file.symlink(file, link)
}

# Now the KOs.

if (!(file.exists("gene2cellAssociation-ko.RData"))) {
  cellExpressionTable <- getCellExpressionSignature(RData = paste0(result_dir, "/all-ko.RData"), objname = "all_ko", cols = 1:2)
  gene2cellAssociationTables <- cellCB:::associateGeneToCells(cellExpressionTable, winCmpToBeExclusive = 2)
  gene2cellAssociation <- cellCB:::tables2maps(gene2cellAssociationTables)
  file=paste0(result_dir, "/gene2cellAssociation-padj0.01-lfc2-ko-", as.character(Sys.Date()), ".RData")
  link="gene2cellAssociation-ko.RData"
  save(gene2cellAssociation, file=file)
  file.symlink(file, link)
}
