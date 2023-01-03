# Create the ligand receptors tables

# Ligand Receptor analysis

library(org.Mm.eg.db)
library(clusterProfiler)
library(ggplot2)
library(cellCB)

# We load the results from the Differentially expressed analyses.

resultsDir <- "../preprocessing-data/rpkm-RDatas-padj0.01-lfc2/"
tableDir = "sif_tables"
dir.create(tableDir, showWarnings = F)

gene2cell_association_dir = "../associate-gene-to-cell/"
load(paste0(gene2cell_association_dir, "gene2cellAssociation.RData"))

# We have created a list of genes expressed by each cell
# We need to know for each cell those gene that vary during the time


cells <- c("fap", "ec", "mp", "inf", "per")
activeAtDaysCell <- lapply(cells, function(cell) {
  cellCB:::getActiveAtDaysFromRData(paste0(resultsDir, "/", cell, "-wt.RData"), paste0(cell, "_wt"))
})
names(activeAtDaysCell) <- cells

activeAtDaysTot <- cellCB:::getActiveAtDaysFromRData(paste0(resultsDir, "/tot-wt.RData"), "tot_wt")

genesBehavioursCell <- lapply(cells, function(cell) {
  extractGenesBehavious(activeAtDaysCell[[cell]], activeAtDaysTot, gene2cellAssociation[[cell]])
})
names(genesBehavioursCell) <- cells

# Let our cell population be fap, we came up with a list of time points and genes that are activated in fap.
# A list of time point were we stored genes that are active in total muscle but that belong also to fap.
# And a list of genes that constitutively are expressed by fap.
# Now for our purpose we isolate the receptors and ligands. Let's source the database.

db <- read.table("../ligands-receptor-rezza/database-mouse-sif.txt",
                 header=F, sep="\t", quote="\"", stringsAsFactors = F)

rezza.db <- list(ligands=unique(db$V2), receptors=unique(db$V3))

# We create for each cell the sets of ligands and receptors active at different times.

ligandReceptorBehavioursCell <- lapply(cells, function(cell)
  extractLigandReceptorBehaviours(genesBehavioursCell[[cell]], rezza.db))
names(ligandReceptorBehavioursCell) <- cells

# We have collected some list of receptor/ligand activation
# Now we can create the timeSpecific network.
# We are going to create a graph to hanlde ligand receptor interaction.

library(igraph)
library(graph)
lr.sif <- as.matrix(db[,2:3])
colnames(lr.sif) <- c("src", "dest")
receptorGraph <- igraph::as_graphnel(igraph::graph(hyperGraphite::edgeList(lr.sif), directed=TRUE))

# Receptor graph contains the graph of the ligand interactions.
# Now we can extract the tables of the active interactions.

# We transform the data and we collect all the ligands and receptor with their timing

ligandsTimingCell <- lapply(names(ligandReceptorBehavioursCell), function(cell) {
    cdata <- ligandReceptorBehavioursCell[[cell]]
    batch <- meltListOfVectors(lapply(cdata$ligandReceptorDayBatch, function(ddata) ddata$ligandOfTheDay))
    names(batch)=c("day", "ligand")
    batch$source <- "batch"
    
    tot <- meltListOfVectors(lapply(cdata$ligandReceptorDayTot, function(ddata) ddata$ligandOfTheDay))
    names(tot)=c("day", "ligand")
    tot$source <- "total"
    
    const <- data.frame(day = "all", ligand= cdata$ligandReceptorConstitutive$ligandOfTheDay)
    const$source <- "constitutive"
    ltc <- rbind(batch, tot, const)
    ltc$cell <- cell
    ltc
  })

names(ligandsTimingCell) <- names(ligandReceptorBehavioursCell)

receptorsTimingCell <- lapply(names(ligandReceptorBehavioursCell), function(cell) {
    cdata <- ligandReceptorBehavioursCell[[cell]]
    batch <- meltListOfVectors(lapply(cdata$ligandReceptorDayBatch, function(ddata) ddata$receptorOfTheDay))
    names(batch)=c("day", "receptor")
    batch$source <- "batch"
    
    tot <- meltListOfVectors(lapply(cdata$ligandReceptorDayTot, function(ddata) ddata$receptorOfTheDay))
    names(tot)=c("day", "receptor")
    tot$day <- paste0("t_",tot$day)
    tot$source <- "total"
    
    const <- data.frame(day = "all", receptor= cdata$ligandReceptorConstitutive$receptorOfTheDay)
    const$source <- "constitutive"
    
    ltc <- rbind(batch, tot, const)
    ltc$cell <- cell
    ltc
  })
names(receptorsTimingCell) <- names(ligandReceptorBehavioursCell)

saveFiles = TRUE
if (saveFiles) {
  if (!file.exists(tableDir))
    dir.create(tableDir)
}

sifTables <- lapply(names(ligandReceptorBehavioursCell), function(cell) {
  removeAutocrine=FALSE
  useConstitutiveLigands=TRUE
  cat(cell, "\n")
  condition="wt"
  
  receptorFrom <- cell
  autocrine <- wiseMergeTimePoints(ligandReceptorBehavioursCell[[cell]])
  ligandsReceptorCumulative <- paraWiseMultiMergeTimePoints(ligandReceptorBehavioursCell, ligCells = cells,
                                                            recCell = receptorFrom,
                                                            removeAutocrine=removeAutocrine,
                                                            useConstitutiveLigands=useConstitutiveLigands)
  
  paracrine <- paraWiseMultiMergeTimePoints(ligandReceptorBehavioursCell, ligCells = cells,
                                                            recCell = receptorFrom,
                                                            removeAutocrine=removeAutocrine,
                                                            useConstitutiveLigands=useConstitutiveLigands)
  
  dailySif <- cellCB:::createSifOFLigandReceptor(ligandsReceptorCumulative, receptorGraph, cell, condition)
  # days intersection
  days <- intersect(names(autocrine), names(dailySif))
  for (d in days) {
    dailySif[[d]]$day <- d
    dailySif[[d]]$feedSystem <- "paracrine"
    
    isAutocrine <- dailySif[[d]]$src %in% autocrine[[d]]$ligandOfTheDay
    isParacrine <- dailySif[[d]]$src %in% paracrine[[d]]$ligandOfTheDay
    dailySif[[d]]$feedSystem[isAutocrine] <- "autocrine"
    dailySif[[d]]$feedSystem[isAutocrine & isParacrine] <- "autocrine;paracrine"
    
    dailySif[[d]]$cell <- cell
    # colnames(dailySif[[d]]) <- c("src","dest","day", "feedSystem", "cell")
    colnames(dailySif[[d]]) <- c("ligand","receptor","day", "feedSystem", "cell")
    dailySif[[d]] <- dailySif[[d]][c("day", "ligand", "receptor", "cell", "feedSystem")]
    if (saveFiles) {
      write.table(dailySif[[d]], 
                  file=paste0(tableDir, "/", paste(c(cell, condition, d), collapse = "-"), ".txt"),
                  sep="\t", quote=F, row.names=F) 
    }
  }
  
  # if (make.video)
  #   renderHTMLvideo2.0(ligandsReceptorCumulative, receptorGraph, cell, condition)
  
  dailySif
})
names(sifTables) <- names(ligandReceptorBehavioursCell)

# We have collected the tables with the interactome network.

ligandMatrix <- do.call(rbind, ligandsTimingCell)
row.names(ligandMatrix) <- NULL
write.table(ligandMatrix, file=paste0(tableDir, "/", "ligandDescription.txt"), sep="\t", quote=F, row.names=F)

receptorMatrix <- do.call(rbind, receptorsTimingCell)
row.names(receptorMatrix) <- NULL
write.table(receptorMatrix, file=paste0(tableDir, "/", "receptorDescription.txt"), sep="\t", quote=F, row.names=F)

save(sifTables, receptorMatrix, ligandMatrix, file=paste0(tableDir, '/sifTables.RData') )
