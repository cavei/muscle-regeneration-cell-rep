compute_pcs_only_no_test <- function(g, expr, method=c("regular", "sparse"), shrink=FALSE, cliques=NULL, maxPCs=5) {
    cliExp <- expr[g, ,drop=FALSE]
    
    if (NROW(expr) == 0) {
      return(NULL)
    }
    cliExp <- t(cliExp) ## check this
    
    if (NCOL(cliExp)!=1) {
      pcs <- houseOfClipUtility::computePCs(cliExp, shrink=shrink, method=method, cliques=cliques, maxPCs=maxPCs)
    } else {
      colnames(cliExp) <- "PC1"
      pcs <- list(x=cliExp, sdev=sd(cliExp), loadings=1)
    }
    
    list(bestPc=NA, cliAlphas = NA, pcs=pcs$x, cliqueLoadings=pcs$loadings, cliqueExp=cliExp)
}
