# Analyze module cell type specificity in mousebrain.org (Zeisel et al. 2018)
# data, at three different levels.
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)
library(preprocessCore)
library(pSI)
source("R/expression/DE/DE-functions.R")

# list files
files = list.files("data/expression/mousebrain.org", pattern = "sc*", 
                   full.names = T)
for (file in files) {
  # read in processed expression data
  data = read.delim(file)
  genes = data[, 1]
  data = data[, -1]
  
  # create matrix
  expr = as.matrix(data)
  
  # quantile normalize to average empirical distribution across all samples
  mat = normalize.quantiles(expr)
  
  # quantile normalize each gene to standard normal distribution
  quantNorm = function(x) qnorm(
    rank(x, ties.method = "average") / (length(x) + 1))
  mat.n = apply(mat, 2, quantNorm)
  
  # convert back to data frame to set dimension names
  exprEnsembl = as.data.frame(mat.n)
  exprEnsembl = cbind(genes, exprEnsembl)
  colnames(exprEnsembl) = c("Gene", colnames(expr))
  
  # load GTEx modules
  modules = read.delim("data/modules/GTEx-modules.tsv")
  
  # map to human consensus orthologs
  orthologs = read.delim("data/orthologs/biomart-ensembl-orthologs-mouse.tsv")
  orthologs = orthologs[orthologs[, 2] != "",] # remove unmapped orthologs
  # pick a single ortholog when multiple match
  orthologs = pick.orthologs(orthologs, modules$gene)
  exprEnsembl[, 1] = orthologs[match(exprEnsembl[, 1], orthologs[, 1]), 2]
  exprEnsembl = exprEnsembl[!is.na(exprEnsembl[, 1]),] # remove NA
  exprEnsembl = exprEnsembl[which(exprEnsembl[, 1] != ""),] # remove blank
  
  # aggregate again, if needed
  exprEnsembl = aggregate(exprEnsembl[, -1], by = list(exprEnsembl[, 1]),
                          FUN = median)
  # re-format as matrix
  geneNames = exprEnsembl[,1]
  exprEnsembl = exprEnsembl[,-1]
  rownames(exprEnsembl) = geneNames
  exprEnsembl = as.matrix(exprEnsembl)
  
  # get 25% percentile as cutoff
  cutoff = quantile(exprEnsembl, 0.25)
  
  # get pSI values from cell type specific data
  pSI.output = suppressWarnings(specificity.index(
    pSI.in = exprEnsembl, e_min = cutoff))
  
  # analyze each module in trn
  results = data.frame(id = character(0), cell = character(0), p = numeric(0))
  for (module in 1:15) {
    genes = modules$gene[modules$module == module]
    fisher = fisher.iteration(pSI.output, genes, background = "data.set", 
                              p.adjust = T)
    cells = rownames(fisher)
    p = fisher[2]
    out = data.frame(module = rep(module, nrow(fisher)), cell = cells, p = p)
    colnames(out) = c("module", "cell", "p")
    results = rbind(results, out)
  }
  # write output
  write.table(results, paste0("data/expression/mousebrain.org/pSI_", gsub(
    "sc_|\\.txt.*$", "", basename(file)), ".txt"), 
    quote = F, row.names = F, sep = "\t")
}
