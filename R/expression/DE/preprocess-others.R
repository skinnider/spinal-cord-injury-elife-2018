# Preprocess generic array data. 
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)
library(affy)
library(limma)
library(genefilter)
library(mouse4302.db)
library(rat2302.db)
source("R/expression/DE/DE-functions.R")

dbs = list(mouse4302.db, rat2302.db, rat2302.db)
arrays = list("A-AFFY-45", "A-AFFY-43", "A-AFFY-43")
expts = list("E-GEOD-5296", "E-GEOD-45006", "E-GEOD-69334")

# load GTEx modules and seed proteins 
modules = read.delim("data/modules/GTEx-modules.tsv")
sci = read.delim("data/literature_review/Table_S1.txt")
seeds = unique(sci$ensembl)

for (i in seq_len(length(arrays))) {
  db = dbs[[i]]
  array = arrays[[i]]
  expt = expts[[i]]
  
  # Read experiment targets
  expt.dir = file.path("data/expression/arrayexpress/cel", array, expt)
  targets = readTargets(file.path("data/expression/DE", expt, 
                                  "targets.tsv"))
  
  # read CEL files 
  samples = file.path(expt.dir, targets$filename)
  abatch = ReadAffy(filenames=samples)
  
  # normalize using MAS5
  eset = mas5(abatch)
  calls = mas5calls(abatch)
  save(eset, file = paste("data/expression/DE/", expt, "/", expt, 
                          "-mas5.Rdata", sep=""))
  save(calls, file = paste("data/expression/DE/", expt, "/", expt, 
                           "-calls.Rdata", sep=""))
  expr = exprs(eset)
  calls = exprs(calls)
  
  # filter genes expressed in < 20% of samples 
  samples = colnames(expr)
  f1 = function(x) (length(which(x == "P")) > 0.8 * length(samples))
  selected = genefilter(calls, f1)
  expr = expr[selected,]
  
  # map probe IDs to Ensembl
  map = select(db, featureNames(eset), columns=c("ENSEMBL"))
  # get only first match for duplicates 
  map = map[match(unique(map$PROBEID), map$PROBEID),] 
  
  # collapse probesets to Ensembl IDs by median
  ensembl.list = map$ENSEMBL[match(rownames(expr), map$PROBEID)]
  exprEnsembl = aggregate(expr, by = list(ensembl.list), FUN = median)
  
  # map to human consensus orthologs
  species = ifelse(expt == "E-GEOD-5296", "mouse", "rat")
  orthologs = read.delim(paste("data/orthologs/biomart-ensembl-orthologs-",
                               species, ".tsv", sep=""))
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
  
  # save the matrix 
  save(exprEnsembl, file = paste0("data/expression/diffexpr/", expt, "/", expt, 
                                 "-human.Rdata"))
}
