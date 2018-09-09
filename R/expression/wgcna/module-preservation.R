# Evaluate spinal cord gene coexpression module preservation across datasets
# (human) and conservation across species (mouse, rat). 
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)
library(WGCNA)
library(magrittr)

# load GTEx expression data 
load("data/modules/GTEx-expression.Rdata")
gtex = expr

# load GTEx modules 
modules = read.table("data/modules/GTEx-modules.tsv", header = T) %>%
  extract(, 2)

arrays = 43:45
results = data.frame()
for (array in arrays) {
  arrayName = paste("A-AFFY", array, sep="-")
  exprFile = paste("data/expression/arrayexpress/expression/", arrayName,
                    "_expr.csv.gz", sep="")
  expr = read.csv(exprFile)

  # preprocess for WGCNA 
  geneNames = expr[, 1]
  expr = expr[, -1] # remove gene names
  expr = t(expr)  # transpose for WGCNA analysis
  colnames(expr) = geneNames
  
  # set up module preservation 
  setLabels = c("GTEx", arrayName)
  multiExpr = list(GTEx = list(data = gtex), Affy = list(data = expr))
  multiColor = list(GTEx = modules)
  
  # run module preservation 
  mp = modulePreservation(multiExpr, multiColor, referenceNetworks = 1,
                          nPermutations = 100, randomSeed = 1, quickCor = 0,
                          verbose = 3)
  
  # produce Z-summary table
  ref = 1
  test = 2
  statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], 
                    mp$preservation$observed[[ref]][[test]][, -1])
  statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], 
                  mp$preservation$Z[[ref]][[test]][, -1])
  table = cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
                 signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2))
  table$module = rownames(table)
  table$array = arrayName
  table = table[,c("array", "module", "Zsummary.pres", "medianRank.pres",
                    "medianRank.qual", "Zsummary.qual")] # reorder
  
  # add to results 
  results = rbind(results, table)
}

# write results data frame 
write.table(results, "data/modules/preservation/preservation.tsv", sep = "\t",
            quote = F, row.names = F)
