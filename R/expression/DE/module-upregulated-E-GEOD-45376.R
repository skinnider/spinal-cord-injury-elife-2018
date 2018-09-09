# Analyze module regulation in RNA-Seq data from E-GEOD-45376.
setwd("~/git/spinal-cord-injury-elife-2018")
source("R/Settings.R")
library(limma)
library(biomaRt)
library(pSI)
source("R/scripts/expression/diffexpr/DiffExprFunctions.R")

# define experiment
expt = "E-GEOD-45376"

# read input 
dir = "data/expression/diffexpr/E-GEOD-45376/raw"
files = list.files(dir, pattern = "*.csv", full.names = T)
exprs = list()
for (file in files) {
  expr = suppressWarnings(read.csv(
    file, row.names = 1, col.names = basename(file), header = T))
  exprs[[length(exprs) + 1]] = expr
}
expr = do.call(cbind, exprs)

# create a matrix
data = as.matrix(expr)

# get probe matches to ENSEMBL gene ID
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
mySymbols = rownames(data) # mySymbols is a vector of MGI symbols.
map.database = getBM(attributes = c("ensembl_gene_id", "mgi_symbol"), 
                      filters = "mgi_symbol", values = mySymbols, mart = mouse)
map = map.database[match(unique(map.database$mgi_symbol), row.names(data)), ]
map = map[-which(is.na(map$mgi_symbol)), ] # remove NA

# collapse probesets to ENSEMBL IDs by median
ensembl.list = map$ensembl_gene_id[match(row.names(data), map$mgi_symbol)]
exprEnsembl = aggregate(data, by = list(ensembl.list), FUN = median)

# load GTEx modules
modules = read.delim("data/modules/GTEx-modules.tsv")

# map to human consensus orthologs
orthologs = read.delim("data/orthologs/biomart-ensembl-orthologs-mouse.tsv")
orthologs = orthologs[orthologs[,2] != "", ] # remove unmapped orthologs 
# pick a single ortholog when multiple match
orthologs = pick.orthologs(orthologs, modules$gene) 
exprEnsembl[, 1] = orthologs[match(exprEnsembl[, 1], orthologs[, 1]), 2]
exprEnsembl = exprEnsembl[!is.na(exprEnsembl[, 1]) & 
                            exprEnsembl[, 1] != "", ] # remove NA and blank

# aggregate again, if needed
exprEnsembl = aggregate(exprEnsembl[, -1], by = list(exprEnsembl[, 1]), 
                         FUN = median)
# re-format as matrix
geneNames = exprEnsembl[, 1]
exprEnsembl = exprEnsembl[, -1]
rownames(exprEnsembl) = geneNames
exprEnsembl = as.matrix(exprEnsembl)

# save the matrix 
save(exprEnsembl, file = paste0("data/expression/DE/", expt, "/", expt, 
                                "-human.Rdata"))

# read targets 
targets = readTargets(file.path("data/expression/DE", expt, "targets.txt"))

# create design matrices
f1 = factor(targets$group) # injury vs. control
f2 = factor(paste0(targets$group, targets$time)) # as above + time
d1 = model.matrix(~0+f1)
colnames(d1) = levels(f1)
d2 = model.matrix(~0+f2)
colnames(d2) = levels(f2)
# define a contrast matrix
cm1 = makeContrasts(injury="injury-uninjured", levels=d1)
cm2 = makeContrasts(Ed2="injury2day-uninjuredcontrol", 
                    Ed7="injury7day-uninjuredcontrol", levels=d2)
# create list of contrast matrices
matrix.list = list(cm1, cm2)
design.list = list(d1, d2)

# finally, do geneSetTest for each module
diff = data.frame(module = integer(), condition = character(), p = numeric(), 
                   direction = character())
for (i in 1:length(matrix.list)) {
  # Extract linear model fit for contrasts
  contrast.matrix = matrix.list[[i]]
  design = design.list[[i]]
  # Fit the linear model
  fit = lmFit(exprEnsembl, design)
  fit2 = contrasts.fit(fit, contrast.matrix)
  fit2 = eBayes(fit2)
  for (module in unique(modules$module)) {
    if (module == 0)
      next 
    moduleGenes = modules$gene[modules$module == module]
    moduleOnChip = intersect(moduleGenes, rownames(exprEnsembl))
    conditions = colnames(fit2$t)
    for (i in seq_len(length(conditions))) {
      condition = conditions[i]
      for (alt in c("up", "down")) {
        test = geneSetTest(moduleOnChip, fit2$t[, i], type = "t", 
                            alternative = alt, ranks.only = T)
        diff = rbind(diff, list(module, condition, test, alt))
      }
    }
  }
  colnames(diff) = c("module", "condition", "p", "direction")
}
diff$p = diff$p * 2 * max(modules$module) # Bonferroni correction
diff$sig = diff$p < 0.05
diff = diff[order(diff$module), ]

# write output 
outFile = file.path("data/expression/DE", expt, "regulation.tsv")
write.table(diff, outFile, sep = "\t", quote = F, row.names = F)
