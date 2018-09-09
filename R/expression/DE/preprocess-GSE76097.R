# Preprocess the Anderson et al. RNA-seq dataset.
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(preprocessCore)
library(org.Mm.eg.db)
source("R/expression/DE/DE-functions.R")

# set expt
expt = "GSE76097"

# read SOFT file
lines = readLines("data/expression/geo/GSE76097/GSE76097_family.soft.gz")

# create design
samples = lines[startsWith(lines, "^SAMPLE")] %>%
  gsub("^.*= ", "", .)
names = lines[startsWith(lines, "!Sample_title")] %>%
  gsub("^.*= ", "", .)
genotype = lines[startsWith(
  lines, "!Sample_characteristics_ch1 = genotype")] %>%
  gsub("^.*: ", "", .)
condition = lines[startsWith(
  lines, "!Sample_characteristics_ch1 = condition")] %>%
  gsub("^.*: ", "", .)
# adjust genotype
genotype = ifelse(genotype == "Wild type", "WT", "STAT3")
# create targets
targets = data.frame(sample = samples, name = names, genotype = genotype,
           condition = condition)
# keep only flow through fraction
targets = targets %>% filter(grepl("FLOW", name))

# read gene count matrix, aggregate ensembl, reshape
expr = read.csv(
  "data/expression/geo/GSE76097/GSE76097_Complete_geneList_flow_only.csv.gz")
expr = expr %>% dplyr::select(X, ends_with("raw_count"))

# map probe IDs to Ensembl
map = select(org.Mm.eg.db, expr$X, keytype = "SYMBOL",
       columns = c("SYMBOL", "ENSEMBL"))
# get only first match for duplicates
map = map[match(unique(map$ENSEMBL), map$ENSEMBL), ]

# collapse probesets to Ensembl IDs by median
ensembl.list = map$ENSEMBL[match(expr[, 1], map$SYMBOL)]
expr = aggregate(expr[, -1], by = list(ensembl.list), FUN = median)
rownames(expr) = expr[, 1]
expr = expr[, -1]

# filter genes expressed only in a single sample
keep = rowSums(expr > 0) > 1
expr = expr[keep, ]

# convert to matrix for processing
mat = as.matrix(expr)

# quantile normalize to average empirical distribution across all samples
mat = normalize.quantiles(mat)

# quantile normalize each gene to standard normal distribution
quantNorm = function(x) qnorm(
 rank(x, ties.method = "average") / (length(x) + 1))
mat.n = apply(mat, 2, quantNorm)

# convert back to data frame to set dimension names
exprEnsembl = as.data.frame(mat.n)
exprEnsembl = cbind(rownames(expr), exprEnsembl)
colnames(exprEnsembl) = c("Gene", colnames(expr))

# load GTEx modules and seed proteins
modules = read.delim("data/modules/GTEx-modules.tsv")

# map to human consensus orthologs
orthologs = read.delim("data/orthologs/biomart-ensembl-orthologs-mouse.tsv")
orthologs = orthologs[orthologs[, 2] != "", ] # remove unmapped orthologs
# pick a single ortholog when multiple match
orthologs = pick.orthologs(orthologs, modules$gene)
exprEnsembl[, 1] = orthologs[match(exprEnsembl[, 1], orthologs[, 1]), 2]
exprEnsembl = exprEnsembl[!is.na(exprEnsembl[, 1]),] # remove NA
exprEnsembl = exprEnsembl[which(exprEnsembl[, 1] != ""),] # remove blank

# aggregate again, if needed
exprEnsembl = aggregate(exprEnsembl[, -1], by = list(exprEnsembl[, 1]),
                       FUN = median)
# re-format as matrix
geneNames = exprEnsembl[, 1]
exprEnsembl = exprEnsembl[, -1]
rownames(exprEnsembl) = geneNames
colnames(exprEnsembl) = targets$sample
exprEnsembl = as.matrix(exprEnsembl)

# save the matrix
save(exprEnsembl, file = paste0(
  "data/expression/geo/GSE76097/", expt, "-human.Rdata"))
# save targets
write.table(targets,
            paste0("data/expression/geo/GSE76097/", expt, "-targets.txt"),
            quote = F, row.names = F, sep = "\t")
