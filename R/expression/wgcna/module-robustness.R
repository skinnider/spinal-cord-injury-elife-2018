# Test how robust enrichment of one or more modules for SCI proteins is to 
# removal of seed proteins. 
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors=F)
library(data.table)
source("R/expression/wgcna/wgcna-functions.R")
source("R/theme.R")
set.seed(0)
bootstraps = 100

# load expression data 
load("data/modules/GTEx-expression.Rdata")
# load coexpression network
load("data/modules/GTEx-network.Rdata")
# read SCI proteins
sci = read.delim("data/literature_review/Table_S1.txt")
sci.genes = unique(sci$ensembl)
# read module genes
modules = read.delim("data/modules/GTEx-modules.tsv")
module.genes = getModuleGenes(modules)

# find significant modules 
enrich = module.enrichment(sci.genes, net, expr, correct="bonferroni")
significant = enrich$module[enrich$enrich.p < 0.05]

# set up results table
results = data.table(condition = character(), module = integer(), 
                      bootstrap = integer(), idx = integer(), p = numeric())

# calculate enrichment in SCI modules for SCI proteins with removal of random 
# proteins from the set of seed genes. 
nGenes = length(sci.genes)
message("Analyzing robustness of module enrichment to removal of ", nGenes, 
        " SCI proteins")
for (i in seq_len(bootstraps)) {
  message("... running bootstrap ", i, " of ", bootstraps)
  genesCopy = sci.genes
  randoms = sample(genesCopy, nGenes)
  for (j in seq_len(nGenes)) {
    gene = randoms[j]
    genesCopy = genesCopy[genesCopy != gene]
    module.df = module.enrichment(genesCopy, net, expr,
                                   correct = "bonferroni")
    df = data.frame(condition = "removal", module = module.df$module, 
                     bootstrap = i, idx = j, p = module.df$enrich.p)
    results = rbind(results, df)
  }
}


# Calculate enrichment in SCI modules for SCI proteins with addition of random
# proteins to the set of seed genes. 
message("Analyzing robustness of module enrichment to addition of ", nGenes, 
        " random proteins")
nonSci = modules$gene[!(modules$gene %in% sci.genes)]
for (i in seq_len(bootstraps)) {
  message("... running bootstrap ", i, " of ", bootstraps)
  genesCopy = sci.genes
  randoms = sample(nonSci, nGenes)
  for (j in seq_len(nGenes)) {
    gene = randoms[j]
    genesCopy = c(genesCopy, gene)
    module.df = module.enrichment(genesCopy, net, expr,
                                   correct = "bonferroni")
    df = data.frame(condition = "addition", module = module.df$module, 
                     bootstrap = i, idx = j, p = module.df$enrich.p)
    results = rbind(results, df)
  }
}

# write
colnames(results) = c("condition", "module", "bootstrap", "idx", "p")
write.table(results, "data/modules/robustness.tsv", sep="\t", quote=F, 
            row.names=F)
