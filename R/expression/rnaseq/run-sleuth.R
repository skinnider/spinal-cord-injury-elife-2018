# Analyze severity-dependent differential expression of spinal cord modules
# in West lab RNA-seq data using sleuth. 
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)
library(sleuth)
library(biomaRt)
library(dplyr)
library(limma)
library(tidyverse)
source("R/expression/DE/DE-functions.R")

# define salmon dir for wasabi results
salmonDir = "data/expression/rnaseq/salmon"
# get wasabi converted files for DE using sleuth
files = list.dirs(salmonDir, recursive = F)
# append wasabi output
was_dirs = file.path(files, "abundance.h5")

# ensure S9 is removed
was_dirs = was_dirs[!grepl("S9", was_dirs)]

# create targets table
s2c = read.delim("data/expression/rnaseq/targets.txt") 
# append file paths to auxillary table
s2c$path = was_dirs

# read transcript to gene map
t2g = read.delim("data/expression/rnaseq/tx2gene/tx2gene-rnor-89.txt",
                 col.names = c("target_id", "ens_gene"))

# prepare sleuth object
so = sleuth_prep(s2c, target_mapping = t2g, aggregation_column = 'ens_gene')

# fit full model
so = sleuth_fit(so, ~ group, 'full')

# fit reduced model
so = sleuth_fit(so, ~ 1, 'reduced')

# perform liklihood ratio test
so = sleuth_lrt(so, 'reduced', 'full')

# perform wald test
so = sleuth_wt(so, 'group')

# get results
results = sleuth_results(so, 'group', test_type = 'wt', show_all = F)

# map to human consensus orthologs
orthologs = read.delim("data/orthologs/biomart-ensembl-orthologs-rat.tsv")
orthologs = orthologs[orthologs[,2 ] != "", ] # remove unmapped orthologs
# pick a single ortholog when multiple match
orthologs = pick.orthologs(orthologs, modules$gene)
matches = orthologs$MajorityOrtholog[match(results$target_id, orthologs[, 1])]
matches[is.na(matches)] = ""
results$ens_human = matches
results = results[results$ens_human != "", ]

# save results
write.table(results, "data/expression/rnaseq/sleuth/sleuth-severity.txt",
            quote = F, row.names = F, sep = "\t")

# also write sleuth's normalized data
norm = so$obs_norm_filt %>%
  dplyr::select(-scaled_reads_per_base) %>%
  spread('sample', 'tpm') %>%
  column_to_rownames('target_id') %>% 
  as.matrix()
# map to human consensus orthologs
orthologs = read.delim("data/orthologs/biomart-ensembl-orthologs-rat.tsv")
orthologs = orthologs[orthologs[, 2] != "", ] # remove unmapped orthologs
# pick a single ortholog when multiple match
orthologs = pick.orthologs(orthologs, modules$gene)
matches = orthologs$MajorityOrtholog[match(rownames(norm), orthologs[, 1])]
matches[is.na(matches)] = ""
rownames(norm) = matches
norm = norm[rownames(norm) != "",]
# sum TPM 
norm %<>%
  as.data.frame() %>% 
  rownames_to_column('gene') %>%
  group_by(gene) %>%
  dplyr::summarise_if(is.numeric, sum) %>%
  ungroup() %>%
  column_to_rownames('gene')
# write 
write.table(norm, "data/expression/rnaseq/sleuth/sleuth-norm-filt.txt",
            quote = F, sep = "\t")

# take most significant result per human gene
sleuth_sub = results %>%
  group_by(ens_human) %>%
  dplyr::slice(1) %>%
  ungroup() 

# load GTEx modules
modules = read.delim("data/modules/GTEx-modules.tsv")

# finally, do gene set enrichment analysis for each module
diff = data.frame()
for (module in unique(modules$module)) {
  if (module == 0)
    next
  moduleGenes = modules$gene[modules$module == module]
  moduleSleuth = intersect(moduleGenes, sleuth_sub$ens_human)
  for (alt in c("up", "down")) {
    test = geneSetTest(moduleSleuth,
                       setNames(sleuth_sub$b, sleuth_sub$ens_human),
                       type = "t", alternative = alt, ranks.only = T)
    diff = rbind(diff, list(module, test, alt))
  }
}
colnames(diff) = c("module", "p", "direction")
diff$p = diff$p * 2 * length(unique(modules$module)) # Bonferroni correction
diff$sig = diff$p < 0.05
diff = diff[order(diff$p), ]

# write output
write.table(diff, 
            "data/expression/rnaseq/sleuth/sleuth-severity-module-regulation.txt",
            sep = "\t", quote = F, row.names = F)
