# Preprocess proteomics data for the validation of gene co-expression modules
# in spinal cord injury.
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)
library(tidyverse)
library(pipeR)

# read raw data
dat = read.delim("data/proteomics/maxquant/proteinGroups.txt.gz")

# load rat uniprot-Ensembl map and consensus orthologs
rat_map = read.delim("data/identifier/rat-uniprot-Ensembl.tsv")
orthologs = read.delim("data/orthologs/biomart-ensembl-orthologs-rat.tsv",
                       col.names = c("Ensembl", "Ortholog"), strip.white = T)
orthologs = orthologs %>%
  unnest(Ortholog = strsplit(Ortholog, " "))
# get sample ids
samples = colnames(dat)[grepl("LFQ.intensity.", colnames(dat))]

# subset to protein IDs and LFQ intensities for each sample
dat = dat[, colnames(dat) %in% c("Protein.IDs", samples)]

# map genes to rat ENSEMBL, then human ENSEMBL, and summarise intensity to
# gene level
dat = dat %>>%
  unnest(UniProt = strsplit(Protein.IDs, ";")) %>>%
  left_join(rat_map) %>>%
  dplyr::filter(!is.na(Ensembl)) %>>%
  left_join(orthologs) %>>%
  dplyr::filter(Ortholog != "") %>>%
  dplyr::select(-Protein.IDs, -UniProt, -Ensembl) %>>%
  (~dat_tmp) %>>%
  dplyr::select(starts_with("LFQ")) %>>%
  rowwise() %>>%
  do(data.frame(variance = var(unlist(.))/mean(unlist(.)))) %>>%
  bind_cols(dat_tmp, .) %>>%
  group_by(Ortholog) %>>%
  dplyr::slice(which.max(variance)) %>>%
  dplyr::select(-variance, -variance1, -Ortholog)

# simplify to data frame and write
dat = as.data.frame(dat)
write.table(dat, "data/proteomics/proteomics-processed.txt", quote = F,
            sep = "\t", row.names = F)
