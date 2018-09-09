# Calculate all shared terms between pairs of literature-curated proteins.
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(flavin)

# read LC dataset
dat = read.delim("data/literature_review/Table_S1.txt")
lc_genes = unique(dat$ensembl)

# read GO terms
go = read.delim("data/GO/hsa-proteins.txt.gz")
all_genes = unique(go$ensembl)
G = length(all_genes)

# are LC and random proteins annotated with different # of terms?
x = go$breadth
y = go$breadth[go$ensembl %in% lc_genes]
wilcox.test(x, y, alternative = 'greater')

results = data.frame()
# iterate over breadths
breadths = c(100, 30, 10)
for (max_breadth in breadths) {
  message("analyzing shared GO terms at breadth = ", max_breadth, " ...")
  go0 = filter(go, breadth < max_breadth)
  # iterate over ontologies
  for (root in c("CC", "BP", "MF")) {
    message("  analyzing shared ", root, " terms ...")
    go1 = filter(go0, ontology == root) %>%
      dplyr::select(ensembl, term) %>% 
      distinct()
    # convert to annotations
    ann = as_annotation_list(go1, "term", "ensembl")
    # get shared terms
    shared = shared_annotations(ann)
    # get SCI pairs
    pairs = as.data.frame(t(combn(lc_genes, 2)))
    # run FET
    proteome.yes = sum(shared[upper.tri(shared)] > 0, na.rm = T)
    proteome.no = (G * (G - 1) / 2) - proteome.yes
    sci.yes = network_shared_annotations(pairs, shared)$shared %>%
      replace(., is.na(.), 0) %>%
      sum(. > 0)
    sci.no = nrow(pairs) - sci.yes
    mat = matrix(c(sci.yes, proteome.yes - sci.yes, sci.no, proteome.no), 
                  nrow = 2)
    f = fisher.test(mat)
    p = f$p.value * 12 # Bonferroni correction
    # calculate enrichment
    enr = (sci.yes / (sci.yes + sci.no)) / 
      (proteome.yes / (proteome.yes + proteome.no))
    # append to results
    results %<>% rbind(data.frame(ontology = root, breadth = max_breadth,
                                  enrichment = enr, p = p))
  }
}

# write 
write.table(results, "data/GO/LC-GO-enrichment.txt", quote = F, 
            sep = "\t", row.names = F)
