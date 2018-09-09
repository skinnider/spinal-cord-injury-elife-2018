# Compare the number of co-complex memberships among LC genes to random 
# expectation.
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(igraph)

# read LC genes
dat = read.delim("data/literature_review/Table_S1.txt")
# map to gene IDs 
map = read_tsv("data/identifier/HUMAN_9606_idmapping.dat.gz",
               col_names = c("uniprot", "db", "id"))
ens = filter(map, db == 'Ensembl') %>%
  dplyr::rename(ensembl = id) %>%
  dplyr::select(-db)
ent = filter(map, db == 'GeneID') %>%
  dplyr::select(-db)
ens_ent = left_join(ens, ent, by = 'uniprot') %>% dplyr::select(-uniprot)
dat %<>% left_join(ens_ent, by = 'ensembl') %>%
  drop_na(id)
ids = na.omit(unique(dat$id))

# read complexes
lines = readLines("data/PPI/complexes/huMAPclusters.txt")
split = strsplit(lines, " ")
humap = split %>%
  map(~ as.data.frame(t(combn(., 2)))) %>%
  bind_rows()

# get number of co-complex memberships
g = graph_from_data_frame(humap, directed = F)
nodes = intersect(ids, names(V(g)))
subgraph = induced_subgraph(g, nodes)
cc = length(E(subgraph))

# randomly shuffle complexes
set.seed(0)
bootstraps = 1e3
pos = which(seq_along(vec) %in% (cumsum(lengths(split)) + 1))
rnd_cc = pbapply::pblapply(seq_len(bootstraps), function(x) {
  # shuffle proteins
  vec = sample(unlist(split))
  # re-assemble into complexes
  complexes = split(vec, cumsum(seq_along(vec) %in% pos))
  # get number of co-complex memberships
  pairs = complexes %>%
    map(~ as.data.frame(t(combn(., 2)))) %>%
    bind_rows()
  g = graph_from_data_frame(pairs, directed = F)
  nodes = intersect(ids, names(V(g)))
  subgraph = induced_subgraph(g, nodes)
  length(E(subgraph))
}) %>% unlist()

# write 
results = cbind(source = "hu.MAP", rbind(
  data.frame(class = 'observed', co_complex = cc),
  data.frame(class = 'shuffled', co_complex = rnd_cc)))
# write
write.table(results, "data/PPI/co-complex.txt", quote = F, 
            sep = "\t", row.names = F)
