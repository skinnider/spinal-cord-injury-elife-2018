# Compare the number of protein-protein interactions (PPIs) and the size of the
# largest connected component (LCC) formed by LC genes in human interactomes.
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(igraph)

# read LC genes
dat = read.delim("data/literature_review/Table_S1.txt")
# map to UniProt 
map = read_tsv("data/identifier/HUMAN_9606_idmapping.dat.gz",
               col_names = c("uniprot", "db", "id"))
ens = filter(map, db == 'Ensembl') %>%
  dplyr::rename(ensembl = id) %>%
  dplyr::select(-db)
dat %<>% left_join(ens, by = 'ensembl')
proteins = na.omit(unique(dat$uniprot))
# also map to Entrez
ent = filter(map, db == 'GeneID') %>%
  dplyr::select(-db)
ens_ent = left_join(ens, ent, by = 'uniprot') %>% dplyr::select(-uniprot)
dat %<>% left_join(ens_ent, by = 'ensembl')
ids = na.omit(unique(dat$id))

# list interactomes
interactome_files = list.files("data/PPI/interactomes", full.names = T,
                               pattern = "*.gz")

# analyze each separately 
bootstraps = 1e3
all_results = data.frame()
for (interactome_file in interactome_files) {
  interactome = gsub("-interactome.*$", "", basename(interactome_file))
  message("analyzing interactome: ", interactome)
  set.seed(0)
  
  # read PPI
  g = graph_from_data_frame(read.delim(interactome_file, comment.char = '#'))
  # calculate observed # of PPIs and LCC
  nodes = unique(intersect(proteins, names(V(g))))
  if (length(nodes) == 0)
    nodes = unique(intersect(ids, names(V(g))))
  subgraph = induced_subgraph(g, nodes)
  obs_PPI = length(E(subgraph))
  clusters = clusters(subgraph)
  obs_LCC = max(clusters$csize)
  
  # rewire networks
  iterations = 6.9077 * length(E(g))
  rewires = pbapply::pblapply(seq_len(bootstraps), function(x) rewire(
    g, keeping_degseq(loops = F, niter = iterations)))
  
  # calculate rewired PPI and LCC
  subgraphs = map(rewires, ~ induced_subgraph(., nodes))
  rnd_PPI = map_int(subgraphs, ~ length(E(.)))
  rnd_LCC = map_dbl(subgraphs, ~ max(clusters(.)$csize))
  
  # calculate empirical P values
  p_LCC = mean(rnd_LCC >= obs_LCC)
  p_PPI = mean(rnd_PPI >= obs_PPI)
  
  # add to results
  results = cbind(interactome = interactome, rbind(
    data.frame(class = 'observed', PPI = obs_PPI, LCC = obs_LCC),
    data.frame(class = 'rewired', PPI = rnd_PPI, LCC = rnd_LCC)))
  all_results %<>% rbind(results)
}

# write
write.table(all_results, "data/PPI/PPI-LCC.txt", quote = F, 
            sep = "\t", row.names = F)
