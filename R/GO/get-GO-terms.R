# Retrieve GO terms from Bioconductor for all literature-curated proteins.
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)
library(GO.db)
library(org.Hs.eg.db)
library(tidyverse)
library(magrittr)

db = org.Hs.eg.db
name = "hsa"

# get GO terms for each protein
go = select(db, keys = keys(db), columns = c("GO", "ENSEMBL"))
go = subset(go, !is.na(GO) & !is.na(ENSEMBL)) # ignore NAs

# remove unknown terms and roots
unknown = c("GO:0000004", "GO:0005554", "GO:0008372")
roots = c("GO:0003673", "GO:0005575", "GO:0008150 ")
go = go[!(go$GO %in% unknown | go$GO %in% roots),]

# remove Entrez ID and re-order
go = go[, c("ENSEMBL", "GO", "EVIDENCE", "ONTOLOGY")]

# calculate breadth
breadth.df = go %>%
  group_by(GO) %>%
  summarise(breadth = n_distinct(ENSEMBL)) %>%
  ungroup()
go %<>% left_join(breadth.df, by = 'GO')

# write
colnames(go) = c("ensembl", "term", "evidence", "ontology", "breadth")
goFile = paste("data/GO/", name, "-proteins.txt", sep="")
write.table(go, goFile, sep = "\t", quote = F, row.names = F)
system(paste("gzip --force", goFile))
