# Map rat transcripts to genes.
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)
library(ensembldb)

# create ensembldb object from archived GTF file 
DB = ensDbFromGtf(
  gtf = "data/expression/rnaseq/tx2gene/Rattus_norvegicus.Rnor_6.0.89.gtf.gz")
EDB = EnsDb(DB)

# make transcripts data frame 
transcripts = transcripts(EDB, return.type = "data.frame")
tx2gene = transcripts[, c("tx_name", "gene_id")]

# save
write.table(tx2gene, "data/expression/rnaseq/tx2gene/tx2gene-rnor-89.txt", 
            row.names = F, quote = F, sep = "\t")
