# Summarise data from mousebrain.org to a data frame where each column
# contains the mean expression level for each cell subtype, at three different
# levels of cell type classification.
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)
library(tidyverse)

# load in main expression matrix
## not committed to github due to size 
df = read.csv("data/expression/mousebrain.org/sc_expr.txt", header = F)

# get list of taxonomy levels
tax_levels = list.files("data/expression/mousebrain.org",
                        pattern = "tax*", full.names = T)

# get genes
genes = read.csv("data/expression/mousebrain.org/sc_genes.txt",
                 header = F) %>% rename(genes = "V1")

for (labels in tax_levels) {
  message("... working on labels: ", labels)
  labs = read.delim(labels, header = F, col.names = "type")
  expr = labs %>% bind_cols(as.data.frame(t(df)))
  sum = expr %>%
    group_by(type) %>%
    summarise_all(funs(mean))
  
  # reshape
  rownames(sum) = sum$type
  
  # transpose back, add in gene names
  final = genes %>% bind_cols(as.data.frame(t(sum[, -1])))
  colnames(final) = c("genes", sum$type)
  write.table(final, paste0("data/expression/mousebrain.org/",
                            gsub(".txt", "", basename(labels)), ".txt"),
              quote = F, row.names = F, sep = "\t")
}
