# Use WGCNA to construct a gene coexpression network in spinal cord, and 
# detect coexpressed gene modules, using data from GTEx.
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(limma)
library(WGCNA)

# set input and parameters
file = "data/expression/gtex/spinal_cord_expr_norm.csv.gz"
beta = 5
minModuleSize = 20
mergeCutHeight = 0.2

# read expression data
expr = read.csv(file) %>%
  column_to_rownames('Gene') %>%
  as.matrix() %>% 
  # remove gene version numbers
  set_rownames(gsub("\\..*$", "", rownames(.))) %>%
  # transpose for WGCNA
  t() 

# create networks with `blockwiseModules`
net = blockwiseModules(expr, power = beta, minModuleSize = minModuleSize, 
                       mergeCutHeight = mergeCutHeight, verbose = 3, 
                       corType = "bicor", maxPOutliers = 0.05,
                       networkType = "signed", numericLabels = T)      

# write network
save(net, file = "data/modules/GTEx-network.Rdata")

# write modules
modules = data.frame(gene = colnames(expr), module = net$colors)
write.table(modules, "data/modules/GTEx-modules.tsv", sep = "\t", quote = F,
            row.names = F)

# write expression matrix
save(expr, file = "data/modules/GTEx-expression.Rdata")
