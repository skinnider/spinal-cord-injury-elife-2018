# Plot the network of module eigengenes in the spinal cord gene 
# coexpression network as a tree.  
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)
library(WGCNA)
library(ggplot2)
library(ggdendro)
source("R/theme.R")

# load GTEx spinal cord gene expression 
load("data/modules/GTEx-expression.Rdata")

# load GTEx modules 
modules = read.delim("data/modules/GTEx-modules.tsv")

# calculate eigengenes
ME = moduleEigengenes(expr, modules$module)

# calculate correlation between M3 and M7
cor = cor.test(ME$eigengenes$ME3, ME$eigengenes$ME7, method = "spearman")
message("Spearman's rho = ", format(cor$estimate, digits = 2), ", P = ", 
        format(cor$p.value, digits = 2))

# Calculate dissimilarity matrix
dissimME = (1 - t(cor(ME$eigengenes))) / 2
hclustdatME = hclust(as.dist(dissimME), method = "average")

# Plot the eigengene dendrogram
# Convert cluster object to use with ggplot
dendr = dendro_data(hclustdatME, type="rectangle") 
p = ggplot() + 
  geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend),
               color="grey50") + 
  geom_text(data=label(dendr), aes(x=x, y=y-0.02, label=label, hjust=0), size=2) +
  coord_flip() + 
  scale_y_reverse(expand=c(0.2, 0)) + 
  sci_theme +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())
p
ggsave("figures/figure-3D.pdf", p, width = 6, height = 4.5, 
       useDingbats = F, units = "cm")
