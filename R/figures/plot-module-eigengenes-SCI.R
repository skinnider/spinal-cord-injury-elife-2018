# Plot the expression of the M3 eigengene in GSE464, GSE69334, and our own
# RNA-seq and proteomics datasets.
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(limma)
library(WGCNA)
source("R/theme.R")

# create container for all data
all_dat = data.frame()

# GSE464
##########
expt = "E-GEOD-464"
expt.dir = paste0("data/expression/DE/", expt)
targets = readTargets(paste0(expt.dir, "/targets.tsv"))
targets$filename = gsub("\\..*$", "", targets$filename)
# set weight drop level for each severity (combine moderate/severe d/t outcomes)
drops = list(normal = 0, mild = 17.5, moderate = 50, severe = 50)
targets$drop = unlist(drops[targets$injury])

# read in normalized, human-mapped data
arrays = c("A-AFFY-19", "A-AFFY-20")
exprs = lapply(arrays, function(array) {
  # Load human data
  exprFile = paste0("data/expression/DE/", expt, "/", array, "_expr.txt")
  expr = read.delim(exprFile)
  colnames(expr) = gsub("\\..*$", "", colnames(expr))
  return(expr)
})
names(exprs) = arrays

# read M3 genes 
modules = read.delim("data/modules/GTEx-modules.tsv")
M = seq_len(max(modules$module))

# calculate module eigengene
samples = unlist(sapply(exprs, colnames))
# create XY dataframe
xy = data.frame(sample = character(), chip = character(), ME3 = numeric(), 
                 drop = numeric())
# calculate expression
for (name in names(exprs)) {
  expr = exprs[[name]]
  expr = expr[,which(names(expr) %in% samples)]
  colors = modules$module[match(rownames(expr), modules$gene)]
  ME = moduleEigengenes(t(expr), colors)
  drops = targets$drop[match(colnames(expr), targets$filename)]
  df = data.frame(sample = colnames(expr), chip = name, 
                   ME3 = ME$eigengenes$ME3, drop = drops)
  xy = rbind(xy, df)
}
xy$array = targets$array[match(xy$sample, targets$filename)]
xy$level = targets$level[match(xy$sample, targets$filename)]
xy$time = targets$time[match(xy$sample, targets$filename)]
xy$drop = factor(xy$drop)
site = xy[xy$level == "T9" & xy$time == "7d",]

# append to all_dat
site %<>% 
  dplyr::select(drop, ME3) %>%
  dplyr::mutate(drop = fct_recode(drop, sham = '0', moderate = '17.5',
                                  severe = '50')) %>%
  dplyr::rename(group = drop) %>%
  dplyr::mutate(facet = 'GSE464')
str(with(site, cor.test(as.integer(group), ME3, method = 'spearman')))
all_dat %<>% rbind(site)

# RNA-seq
##########
# read normalized and filtered expression
expr = read.delim("data/expression/rnaseq/sleuth/sleuth-norm-filt.txt")

# get targets
targets = data.frame(sample = colnames(expr),
                     group = gsub(".*\\.", "", gsub(
                       "\\.Spinal.*", "", colnames(expr))))
targets$label = ifelse(targets$group == 0, "sham", ifelse(
  targets$group == 100, "moderate", "severe"))

# read modules
modules = read.delim("data/modules/GTEx-modules.tsv")
# filter to genes in modules
overlap = intersect(modules$gene, rownames(expr))
modules %<>% dplyr::filter(gene %in% overlap)
expr = expr[overlap, ]

# calculate eigengenes
colors = modules$module[match(rownames(expr), modules$gene)]
ME = moduleEigengenes(t(expr), colors)

# create frame to plot
df = data.frame(group = factor(targets$label, levels = c(
  "sham", "moderate", "severe")),
  ME3 = ME$eigengenes$ME3)
str(with(df, cor.test(as.integer(group), ME3, method = 'spearman')))

# add to all data
df$facet = "RNA-seq"
all_dat %<>% rbind(df)

# proteomics
##########
# read human-mapped data 
expr = read.delim("data/proteomics/proteomics-processed.txt")
genes = expr$Ortholog
expr = expr[,- 1]
rownames(expr) = genes

# read targets
targets = read.delim("data/proteomics/targets.txt")

# read M3 genes
modules = read.delim("data/modules/GTEx-modules.tsv")

# get eigengene
colors = modules$module[match(rownames(expr), modules$gene)]
ME = moduleEigengenes(t(expr), colors)

# create frame to plot
df = data.frame(group = factor(
  targets$label, levels = c("sham", "moderate", "severe")),
  ME3 = ME$eigengenes$ME3)

# add to all data
df$facet = "Proteomics"
all_dat %<>% rbind(df)

# NT-3
##########
expt = "E-GEOD-69334" 
expt.dir = file.path("data/expression/DE", expt)
targets = readTargets(file.path(expt.dir, "targets.tsv"))
# load expression data
load(paste0("data/expression/DE/", expt, "/", expt, "-human.Rdata"))
expr = exprEnsembl # Rename expression dataframe 
# read modules
modules = read.delim("data/modules/GTEx-modules.tsv")
M = seq_len(max(modules$module))
# calculate module eigengenes
colors = modules$module[match(rownames(expr), modules$gene)]
ME = moduleEigengenes(t(expr), colors)$eigengenes
# get eigengene data for ME3 
ME$sample = colnames(expr)
ME$time = targets$time[match(ME$sample, targets$filename)]
ME$group = targets$group[match(ME$sample, targets$filename)]
ME$tissue = targets$tissue[match(ME$sample, targets$filename)]
ME3 = ME[, c("sample", "time", "group", "tissue", "ME3")]
ME3 = ME3[ME3$group != "Uninjured",]
ME3 = ME3[ME3$tissue == "Middle",]
ME3 = ME3[ME3$group != "Tube",]
df = ME3 %>%
  dplyr::select(group, ME3) %>%
  dplyr::mutate(facet = "NT-3", group = as.character(forcats::fct_recode(
    group, 'NT-3' = 'chitosan')))
x = df$ME3[df$group == 'NT-3']
y = df$ME3[df$group == 'Control']
wilcox.test(x, y)
all_dat %<>% rbind(df)

# plot
all_dat$group %<>% as.character() %>%
  Hmisc::capitalize()
all_dat$group %<>% factor(levels = c(
  'Sham', 'Moderate', 'Severe', 'Control', 'NT-3'))
all_dat$facet %<>% factor(levels = c('GSE464', 'RNA-seq', 'Proteomics', 'NT-3'))
p = ggplot(all_dat, aes(x = group, y = ME3, color = group, fill = group)) +
  facet_wrap(~ facet, scales = 'free', ncol = 4) + 
  geom_boxplot(alpha = 0.4, outlier.size = 0, outlier.shape = NA) +
  geom_jitter(height = 0, width = 0.2, size = 0.4) +
  scale_x_discrete("Injury") + #, labels = c("Sham", "Mild", "Severe", "Control", "NT-3")) +
  scale_color_manual(values = ryb8[c(1, 3, 8, 8, 1)], guide = F) +
  scale_fill_manual(values = ryb8[c(1, 3, 8, 8, 1)], guide = F) +
  sci_theme + 
  theme(axis.line = element_line(color = 'grey50'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        strip.background = element_rect(fill = 'grey40', color = 'grey40'),
        strip.text = element_text(color = 'white', size = 8))
p
ggsave("figures/figure-5E-H.pdf", p, width = 13.5, height = 4.5, units = 'cm',
       useDingbats = F)
