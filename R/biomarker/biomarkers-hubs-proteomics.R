# Analyze the capacity of M3 hub genes to stratify rats based on injury
# severity in the proteomics dataset. 
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)
library(limma)
library(MASS)
library(ROCR)
library(magrittr)

# read proteomics data
dat = read.delim("data/proteomics/proteomics-processed.txt")
genes = dat$Ortholog
dat = dat[, -1]
rownames(dat) = genes
# convert to matrix
expr = as.matrix(dat)
# remove proteins found in two or fewer samples
detections = rowSums(expr > 0)
expr = expr[detections > 2, ]

# get targets
targets = read.delim("data/proteomics/targets.txt")

# read modules
modules = read.delim("data/modules/GTEx-modules.tsv")

# read M3 hub genes, defined as the top 10% by ME3 connectivity 
modules = read.delim("data/modules/GTEx-modules.tsv")
M3_genes = modules$gene[modules$module == 3]
kME = read.delim("data/modules/kME.txt")
M3_idxs = rownames(kME) %in% M3_genes
M3_kME = kME$M3[M3_idxs]
pct_hubs = 0.9
hub_idxs = which(M3_kME >= quantile(M3_kME, probs = pct_hubs))
hub_names = rownames(kME[M3_idxs,][hub_idxs,])
hub_kMEs = kME[M3_idxs,][hub_idxs,]$M3
hubs = setNames(hub_kMEs, hub_names)
# subset to hubs detected in the proteomics data
hubs = hubs[names(hubs) %in% rownames(expr)]
hubs = hubs[order(-hubs)]

# perform linear discriminant analysis for mild/severe injuries 
results = data.frame()
rocs = data.frame()
for (hub in names(hubs)) {
  message("analyzing hub ", hub, " ...")
  # get rank
  rank = which(names(hubs) == hub)
  labels = targets$group
  input = data.frame(intensity = expr[hub, ], label = labels)
  # train model
  fit = MASS::lda(formula = label ~ intensity, data = input, 
                  na.action = "na.omit", CV = T)
  # get results
  ct = table(input$label, fit$class)
  # calculate total % correctly classified
  correct = sum(diag(prop.table(ct)))
  # calculate AUCs
  for (outcome in c("100", "200")) {
    binary_labels = as.character(labels) == outcome
    pred = prediction(fit$posterior[, outcome], binary_labels) 
    auc = performance(pred, "auc")@y.values
    roc = performance(pred, "tpr", "fpr")
    # add to results
    results %<>% rbind(list(hub = hub, rank = rank, severity = outcome, 
                            accuracy = correct, auc = unlist(auc)))
    # add to ROCs
    rocs %<>% rbind(data.frame(hub = hub, rank = rank, severity = outcome, 
                               x = unlist(roc@x.values), 
                               y = unlist(roc@y.values)))
  }
}

# write results
write.table(results, "data/biomarker/hubs-proteomics.txt",
            quote = F, row.names = F, sep = "\t")
