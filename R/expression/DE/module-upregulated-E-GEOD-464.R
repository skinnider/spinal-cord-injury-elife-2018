# Test whether M3 and M7 genes are differentially expressed following spinal 
# cord injury in E-GEOD-464 using limma.
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)
library(genefilter)
library(limma)
library(rgu34a.db)
library(rgu34b.db)
library(rgu34c.db)
source("R/expression/DE/DE-functions.R")

# set expt directory and get targets
expt = "E-GEOD-464" 
expt.dir = paste0("data/expression/DE/", expt)
targets = readTargets(paste0(expt.dir, "/targets.tsv"))

# set weight drop level for each severity
drops = list(normal = 0, mild = 17.5, moderate = 50, severe = 50) # M/S same
targets$drop = unlist(drops[targets$injury])

# read in normalized mas5 data and calls for each array
arrays = c("A-AFFY-18", "A-AFFY-19", "A-AFFY-20")
dbs = c(rgu34a.db, rgu34b.db, rgu34c.db)
for (i in 1:length(arrays)) {
  ## map probesets 
  db = dbs[[i]]
  array = arrays[[i]]
  message("Processing data for array ", array)
  
  # read MAS5-normalized data 
  expr = read.delim(paste0(expt.dir, "/", array, "_mas5.txt.gz"))
  rownames(expr) = expr[, 2]
  expr = expr[, -c(1:2)]
  # read MAS5 P/M/A calls 
  calls = read.delim(paste0(expt.dir, "/", array, "_mas5-calls_pval.txt.gz"))
  rownames(calls) = calls[, 2]
  calls = calls[, -c(1:2)]
  # filter genes expressed in < 20% of samples 
  f1 = function(x) (length(which(x < 0.04)) > 0.8 * length(expr))
  selected = genefilter(calls, f1)
  expr = expr[selected, ]
  
  # get probe matches to ENSEMBL gene ID
  map = select(db, rownames(expr), columns = "ENSEMBL") 
  ensembl.list = map$ENSEMBL[match(rownames(expr), map$PROBEID)]
  ensembl.list[which(is.na(ensembl.list))] = ""
  exprEnsembl = aggregate(expr, by = list(ensembl.list), FUN = median)
  exprEnsembl = exprEnsembl[exprEnsembl[,1] != "", ] # remove NA
  
  ## ortholog mapping 
  # load GTEx modules 
  modules = read.delim("data/modules/GTEx-modules.tsv")
  # map to human consensus orthologs
  orthologs = read.delim("data/orthologs/biomart-ensembl-orthologs-rat.tsv")
  orthologs = orthologs[orthologs[,2] != "", ] # remove unmapped orthologs 
  # pick a single ortholog when multiple match
  orthologs = pick.orthologs(orthologs, modules$gene) 
  exprEnsembl[, 1] = orthologs[match(exprEnsembl[, 1], orthologs[, 1]), 2]
  exprEnsembl = exprEnsembl[!is.na(exprEnsembl[,1]), ] # remove NA
  exprEnsembl = exprEnsembl[which(exprEnsembl[,1] != ""), ] # remove blank
  
  # aggregate again, if needed
  exprEnsembl = aggregate(exprEnsembl[, -1], by=list(exprEnsembl[, 1]), 
                          FUN = median)
  # re-format as matrix
  geneNames = exprEnsembl[, 1]
  exprEnsembl = exprEnsembl[, -1]
  rownames(exprEnsembl) = geneNames
  exprEnsembl = as.matrix(exprEnsembl)
  
  # write human data
  exprFile = paste0("data/expression/DE/", expt, "/", array, "_expr.txt")
  write.table(exprEnsembl, exprFile, sep = "\t", quote = F, row.names = T)
  
  ## before geneSetTest: which genes are significantly correlated to severity?
  colnames(exprEnsembl) = gsub("\\..*$", "", colnames(exprEnsembl))
  targets$filename = gsub("\\..*$", "", targets$filename)
  heights = targets$drop[match(colnames(exprEnsembl), targets$filename)]
  height.cor = suppressWarnings(apply(exprEnsembl, 1, cor.test, heights, 
                                      method = "spearman"))
  height.p = sapply(height.cor, function(x) x$p.value)
  # calculate Delahaye-Duriez et al's gene-level score
  height.rho = sapply(height.cor, function(x) x$estimate)
  height.score = -log10(height.p) * height.rho
  height.df = data.frame(p = height.p, score = height.score)
  # save correlation
  pFile = paste("data/expression/diffexpr/", expt, "/", array, 
                "_height_scores.txt", sep="")
  write.table(height.df, pFile, sep="\t", row.names=T, quote=F, col.names = T)
  assign(paste("P", i, sep = ""), height.df)
  
  ## limma processing
  # create design matrix: injury vs. control at T9
  targets.subset = targets[targets$array == array, ]
  targets.subset$binary = unlist(sapply(targets.subset$injury, function(x) 
    ifelse(x == "normal", "normal", "injury")))
  f1 = paste0(targets.subset$binary, targets.subset$level)
  f1 = factor(f1)
  design1 = model.matrix(~0 + f1)
  colnames(design1) = levels(f1)
  # create second design matrix: injury vs. control at T9 for 1d, 3d, and 7d
  f2 = paste0(targets.subset$binary, targets.subset$level, targets.subset$time)
  f2 = factor(f2)
  design2 = model.matrix(~0 + f2)
  colnames(design2) = levels(f2)
  # fit linear models
  fit = lmFit(exprEnsembl, design1)
  fitTime = lmFit(exprEnsembl, design2)
  # define contrast matrices
  cm1 = makeContrasts(injuryT9 = "injuryT9-normalT9", levels = design1)
  cm2 = makeContrasts(
    E1d = "injuryT91d-normalT91d", 
    E3d = "injuryT93d-normalT93d", 
    E7d = "injuryT97d-normalT97d", 
    E14d = "injuryT914d-normalT914d",
    levels = design2)
  # extract linear model fits for contrasts
  fit = contrasts.fit(fit, cm1)
  fit = eBayes(fit)
  assign(paste0("fit", i), fit)
  fitTime = contrasts.fit(fitTime, cm2)
  fitTime = eBayes(fitTime)
  assign(paste0("fitTime", i), fitTime)
}

## module upregulation analysis 
# combine fits from different chips
fitFinal = rbind(fit1$t, fit2$t, fit3$t)
fitTime = rbind(fitTime1$t, fitTime2$t, fitTime3$t)
fit.list = list(fitFinal, fitTime)
# finally, do gene set enrichment analysis for each module
diff = data.frame(module = integer(), condition = character(), p = numeric(), 
                  direction = character())
for (k in 1:length(fit.list)) {
  fit = fit.list[[k]]
  for (module in unique(modules$module)) {
    if (module == 0)
      next 
    moduleGenes = modules$gene[modules$module == module]
    moduleOnChip = intersect(moduleGenes, rownames(fit))
    conditions = colnames(fit)
    for (i in seq_len(length(conditions))) {
      condition = conditions[i]
      for (alt in c("up", "down")) {
        test = geneSetTest(moduleOnChip, fit[, i], type = "t", 
                           alternative = alt, ranks.only = T)
        diff = rbind(diff, list(module, condition, test, alt))
      }
    }
  } 
  colnames(diff) = c("module", "condition", "p", "direction")
}
diff$p = diff$p * max(modules$module) # Bonferroni correction
diff$sig = diff$p < 0.05
diff = diff[order(diff$module), ]

# write output 
outFile = paste("data/expression/DE", expt, "regulation.tsv", sep = "/")
write.table(diff, outFile, quote = F, row.names = F, sep = "\t")

# module correlation analysis
# combine correlation data with height from different chips
corrs = lapply(1:3, function(i) read.delim(
  paste0("data/expression/DE/", expt, "/", arrays[[i]], 
         "_height_scores.txt")))
corrs = do.call(rbind, corrs)
# analyze module enrichment for correlated or anti-correlated genes 
enrich = data.frame(module = integer(), condition = character(), p = numeric(), 
                    direction = character())
scores = corrs$score
names(scores) = rownames(corrs)
for (module in unique(modules$module)) {
  if (module == 0)
    next 
  moduleGenes = modules$gene[modules$module == module]
  moduleOnChip = intersect(moduleGenes, rownames(corrs))
  for (alt in c("up", "down")) {
    test = geneSetTest(moduleOnChip, scores, type = "t",
                       alternative = alt, ranks.only = T)
    dir = ifelse(alt == "up", "correlated", "anti-correlated")
    enrich = rbind(enrich, list(module, test, dir))
  }
}
colnames(enrich) = c("module", "p", "direction")
enrich$p = p.adjust(enrich$p, method = "bonferroni")
enrich$sig = enrich$p < 0.05
enrich = enrich[order(enrich$module), ]

# write correlation output 
outFile = file.path("data/expression/DE", expt, "correlation-enrichment.tsv")
write.table(enrich, outFile, sep = "\t", quote = F, row.names = F)
