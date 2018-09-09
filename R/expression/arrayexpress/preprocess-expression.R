# Preprocess all microarray samples by normalizing using the MAS5 algorithm,
# then adjusting for batch effects using ComBat.
setwd("~/git/spinal-cord-injury-elife-2018")
library(affy)
library(sva)
library(rat2302.db)
library(mouse4302.db)
library(hgu133plus2.db)
library(limma)

arrays = c("A-AFFY-43", "A-AFFY-45", "A-AFFY-44")
dbs = list(rat2302.db, mouse4302.db, hgu133plus2.db)
dbs = setNames(dbs, arrays)

for (array in arrays) {
  message("preprocessing expression files for array ", array, " ...")
  
  # Get the paths of all CEL files
  message(".. reading CEL files...")
  files = character()
  cel.dir = paste("data/expression/arrayexpress/cel/", array, sep = "")
  accession.dirs = list.dirs(cel.dir, recursive = F)
  for (accession.dir in accession.dirs) {
    cels = list.files(accession.dir, pattern = "*.cel", full.names = T, 
                      ignore.case = T)
    files = c(files, cels)
  }
  
  # create AffyBatch
  affy.batch = read.affybatch(files)
  
  # run MAS5
  message(".. normalizing CEL files using MAS5 ...")
  eset = mas5(affy.batch)
  expr = exprs(eset)
  calls = mas5calls(affy.batch)
  calls = exprs(calls)
  
  # read batches
  samples.file = paste("data/expression/arrayexpress/samples/", 
                       array, "-samples.csv", sep = "")
  batches = read.csv(samples.file, header = F, 
                     col.names = c("batch", "sample"))
  
  # create batches vector for ComBat 
  samples = colnames(expr)
  batch = unlist(sapply(samples, function(s) 
    ifelse(s %in% batches$sample, batches$batch[batches$sample == s], "")))
  batch = as.factor(batch)
  
  # run ComBat 
  message(".. adjusting for batch effects using ComBat ...")
  mod = model.matrix(~ 1, data = pData(eset))
  expr = ComBat(dat = expr, batch = batch, mod = mod, par.prior = T, 
                prior.plots = F)
  
  # filter probes present in < 20% of the data
  message(".. filtering low-variance probes ...")
  f1 = function(x) (length(which(x == "P")) > 0.8 * length(samples))
  selected = genefilter(calls, f1)
  expr = expr[selected, ]
  
  # map probe IDs to Ensembl 
  message(".. mapping probes to Ensembl accessions ...")
  db = dbs[[array]]
  map = select(db, rownames(expr), columns = c("ENSEMBL"))
  map = map[match(unique(map$PROBEID), map$PROBEID), ]
  
  # print some messages about Ensembl mapping 
  nProbes = nrow(expr)
  message(".... mapping " , nProbes, " probes to Ensembl ...")
  nMapped = nProbes - sum(is.na(map$ENSEMBL))
  message(".... ", nMapped, " probes mapped to a Ensembl accession")
  
  # collapse probesets to Ensembl IDs by median
  ensembl.list = map$ENSEMBL[which(map$PROBEID == row.names(expr))]
  expr.ENSEMBL = aggregate(expr, by = list(ensembl.list), FUN=median)
  
  # read orthologs
  rat.map = read.delim(
    "data/expression/arrayexpress/orthologs/A-AFFY-43-ensembl-orthologs.tsv", 
    strip.white = T)
  mouse.map = read.delim(
    "data/expression/arrayexpress/orthologs/A-AFFY-45-ensembl-orthologs.tsv", 
    strip.white = T)
  if (array == "A-AFFY-43" | array == "A-AFFY-45") {
    if (array == "A-AFFY-43") {
      ortholog.map = rat.map
    } else if (array == "A-AFFY-45") {
      ortholog.map = mouse.map
    }
    # read GTEx data
    gtex = read.csv(
      "data/expression/gtex/expression/filtered/spinal_cord_expr_filtered.csv.gz")
    # map orthologs
    ortholog.map = ortholog.map[which(ortholog.map$MajorityOrtholog != ""), ]
    ortholog.map$MajorityOrtholog = unlist(
      sapply(ortholog.map$MajorityOrtholog,function(x) {
        accessions = unlist(strsplit(x, " "))
        if (length(accessions) > 1) {
          y = accessions[which(accessions %in% gsub("\\..*$","", gtex$Name))]
          if (length(y) == 1) 
            return(y)
          if (length(y) > 1 | length(y) == "0") 
            accessions = accessions[order(accessions)][1]
          return(accessions)
        } else if (length(accessions) == 1)  {
          return(accessions[1])
        } else if (length(accessions) == 0) {
          return("")
        }
      }))
    for (i in 1:length(row.names(expr.ENSEMBL))) {
      name = row.names(expr.ENSEMBL)[i]
      ortholog = ortholog.map$MajorityOrtholog[which(ortholog.map$ENSEMBL == name)]
      if (length(ortholog) == 0) {
        ortholog = "not-mapped"
      }
      row.names(expr.ENSEMBL)[i] = ortholog  
    }
    # get final dataset with non-mapped probes removed
    expr.FINAL = expr.ENSEMBL[-which(row.names(expr.ENSEMBL) == "not-mapped"), ]
        
    # collapse probesets to ENSEMBL IDs by median
    expr.FINAL = aggregate(expr.FINAL, by = row.names(expr.FINAL), FUN = median)
    
    # write
    df = data.frame(expr.FINAL, stringsAsFactors = F)
    write.csv(df, file = paste0(
      "data/expression/arrayexpress/expression/", array, "_expr.csv"), 
      row.names = T, quote = F)
  }
}
