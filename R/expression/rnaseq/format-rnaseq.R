# Format gene expression estimates from RNA-seq dataset output by Salmon.
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)
library(org.Rn.eg.db)
library(tximport)

# define salmon dir
salmonDir <- "data/wetlab/salmon"

# get files to write counts for
files <- list.dirs(salmonDir, recursive = F)

# read tx2gene data frame
tx2gene <- read.delim("data/wetlab/tx2gene/tx2gene-rnor-89.txt")


for (file in files) {
  # get quant file
  quant <- file.path(file, "quant.sf")
  # convert to gene-level with tximport
  import <- tximport(quant, type = "salmon", tx2gene = tx2gene, 
                     ignoreTxVersion = T, 
                     countsFromAbundance = "lengthScaledTPM")
  # get id
  id <- gsub(".*/", "", file)
  id <- gsub("\\..*", "", id)
  # append "S" due to numeric beginning
  id <- paste0("S", id)
  # add colname
  colnames(import$counts) <- id
  colnames(import$abundance) <- id
  # write counts
  write.table(import$counts, paste0("data/wetlab/counts/", id, "-counts.txt"),
                                    sep = "\t", quote = F)
  # write tpm
  write.table(import$abundance, paste0("data/wetlab/tpm/", id, "-tpm.txt"),
              sep = "\t", quote = F)
  # gzip
  system(paste0(
    "gzip -f ~/git/spinal-cord-injury-elife-2018/data/wetlab/counts/", id, 
    "-counts.txt"))
  system(paste0(
    "gzip -f ~/git/spinal-cord-injury-elife-2018/data/wetlab/tpm/", id, 
    "-tpm.txt"))
}

# merge all counts 
counts <- list.files("data/wetlab/counts", pattern = ".txt.gz", full.names = T)
list <- lapply(counts, function(x) {
  tmp <- read.delim(x)
  return(tmp)
})

merged <- do.call(cbind, list)

# write and gzip
write.table(merged, "data/wetlab/counts/merged-counts.txt", sep = "\t",
            quote = F)
system(
  "gzip -f ~/git/spinal-cord-injury-elife-2018/data/wetlab/counts/merged-counts.txt")

# merge all tpm 
tpm <- list.files("data/wetlab/tpm", pattern = ".txt.gz", full.names = T)
list <- lapply(tpm, function(x) {
  tmp <- read.delim(x)
  return(tmp)
})

merged <- do.call(cbind, list)

# write and gzip
write.table(merged, "data/wetlab/tpm/merged-tpm.txt", sep = "\t",
            quote = F)
system(
  "gzip -f ~/git/spinal-cord-injury-elife-2018/data/wetlab/tpm/merged-tpm.txt")

