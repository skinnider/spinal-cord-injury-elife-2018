# Quantify gene expression in raw RNA-seq data from SCI severity experiment 
# using Salmon.
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)

# define fastq directory
fastqDir = "data/expression/rnaseq/fastq"
# define Salmon output directory
salmonDir = "data/expression/rnaseq/salmon"
# define Salmon binary 
salmonbin = "~/Salmon-0.8.2_macOX_10.12/bin/salmon"
# define index file
index = "data/expression/rnaseq/index/Rattus_norvegicus.Rnor_6.0.89.cdna.ncrna.0.8.2.sidx"
# define libtype: let Salmon infer
libtype = "A"
# fix salmon path
system("export DYLD_FALLBACK_LIBRARY_PATH=~/Salmon-0.8.2_macOX_10.12/lib")

# define runs to process
runs = list.files("data/expression/rnaseq/fastq", pattern = "*.fastq.gz",
                  full.names = T)
# get only first pair (second to be appended during salmon call)
runs = runs[grepl("R1", runs)]

# run Salmon on each run
for (i in seq_len(length(runs))) {
  run = runs[i]
  message(".. analyzing run ", i, " ...")
  salmonRunDir = file.path(salmonDir, basename(run))
  salmonOutput = file.path(salmonRunDir, "quant.sf")
  if (!file.exists(salmonOutput)) {
    # run Salmon
    file1 = run
    file2 = file.path(fastqDir, paste0(gsub("_R.*", "", basename(run)),
                                       "_R2_001.fastq.gz"))
    salmon = sprintf(
      "bash -c '%s quant --numBootstraps 100 -p 10 -l %s -i %s -1 %s -2 %s -o %s'",
      salmonbin, libtype, index, file1, file2, salmonRunDir)
    system(salmon)
  } else {
    message("Salmon has already been run for ", run)
  }
}
