# Run FastQC on raw RNA-seq data from SCI severity experiment.
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)

# read runs 
runs = list.files("data/expression/rnaseq/fastq", pattern = "*.fastq.gz",
                  full.names = T)

# define fastqc binary 
fastqcBin = "~/FastQC/fastqc"

# run FastQC on each library
for (i in seq_len(length(runs))) {
  run = runs[i]
  message(".. analyzing run ", i, " (", which(runs == run), " of ",
          length(runs), ")")
  # get ID and pair
  id = gsub("_R.*", "", basename(run))
  for (pair in c("_R1_001", "_R2_001")) {
    # make directory
    runDir = file.path("data/expression/rnaseq/fastq/fastqc", id)
    if (!dir.exists(runDir)) 
      dir.create(runDir)
    # run fastqc
    if (!file.exists(paste0(runDir, "/", id, pair, "_fastqc.html"))) {
      fastqc = sprintf(
        "bash -c '%s -o %s -f fastq %s'",
        fastqcBin, runDir, run)
      system(fastqc)
    } else {
      message("FastQC has already been run for ", run, pair)
    }
  }
}
