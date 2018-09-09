# Convert salmon output to sleuth-compatible files using wasabi.
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)
library(wasabi)

# define salmon dir
salmonDir = "data/expression/rnaseq/salmon"

# get salmon files to convert to sleuth compatable versions
files = list.dirs(salmonDir, recursive = F)

for (file in files) {
  # prepare salmon for sleuth
  prepare_fish_for_sleuth(file)
}
  