# Make index files for Salmon.
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)

rat_cdna_fa = "data/expression/rnaseq/index/Rattus_norvegicus.Rnor_6.0.cdna.all.fa"
rat_ncrna_fa = "data/expression/rnaseq/index/Rattus_norvegicus.Rnor_6.0.ncrna.fa"
rat_ensembl_version = 89

salmonbin = "~/Salmon-0.8.2_macOX_10.12/bin/salmon"
salmonversion = "0.8.2"

cmd = sprintf("bash -c '%s index -t %s -i %s --type quasi'",
              salmonbin,
              paste0("<(cat ", rat_cdna_fa, " ", rat_ncrna_fa, ")"),
              gsub("cdna.all.fa$", paste0(rat_ensembl_version, 
                                          ".cdna.ncrna.", salmonversion,
                                          ".sidx"), rat_cdna_fa))
message(cmd)
system(cmd)
