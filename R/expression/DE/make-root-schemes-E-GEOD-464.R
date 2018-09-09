# Make the root schemes to analyze E-GEOD-464.
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)

xpsDir = "data/expression/xps"
libDir = file.path(xpsDir, "libraryfiles")
annDir = file.path(xpsDir, "Annotation")
schemeDir = file.path(xpsDir, "schemes")

chips = c("RG_U34A", "RG_U34B", "RG_U34C")
for (chip in chips) {
  schemeFilename = paste0(chip, ".cdf")
  probeFilename = paste0(chip, ".probe.tab")
  annotFilename = paste0(chip, ".na36.annot.csv")
  scheme = import.expr.scheme(chip, filedir = schemeDir,
                              schemefile = file.path(libDir, schemeFilename),
                              probefile = file.path(libDir, probeFilename),
                              annotfile = file.path(annDir, annotFilename))
}
