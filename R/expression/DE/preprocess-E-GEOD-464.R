# Use the XPS package to normalize and filter transcriptomic data from a mouse
# model of severity-dependent gene expression. 
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)
library(xps)

# read targets 
expt = "E-GEOD-464"
array = "A-AFFY-18-20"
celDir = file.path("data/expression/arrayexpress/cel", array, expt)
rootFileDir = file.path("data/expression/DE", expt)
targetsFile = paste("data/expression/DE", expt, "targets.tsv", sep = "/")
targets = read.delim(targetsFile)

# set up xps schemes
schemeNames = c("RG_U34A", "RG_U34B", "RG_U34C")
schemes = paste0("data/expression/xps/schemes/", schemeNames, ".root")

# normalize each array separately
for (i in 18:20) {
  # subset to samples from this array 
  a = paste("A-AFFY-", i, sep="")
  subset = targets[targets$array == a,]
  # read scheme
  scheme = root.scheme(schemes[i - 17])  
  # read raw data
  rootFile = paste0(a, ".root")
  data = import.data(scheme, filename = rootFile, filedir = rootFileDir, 
                      celdir = celDir, celfiles = subset$filename)
  # normalize with mas5 
  mas5File = paste0(a, "_mas5.root")
  data.mas5 = mas5(data, filename = mas5File, filedir = rootFileDir, 
                   normalize = T, sc = 500, update = T)
  # extract expression
  # expr.mas5 = validData(data.mas5)
  # calculate calls
  callFile = paste0(a, "_mas5-calls.root")
  call.mas5 = mas5.call(data, filename = callFile, filedir = rootFileDir)
}
