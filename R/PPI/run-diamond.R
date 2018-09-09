# Run the DIAMOnD disease gene recovery algorithm given:
# 1) input network (3 column, A, B, evidence)
# 2) input proteins (single column headless)
# 3) bootstrap ID (integer)
# 4) random or sampled proteins ('random' or 'sampled')
# 5) diamond script location ('/path/to/diamond.py')
# 6) output directory for temporary files and output files
options(stringsAsFactors = F)

# parse arguments
args = commandArgs(trailingOnly = T)
net_in = args[1]
proteins_in = args[2]
i = as.integer(args[3])
type = args[4]
diamond = args[5]
outdir = args[6]

# load in proteins and interactome, subset proteins to those within network
net = read.delim(net_in, col.names = c("A", "B", "evidence"), 
                 comment.char = "#")
proteins = unique(read.delim(proteins_in, header = F)[, 1])
proteins = proteins[proteins %in% c(net$A, net$B)]

# split proteins into 80 / 20
set.seed(i)
train = sample(proteins, (0.8 * length(proteins)))
test = proteins[!proteins %in% train]
if (type == "random") {
  train = sample(unique(c(net$A, net$B)), (0.8 * length(proteins)))
}

# create output directories
tmp = file.path(outdir, "tmp")
if (!dir.exists(tmp)) {
  dir.create(tmp)
}
output = file.path(outdir, "result")
if (!dir.exists(output)) {
  dir.create(output)
}
# define output
basename = paste0(i, "_", type, ".txt")
out_net = file.path(tmp, paste0(basename(net_in), "_", basename))
out_seeds = file.path(tmp, paste0(basename(proteins_in), "_", basename))
out_file = file.path(output, paste0(basename(net_in), "_results_", basename))

# write input files for diamond
write.table(net, out_net, quote = F, row.names = F, sep = "\t")
write.table(train, out_seeds, quote = F, row.names = F, sep = "\t", 
            col.names = F)

# run diamond
command = "python"
arguments = c(diamond, out_net, out_seeds, 1000, out_file)
system2(command, args = arguments)
