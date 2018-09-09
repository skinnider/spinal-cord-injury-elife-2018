# functions used to analyze differential gene expression

enr <- function(seeds, modules, correct=p.adjust.methods) {
  # Count modules
  nMods <- length(unique(modules$module))
  # Calculate enrichment
  module.df <- data.frame(module = seq_len(nMods-1), enrich.p = 0)
  for (i in seq_len(nMods - 1)) {
    moduleGenes <- modules$gene[modules$module == i]
    otherGenes <- modules$gene[modules$module != i]
    in.module <- length(intersect(seeds, moduleGenes))
    not.in.module <- length(moduleGenes) - in.module
    in.other <- length(intersect(seeds, otherGenes))
    not.in.other <- length(otherGenes) - in.other
    mat <- matrix(c(in.module, in.other, not.in.module, not.in.other), nrow=2, 
                  dimnames=list(c("Module", "Outside"), c("SCI", "Non-SCI")))
    f <- fisher.test(mat, alternative="greater")
    module.df[i,2] <- f$p.value
  }
  # Apply MHT correction
  module.df$enrich.p <- p.adjust(module.df$enrich.p, method=correct)
  return(module.df)
}

pick.orthologs <- function(orthologs, targets) {
  orthologs[,2] <- unlist(sapply(orthologs[,2], function(x) {
    accessions <- unlist(strsplit(x, " "))
    if (length(accessions) > 1) {
      # Try to match to the target genes first 
      y <- accessions[which(accessions %in% targets)]
      if (length(y) == 1) {
        return(y)
      } else {
        return(accessions[order(accessions)][1]) # alphabetize 
      }
      return(accessions)
    } else if (length(accessions) == 1)  {
      return(accessions[1])
    } else if (length(accessions) == 0) {
      return("")
    }
  }))
  return(orthologs)
}
