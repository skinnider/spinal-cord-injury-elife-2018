getModuleGenes = function(modules) {
  moduleNumbers = unique(modules$module)
  moduleGenes = lapply(moduleNumbers, function(i) 
    modules$gene[modules$module == i])
  df = data.frame(module = moduleNumbers, genes = 
                     sapply(moduleGenes, paste, collapse=" "))
  df = df[order(df$module),]
  return(df)
}

module.enrichment = function(sci.genes, network, exprData, 
                              correct=p.adjust.methods) {
  # Get module genes
  modules = list()
  nMods = length(network$MEs)
  for (i in seq_len(nMods - 1)) 
    modules[[i]] = colnames(exprData)[network$colors == i]
  # Calculate enrichment
  module.df = data.frame(module = seq_len(nMods-1), enrich.p = 0)
  for (i in seq_len(nMods - 1)) {
    module = modules[[i]]
    others = unlist(modules[-i])
    in.module = length(intersect(sci.genes, module))
    not.in.module = length(module) - in.module
    in.other = length(intersect(sci.genes, others))
    not.in.other = length(others) - in.other
    mat = matrix(c(in.module, in.other, not.in.module, not.in.other), nrow=2, 
                  dimnames=list(c("Module", "Outside"), c("SCI", "Non-SCI")))
    f = fisher.test(mat, alternative="greater")
    module.df[i,2] = f$p.value
  }
  # Apply MHT correction
  module.df$enrich.p = p.adjust(module.df$enrich.p, method=correct)
  return(module.df)
}
