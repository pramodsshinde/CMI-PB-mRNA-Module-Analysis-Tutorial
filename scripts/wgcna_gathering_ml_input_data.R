library(stringr)
library(plyr)
options(stringsAsFactors=FALSE); 
argv = commandArgs(trailingOnly=TRUE)

topGenes = list()
for (day in c(0, 1, 3, 7)){

  cat(day, "\n")

  # Setting the path to the data files 
  ge_rdata = sprintf('output/wgcna/runs/wgcna_day_%s.filteredExpression.RData', day)
  net_rdata = sprintf('output/wgcna/runs/wgcna_day_%s.networkConstruction-auto.RData', day)
  moduleTraitCor_fn = sprintf('output/wgcna/runs/wgcna_day_%s.module_vs_igg_corr.tsv', day)
  moduleTraitPvalue_fn = sprintf('output/wgcna/runs/wgcna_day_%s.module_vs_igg_pvals.tsv', day)

  ### Locating modules which correlated significantly with each ab titer 
  sigModules = read.table(moduleTraitPvalue_fn)
  sigModules = sigModules < 0.05

  # storing each of the top modules with in topModules 
  topModules = list()
  abExperiments = colnames(sigModules)
  modules = rownames(sigModules)
  for (ab in abExperiments){
    currTop = modules[sigModules[, ab]]
    topModules[[ab]] = currTop
  }

  ### Extracting the genes belonging to each of the top modules 
  ge_names = load(ge_rdata)
  net_names = load(net_rdata)

  moduleInfo = cbind(moduleLabels, moduleColors)
  moduleInfo = as.data.frame(cbind(moduleLabels, moduleColors))

  # Collecting the genes belonging to the top modules 
  for (ab in abExperiments){
    currTop = gsub("ME", "", topModules[[ab]])
    bools = moduleInfo$moduleColors %in% currTop
    currGenes = moduleInfo[bools, ]
    currGenes = rownames(currGenes)

    key = sprintf("%s.%s", ab, day)
    topGenes[[key]] = datExpr[, currGenes]
  }
}

### concating all IgG-specific data together
### (i.e) cbind(IgG day0, day1, day3, day7)
final_datasets = list()
for (ab in abExperiments){

  cat(ab, "\n")

  for (day in c(0, 1, 3, 7)){

    cat("\t", day, "\n")

    key = sprintf("%s.%s", ab, day)

    if (ncol(topGenes[[key]]) == 0){
      next
    }
    else if (day == 0){
      tdf = topGenes[[key]]
      colnames(tdf) = paste0(colnames(tdf), sprintf('.day%s', day))
      final_datasets[[ab]] = tdf
    }
    else {
      tdf = topGenes[[key]]
      colnames(tdf) = paste0(colnames(tdf), sprintf('.day%s', day))
      final_datasets[[ab]] = merge(final_datasets[[ab]], tdf, 
                                 by='row.names')
      rownames(final_datasets[[ab]]) = final_datasets[[ab]]$Row.names
      final_datasets[[ab]] = subset(final_datasets[[ab]], select=-c(Row.names))
    }
  }
}

# # count the presence of a day 
# days = gsub('ENSG.*\\.[0-9]*\\.', '', colnames(final_datasets[["IgG_FHA"]]))
# count(days)

### Save the final datasets as input to ML algorithms 
for (ab in abExperiments){
  fn = "output/wgcna/runs/wgcna_gene_expression_members_of_sigmod_%s.tsv"
  fn = sprintf(fn, tolower(ab))
  write.table(final_datasets[[ab]], file = fn, quote = FALSE, sep = "\t", 
            row.names = TRUE, col.names = TRUE) 
}



