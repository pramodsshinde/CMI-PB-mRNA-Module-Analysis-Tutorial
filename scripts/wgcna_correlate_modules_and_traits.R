library(WGCNA)
library(dplyr)
library(stringr)

options(stringsAsFactors=FALSE); 
argv = commandArgs(trailingOnly=TRUE)

ss_table_fn = argv[1]
ge_rdata = argv[2]
net_rdata = argv[3]
datTraits_fn = argv[4]
png_fn = argv[5]
moduleTraitCor_fn = argv[6]
moduleTraitPvalue_fn = argv[7]

test = TRUE
test = FALSE 
if (test) {
  ss_table_fn = '../output/database_dump/subject_sample_days.tsv'
  ge_rdata = '../output/wgcna/runs/wgcna_day_3-filteredExpression.RData'
  net_rdata = '../output/wgcna/runs/wgcna_day_3-networkConstruction-auto.RData'
  datTraits_fn = '../output/wgcna/input/ab_titers_as_clinical_phenotypes.tsv'
  png_fn = '../output/wgcna/runs/wgcna_day_3-module_vg_igg_corr.png'
  moduleTraitCor_fn = '../output/wgcna/runs/wgcna_day_3-module_vg_igg_corr.tsv'
  moduleTraitPvalue_fn = '../output/wgcna/runs/wgcna_day_3-module_vg_igg_pvals.tsv'
}

day = str_match(ge_rdata, pattern = 'wgcna_day_([:digit:]*)')[1,2]

# Loading sample to subject/day
ss_table = read.csv(ss_table_fn, sep='\t')
ss_table = ss_table[(ss_table$planned_days_relative_to_boost == day),]

# Loading RData from run_wgnca.R
ge_names = load(ge_rdata) # includes datExpr 
net_names = load(net_rdata) # includes module info

# Adding meta info to  expression data 
datExpr$subject_id = rownames(datExpr)
datExpr = merge(ss_table, datExpr, by='subject_id')

# Load clinical traits data 
datTraits = read.csv(datTraits_fn, sep="\t")
datTraits = datTraits %>% filter(subject_id %in% datExpr$subject_id)

# Remove extra columns 
drop = c("sample_id", "subject_id", "planned_days_relative_to_boost")
datExpr = datExpr[, !(colnames(datExpr) %in% drop)]

drop = c("sample_id", "subject_id", "planned_days_relative_to_boost")
datTraits = datTraits[, !(colnames(datTraits) %in% drop)]

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

print(MEs)
print(datTraits)

moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

# Visualize the heatmap of correlations 
png(file = png_fn, width = 12, height = 9, units = 'in', res=300);
par(mar = c(6, 8.5, 3, 3));

# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                           signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

# Save the module trait correlation and pvalues 
write.table(x=moduleTraitCor, file=moduleTraitCor_fn, sep="\t", quote = FALSE)
write.table(x=moduleTraitPvalue, file=moduleTraitPvalue_fn, sep='\t', quote = FALSE)



