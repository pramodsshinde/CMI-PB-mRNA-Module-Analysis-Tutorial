#!/share/apps/R/3.3.3/bin Rscript

library(WGCNA);
library(stringr);
options(stringsAsFactors = FALSE);

# Getting the commandline values 
argv = commandArgs(trailingOnly=TRUE)

# Setting global values from user input 
fltExpr = argv[1]
netFile = argv[2]
modFile = argv[3]
eigenExprFile = argv[4]
modMapFile = argv[5]

test = TRUE
test = FALSE
if (test){
    fltExpr = "../nb_output/wgcna/runs/wgcna_day_0.filteredExpression.RData"
    netFile = "../nb_output/wgcna/runs/wgcna_day_0.networkConstruction-auto.RData"
    modFile = "../nb_output/wgcna/runs/wgcna_day_0.modules.tsv"
    eigenExprFile = "../nb_output/wgcna/runs/wgcna_day_0.eigengenes.expression.tsv"
    modMapFile = "../nb_output/wgcna/runs/wgcna_day_0.modulemap.tsv"
}

workingDir = dirname(netFile)
sample = str_split(basename(netFile), '\\.')[[1]][1]

lnames = load(fltExpr)

###### Step 1: Determine the modules
## Allow multi-threading within WGCNA. This helps speed up certain calculations.
## At present this call is necessary for the code to work.
##enableWGCNAThreads() # needed for Rstudio 
allowWGCNAThreads() # needed for Rscript 

# Determine the modules 
net = blockwiseModules(datExpr, power = 6,
    TOMType = "unsigned", minModuleSize = 30,
    reassignThreshold = 0, mergeCutHeight = 0.25,
    numericLabels = TRUE, pamRespectsDendro = FALSE,
    saveTOMs = TRUE,
    saveTOMFileBase = paste0(workingDir, '/', sample, "-TOM"),
    verbose = 3)

###### Step 2: Plot the dendrogram and the module colors underneath

# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
fn = paste0(workingDir, '/', sample, ".dendogramModules.png")
png(file = fn, width = 12, height = 9, units="in", res=300);
par(mar = c(6, 8.5, 3, 3));
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
    "Module colors", dendroLabels = FALSE, hang = 0.03,
    addGuide = TRUE, guideHang = 0.05)
dev.off()

###### Step 3: Save module information 
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
rownames(MEs) = rownames(datExpr)
geneTree = net$dendrograms[[1]];

# Save the modules info into an RData format (further analysis in R)
save(MEs, moduleLabels, moduleColors, geneTree, file = netFile)

# Save the modules/gene pairs into tsv format (further analysis in Python)
geneMods <- do.call(rbind, Map(data.frame, gene=colnames(datExpr), 
                                        moduleColor=moduleColors))
write.table(geneMods,file=modFile, quote = FALSE, row.names = FALSE, sep = '\t')

# Save the eigengene information into tsv format (further analysis in Python)
write.table(MEs,file=eigenExprFile, quote = FALSE, row.names = TRUE, sep = '\t')

# Save a map of module colors and labels (further analysis in Python)
moduleMap = do.call(rbind, Map(data.frame, color=moduleColors, label=moduleLabels))
moduleMap = moduleMap[!duplicated(moduleMap[c(1,2)]),]
moduleMap = as.data.frame(moduleMap)
write.table(moduleMap,file=modMapFile, quote = FALSE, row.names = FALSE, sep = '\t')










