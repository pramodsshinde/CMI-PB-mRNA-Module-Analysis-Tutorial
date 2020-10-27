#!/share/apps/R/3.3.3/bin Rscript

library(WGCNA);
options(stringsAsFactors = FALSE);

# Getting the commandline values 
argv = commandArgs(trailingOnly=TRUE)

# Setting global values from user input 
dendroConst = argv[1]
cutoff_txt = argv[2]
fltExpr = argv[3]

test = TRUE
test = FALSE
if (test){
    dendroConst = '../nb_output/wgcna/runs/wgcna_day_0.dendogramConstructed.RData'
    cutoff_txt = '../nb_output/wgcna/manual_inputs/wgcna_day_0.sample_dendrogram_cutoffs.txt'
    fltExpr = '../nb_output/wgcna/runs/wgcna_day_0.filteredExpression.RData'
}

# Setting cut height; will be automatic later 
cutH = read.csv(cutoff_txt, header = FALSE)[1,1]

# Loading expression and sample dendrogram
lnames = load(dendroConst)

# Determine cluster under the line
print("# Determine cluster under the line")
clust = cutreeStatic(sampleTree, cutHeight = cutH, minSize = 5)
table(clust)

# Keeping large clusters
print("# Keeping large clusters")
keepSamples = (clust != 0)
datExpr = datExpr0[keepSamples, ]
save(datExpr, file = fltExpr)
