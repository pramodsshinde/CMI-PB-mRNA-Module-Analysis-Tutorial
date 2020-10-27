library(WGCNA);
library(stringr)
options(stringsAsFactors = FALSE);

# Getting the commandline values 
argv = commandArgs(trailingOnly=TRUE)

# Setting global values from user input 
input_file = argv[1]
dendrogram = argv[2]
dendroConst = argv[3]

test = TRUE
test = FALSE
if (test){
    input_file = '../nb_output/wgcna/input/wgcna_input_rnaseq_day_0.tsv'
    dendrogram = '../nb_output/wgcna/runs/wgcna_day_0.sampleClustering.png'
    dendroConst = '../nb_output/wgcna/runs/wgcna_day_0.dendrogramConstructed.RData'
}

day = str_match(input_file, pattern = 'wgcna_input_rnaseq_day_([:digit:]*)')[1,2]

###### Step 1: Filter out bad genes and samples 

# Read in the input data set
datExpr0 = read.csv(input_file, row.names = 1, sep = '\t');

# Set the row names to subject id 
row.names(datExpr0) = datExpr0$subject_id

# Remove unnecessary columns 
datExpr0 = subset(datExpr0, select=-c(subject_id, planned_days_relative_to_boost))

###### Step 2: Plot a tree of the samples to investigate outliers
# Plot the sample tree
sampleTree = hclust(dist(datExpr0), method = "average");
png(file = dendrogram, width = 12, height = 9, units = 'in', res=300);
par(cex = 0.6);
par(mar = c(5,8,8,3))
title = paste0("Sample clustering to detect outliers: Day ", day)
plot(sampleTree, main = title, sub="", xlab="", 
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()

# Saving the current expression data and sampleTree
save(datExpr0, sampleTree, file = dendroConst)


