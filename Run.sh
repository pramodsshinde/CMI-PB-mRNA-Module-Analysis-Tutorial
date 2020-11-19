# Set your R version
PATH=/home/jreyna/software/cmi-pb-pertussis-test/R-4.0.3/bin:$PATH

# Dumping the main tables of the CMI-PB database
#snakemake --cores 1 -p dump_cmi_pb_tables


# Parsing RNA-seq and splitting across multiple days
#snakemake --cores 1 -p output/wgcna/input


# Converting the Ab titers information into WGCNA phenotype info # complete
#snakemake --cores 1 -p output/wgcna/phenotypes/ab_titers_as_clinical_phenotypes.tsv


# Constructing dendrogram data
#snakemake --cores 1 -p output/wgcna/runs/wgcna_day_0.dendrogramConstructed.RData \
#    output/wgcna/runs/wgcna_day_1.dendrogramConstructed.RData \
#    output/wgcna/runs/wgcna_day_3.dendrogramConstructed.RData \
#    output/wgcna/runs/wgcna_day_7.dendrogramConstructed.RData


# Filtering samples from gene expression data
#snakemake --cores 1 -p -f output/wgcna/runs/wgcna_day_0.filteredExpression.RData \
#    output/wgcna/runs/wgcna_day_1.filteredExpression.RData \
#    output/wgcna/runs/wgcna_day_3.filteredExpression.RData \
#    output/wgcna/runs/wgcna_day_7.filteredExpression.RData


# Generating modules
#snakemake --cores 1 -p output/wgcna/runs/wgcna_day_0.networkConstruction-auto.RData \
#    output/wgcna/runs/wgcna_day_1.networkConstruction-auto.RData \
#    output/wgcna/runs/wgcna_day_3.networkConstruction-auto.RData \
#    output/wgcna/runs/wgcna_day_7.networkConstruction-auto.RData


# Correlating the modules with ab titer information
#snakemake --cores 4 -p output/wgcna/runs/wgcna_day_0.module_vs_igg_pvals.tsv \
#        output/wgcna/runs/wgcna_day_1.module_vs_igg_pvals.tsv \
#        output/wgcna/runs/wgcna_day_3.module_vs_igg_pvals.tsv \
#        output/wgcna/runs/wgcna_day_7.module_vs_igg_pvals.tsv


# Gathering all the WGCNA based information for ML input
#snakemake --cores 1 -p wgcna_gathering_ml_input_data


## Generating simple models using SKlearn
snakemake --cores 3 -p output/wgcna/runs/sklearn_model_summary_igg_fha.png \
    output/wgcna/runs/sklearn_model_summary_igg_prn.png \
    output/wgcna/runs/sklearn_model_summary_igg_pt.png
