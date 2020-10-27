"""
Generate the WGCNA input data from the RNA-seq dump file.
sys.argv[1] is the CMI-PB RNA-seq dump file 
"""

import os 
import sys 
import pandas as pd 
import pypts 

rnaseq_dump = sys.argv[1] 
outdir = sys.argv[2]

# ### Making the WGCNA data tables
rnaseq = pd.read_table(rnaseq_dump) # dumped from load_api_data('https://staging.cmi-pb.org:443/db/rnaseq')
rnaseq = rnaseq.pivot(index='sample_id', columns='versioned_ensembl_gene_id', values='tpm')
print('Total number of genes loaded:', rnaseq.shape[1])

# filter inclusion for protein coding genes
gene_info = pypts.load_api_data('https://staging.cmi-pb.org:443/db/gene')
gene_info.set_index('versioned_ensembl_gene_id', inplace=True)
gene_info.drop(['description', 'ncbi_id', 'summary', 'synonyms', 'strand'], axis=1, inplace=True)

protein_bools = gene_info.biotype.str.contains('Protein coding')
gene_info = gene_info[protein_bools]
gene_ids = gene_info.index.unique().tolist()
gene_ids = set(rnaseq.columns.tolist()).intersection(gene_ids)
rnaseq = rnaseq.loc[:, gene_ids]
print('Total number of protein-coding genes:', rnaseq.shape[1])

# ### Merge with master
# loading the master data 
master = pypts.load_subject_sample_map()
master = master[['subject_id', 'sample_id', 'planned_days_relative_to_boost']]

# merge with the master for further filtering 
rnaseq = rnaseq.merge(master, left_index=True, right_on='sample_id')

# filtering whole genes (columns) if the number of individuals with greater than thresh tpm is 5 
def count_gt_thresh(l, t):
    """
    Params:
    l: list
        A list of values for which threshold will be used to count less than. 
    t: int 
        A threshold value 
    """
    
    count = 0
    for x in l: 
        if x > t:
            count += 1 
    return(count)

tpm_thresh = 2 
obs_thresh = 5 
master_cols = ['subject_id', 'sample_id', 'planned_days_relative_to_boost']
gene_thresh_cnts = rnaseq.drop(master_cols, errors='ignore', axis=1).apply(count_gt_thresh, t=tpm_thresh)
viable_genes = gene_thresh_cnts[gene_thresh_cnts >= obs_thresh]
rnaseq = rnaseq.loc[:, viable_genes.index.tolist()]
print('Total number of genes after filtering for those expressed:', rnaseq.shape[1])

# adding the meta data after filtering 
rnaseq = pd.merge(master, rnaseq, left_on='sample_id', right_index=True)#.drop('sample_id', inplace=False, axis=1)

# saving the complete RNA-seq file 
fn = os.path.join(outdir, 'wgcna_rnaseq_matrix.tsv')
rnaseq.to_csv(fn, sep='\t', header=True, index=True)

# ### Splitting the data to day by day 
# setting the days that will be included and analyzed 
days = [0, 1, 3, 7]

# Obtaining an intersect of subject ids from days 0, 1, 3, and 7
grps = rnaseq.groupby('planned_days_relative_to_boost')
day0_data = grps.get_group(days[0])
overlapping_subjects = set(day0_data.subject_id)
for day in days[1:]:
    day_df = grps.get_group(day)
    day_subjects = set(day_df.subject_id)   
    overlapping_subjects.intersection_update(day_subjects)
print('Total number of subject within the RNA-seq data:', len(overlapping_subjects))
for s in overlapping_subjects:
    print(s)

grps = rnaseq.groupby('planned_days_relative_to_boost')
for day in days:
    
    print('Generating WGNCA input data for day:', day)

    day_df = grps.get_group(day)
    day_df = day_df[day_df.subject_id.isin(overlapping_subjects)]

    # moving subject_id and planned_days_relative_to_boost to the leftmost column 
    reordered = ['sample_id', 'subject_id', 'planned_days_relative_to_boost'] 
    reordered += [x for x in day_df.columns.tolist() if x not in reordered]
    day_df = day_df[reordered]

    # saving the rna file for this day 
    fn = os.path.join(outdir, 'wgcna_input_rnaseq_day_{}.tsv'.format(day))
    day_df.to_csv(fn, sep='\t', header=True, index=False)
    
# ## Obtaining a summary of the samples per day 
samples_per_day = rnaseq.planned_days_relative_to_boost.value_counts().to_frame()
fn = os.path.join(outdir, 'final_number_of_samples_per_day_uncensored.tsv')
samples_per_day.to_csv(fn, sep='\t')


