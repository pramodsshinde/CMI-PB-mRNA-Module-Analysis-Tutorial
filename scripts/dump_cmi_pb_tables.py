import os 
import sys
import pandas as pd 
import pypts 

outdir = sys.argv[1]

# dumping a table that maps between samples and subject/planned days  
master = pypts.load_subject_sample_map()
ss_table = master[['subject_id', 'sample_id', 'planned_days_relative_to_boost']]
fn = os.path.join(outdir, 'subject_sample_days.tsv')
ss_table.to_csv(fn, sep='\t', index=None)

# dumping ab titer data 
ab_titers = pypts.load_api_data('https://staging.cmi-pb.org:443/db/ab_titer')
ab_titers = ab_titers[ab_titers.sample_id.isin(ss_table.sample_id.tolist())]
fn = os.path.join(outdir, 'ab_titers_dump.tsv')
ab_titers.to_csv(fn, sep='\t', index=None)

# dumping cytof data
cytof = pypts.load_api_data('https://staging.cmi-pb.org:443/db/cytof')
cytof = cytof[cytof.sample_id.isin(ss_table.sample_id.tolist())]
fn = os.path.join(outdir, 'cytof_dump.tsv')
cytof.to_csv(fn, sep='\t', index=None)

# dumping olink data
olink = pypts.load_api_data('https://staging.cmi-pb.org:443/db/olink_prot_exp')
olink = olink[olink.sample_id.isin(ss_table.sample_id.tolist())]
fn = os.path.join(outdir, 'olink_dump.tsv')
olink.to_csv(fn, sep='\t', index=None)

# dumping rnaseq data
rnaseq = pypts.load_api_data('https://staging.cmi-pb.org:443/db/rnaseq')
fn = os.path.join(outdir, 'rnaseq_dump.tsv')
rnaseq.to_csv(fn, sep='\t')
