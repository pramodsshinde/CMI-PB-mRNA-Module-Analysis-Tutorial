import os 
import sys 
import pandas as pd 
import pypts 

days = [0, 1, 3, 7, 14]

print("# Loading subject-sample map")
master = pypts.load_subject_sample_map()
master = master[master.planned_days_relative_to_boost < 15]
ss_table = master[['subject_id', 'sample_id', 'planned_days_relative_to_boost']]

## Loading Ab titer data 
prediction_day = 14

# loading ab titer data 
print("# loading ab titer data")
ab_titers = pypts.load_api_data('https://staging.cmi-pb.org:443/db/ab_titer')
ab_titers.isotype = ab_titers.isotype.str.strip()
ab_titers = ab_titers[ab_titers.unit == 'IU/ML']
ab_titers = ab_titers.merge(ss_table, on='sample_id')
ab_data_day14 = ab_titers[(ab_titers.isotype == 'IgG') & \
                          (ab_titers.planned_days_relative_to_boost == prediction_day) & \
                         (ab_titers.antigen.isin(['PT', 'FHA', 'PRN']))]

ab_data_day14.loc[:, 'isotype_antigen'] =  ab_data_day14['isotype'] + '_' + ab_data_day14['antigen']
ab_data_day14 = ab_data_day14.pivot(index='sample_id', columns='isotype_antigen', values='ab_titer')

# adding subject and planned days relative to boost 
ab_data_day14 = pd.merge(ss_table, ab_data_day14, left_on='sample_id', right_index=True) 

print("# save the ab titer data")
fn = os.path.join(sys.argv[1])
ab_data_day14.to_csv(fn, index=None, sep='\t')
