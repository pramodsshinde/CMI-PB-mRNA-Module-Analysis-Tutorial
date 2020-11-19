import sys
import pandas as pd
import scipy.stats as ss
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression, ElasticNet
from sklearn.svm import LinearSVR
from sklearn.model_selection import LeavePOut
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('whitegrid')

# Loading a scaler to center each column in the expression
# and ab titers table
scaler = StandardScaler()

# Loading a leave 5 out data splitter for the cross validation step
lpo = LeavePOut(5)

# Gathering command line input
top_genes = sys.argv[1]
ab_data = sys.argv[2]
summary_fn = sys.argv[3]
igg_assay = sys.argv[4]
trials = int(sys.argv[5])

# Internal code to set test files
test = True
test = False
if test is True:
    top_genes = 'output/wgcna/runs/' + \
                'wgcna_gene_expression_members_of_sigmod_igg_fha.tsv'
    ab_data = 'output/wgcna/phenotypes/ab_titers_as_clinical_phenotypes.tsv'
    summary_fn = 'output/wgcna/runs/sklearn_model_summary.png'
    igg_assay = 'IgG-FHA'
    trials = 10

# create the igg assay lable which doesn't use -
igg_assay_label = igg_assay.replace('-', '_')

# Loading and processing  IgG specific data
expr = pd.read_table(top_genes)
rows, cols = expr.index.tolist(), expr.columns.tolist()
expr = scaler.fit_transform(expr)
expr = pd.DataFrame(expr, index=rows, columns=cols)

# Loading and processing Ab titer data
ab_titers = pd.read_table(ab_data)
ab_titers.drop(['sample_id', 'planned_days_relative_to_boost'],
               axis=1, inplace=True)
ab_titers.set_index('subject_id', inplace=True)
rows, cols = ab_titers.index.tolist(), ab_titers.columns.tolist()
ab_titers = scaler.fit_transform(ab_titers)
ab_titers = pd.DataFrame(ab_titers, index=rows, columns=cols)
ab_titers = ab_titers.loc[expr.index.tolist(), igg_assay_label]

# Determining the 5 leave out split
indices = lpo.split(expr)

# Setting a summary list to save the results from all models
summary = []

# ## Testing a least squares
ls_model = LinearRegression()
results = []

for i, (train_idxs, test_idxs) in enumerate(indices):

    train_expr = expr.iloc[train_idxs]
    train_abt = ab_titers.iloc[train_idxs]

    test_expr = expr.iloc[test_idxs]
    test_abt = ab_titers.iloc[test_idxs]

    ls_model.fit(train_expr, train_abt)

    pred_abts = ls_model.predict(test_expr)
    spearman_cor = ss.spearmanr(test_abt, pred_abts)

    results.append([i, spearman_cor.correlation, spearman_cor.pvalue])

    if i == trials:
        break

res = pd.DataFrame(results, columns=['iter', 'spearmanr', 'pvalue'])
res['model'] = 'Least Squares'
summary.append(res)

# ## Testing an ElasticNet model
en_model = ElasticNet()
results = []
for i, (train_idxs, test_idxs) in enumerate(indices):

    train_expr = expr.iloc[train_idxs]
    train_abt = ab_titers.iloc[train_idxs]

    test_expr = expr.iloc[test_idxs]
    test_abt = ab_titers.iloc[test_idxs]

    en_model.fit(train_expr, train_abt)

    pred_abts = en_model.predict(test_expr)
    spearman_cor = ss.spearmanr(test_abt, pred_abts)

    results.append([i, spearman_cor.correlation, spearman_cor.pvalue])

    if i == trials:
        break

res = pd.DataFrame(results, columns=['iter', 'spearmanr', 'pvalue'])
res['model'] = 'Elastic Net'
summary.append(res)

# ## Testing a Lasso model
lsvr_model = LinearSVR()
results = []
for i, (train_idxs, test_idxs) in enumerate(indices):

    train_expr = expr.iloc[train_idxs]
    train_abt = ab_titers.iloc[train_idxs]

    test_expr = expr.iloc[test_idxs]
    test_abt = ab_titers.iloc[test_idxs]

    lsvr_model.fit(train_expr, train_abt)

    pred_abts = lsvr_model.predict(test_expr)
    spearman_cor = ss.spearmanr(test_abt, pred_abts)

    results.append([i, spearman_cor.correlation, spearman_cor.pvalue])

    if i == trials:
        break

res = pd.DataFrame(results, columns=['iter', 'spearmanr', 'pvalue'])
res['model'] = 'LinearSVR'
summary.append(res)

# Concat all the results into a single dataframe
summary = pd.concat(summary)

# Plot the distribution of spearmen cc's across the models
fig, ax = plt.subplots()
suptitle = 'Spearman correlation for {} models'.format(igg_assay)
suptitle = fig.suptitle(suptitle, fontsize=16, y=1.01)
sns.boxplot(x='model', y='spearmanr', data=summary, ax=ax)
n = ab_titers.shape[0]
title = '(test v. preds; n={}; leave 5 out; {} trials)'.format(n, trials)
ax.set_title(title, fontsize=12, pad=10)
ax.set_xlabel('Model', labelpad=10)
ax.set_ylabel('Spearman CC')
fig.savefig(summary_fn, bbox_inches='tight',
            bbox_extra_artists=[suptitle], dpi=300)
