import sys
import pandas as pd
import scipy.stats as ss
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression, ElasticNet, Lasso
from sklearn.model_selection import LeavePOut
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('whitegrid')

# Loading a scaler to center each column in the expression
# and ab titers table
scaler = StandardScaler()

# Loading a leave 5 out data splitter for the cross validation step
lpo = LeavePOut(5)

# Setting the number of trials
trials = 100

# Gathering command line input
top_genes = sys.argv[1]
ab_data = sys.argv[2]
ls_fn = sys.argv[3]
en_fn = sys.argv[4]
lasso_fn = sys.argv[5]
summary_fn = sys.argv[6]

# Internal code to set test files
test = True
if test is True:
    top_genes = 'output/wgcna/runs/' + \
                'wgcna_gene_expression_members_of_sigmod_igg_fha.tsv'
    ab_data = 'output/wgcna/phenotypes/ab_titers_as_clinical_phenotypes.tsv'
    ls_fn = 'output/wgcna/runs/sklearn_least_squares_model.png'
    en_fn = 'output/wgcna/runs/sklearn_elasticnet_model.png'
    lasso_fn = 'output/wgcna/runs/sklearn_lasso_model.png'
    summary_fn = 'output/wgcna/runs/sklearn_model_summary.png'


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

# Determining the 5 leave out split
indices = lpo.split(expr)

# Setting a summary list to save the results from all models
summary = []

# ## Testing a least squares
ls_model = LinearRegression()
results = []
for ab in ab_titers.columns.tolist():

    curr_abts = ab_titers[ab].loc[expr.index.tolist()]

    for i, (train_idxs, test_idxs) in enumerate(indices):

        train_expr = expr.iloc[train_idxs]
        train_abt = curr_abts.iloc[train_idxs]

        test_expr = expr.iloc[test_idxs]
        test_abt = curr_abts.iloc[test_idxs]

        ls_model.fit(train_expr, train_abt)

        pred_abts = ls_model.predict(test_expr)
        spearman_cor = ss.spearmanr(test_abt, pred_abts)

        results.append([ab, i, spearman_cor.correlation, spearman_cor.pvalue])

        if i == trials:
            break

res = pd.DataFrame(results, columns=['ab', 'iter', 'spearmanr', 'pvalue'])
res['model'] = 'Least Squares'
summary.append(res)

# Plot the distribution of spearmen cc's
fig, ax = plt.subplots()
fig.suptitle('Least Squares', fontsize=20, y=1.1)
sns.boxplot(x='ab', y='spearmanr', data=res, ax=ax)
title = 'Spearman correlation of antibody titer assays'
title += '\n(test v. preds; {} trials; leave 5 out)'.format(trials)
ax.set_title(title, pad=10)
ax.set_xlabel('Ab titer assay', labelpad=10)
ax.set_ylabel('Spearman CC')

# Replace _ in x labels
new_xlabs = [x.get_text().replace('_', '-') for x in ax.get_xticklabels()]
ax.set_xticklabels(new_xlabs)
fig.savefig(ls_fn, dpi=300)


### Testing an ElasticNet model
en_model = ElasticNet()
results = []
for ab in ab_titers.columns.tolist():

    curr_abts = ab_titers[ab].loc[expr.index.tolist()]

    for i, (train_idxs, test_idxs) in enumerate(indices):

        train_expr = expr.iloc[train_idxs]
        train_abt = curr_abts.iloc[train_idxs]

        test_expr = expr.iloc[test_idxs]
        test_abt = curr_abts.iloc[test_idxs]

        en_model.fit(train_expr, train_abt)

        pred_abts = en_model.predict(test_expr)
        spearman_cor = ss.spearmanr(test_abt, pred_abts)

        results.append([ab, i, spearman_cor.correlation, spearman_cor.pvalue])

        if i == trials:
            break

res = pd.DataFrame(results,
        columns=['ab', 'iter', 'spearmanr', 'pvalue'])
res['model'] = 'Elastic Net'
summary.append(res)

# Plot the distribution of spearmen cc's
fig, ax = plt.subplots()
fig.suptitle('ElasticNet', fontsize=20, y=1.1)
sns.boxplot(x='ab', y='spearmanr', data=res, ax=ax)
title = 'Spearman correlation of antibody titer assays'
title += '\n(test v. preds; {} trials; leave 5 out)'.format(trials)
ax.set_title(title, pad=10)
ax.set_xlabel('Ab titer assay', labelpad=10)
ax.set_ylabel('Spearman CC')

# Replace _ in x labels
new_xlabs = [x.get_text().replace('_', '-') for x in ax.get_xticklabels()]
ax.set_xticklabels(new_xlabs)
fig.savefig(en_fn, dpi=300)


# ## Testing a Lasso model
lasso_model = Lasso()
results = []
for ab in ab_titers.columns.tolist():

    curr_abts = ab_titers[ab].loc[expr.index.tolist()]

    for i, (train_idxs, test_idxs) in enumerate(indices):

        train_expr = expr.iloc[train_idxs]
        train_abt = curr_abts.iloc[train_idxs]

        test_expr = expr.iloc[test_idxs]
        test_abt = curr_abts.iloc[test_idxs]

        lasso_model.fit(train_expr, train_abt)

        pred_abts = lasso_model.predict(test_expr)
        spearman_cor = ss.spearmanr(test_abt, pred_abts)

        results.append([ab, i, spearman_cor.correlation, spearman_cor.pvalue])

        if i == trials:
            break

res = pd.DataFrame(results, columns=['ab', 'iter', 'spearmanr', 'pvalue'])
res['model'] = 'Lasso'
summary.append(res)

# Plot the distribution of spearmen cc's
fig, ax = plt.subplots()
fig.suptitle('Lasso', fontsize=20, y=1.1)
sns.boxplot(x='ab', y='spearmanr', data=res, ax=ax)
title = 'Spearman correlation of antibody titer assays'
title += '\n(test v. preds; {} trials; leave 5 out)'.format(trials)
ax.set_title(title, pad=10)
ax.set_xlabel('Ab titer assay', labelpad=10)
ax.set_ylabel('Spearman CC')

# Replace _ in x labels
new_xlabs = [x.get_text().replace('_', '-') for x in ax.get_xticklabels()]
ax.set_xticklabels(new_xlabs)
fig.savefig(lasso_fn, dpi=300)


# Plot the distribution of spearmen cc's across the models
summary = pd.concat(summary)

fig, ax = plt.subplots()
# fig.suptitle('Lasso', fontsize=20, y=1.1)
sns.boxplot(x='model', y='spearmanr', data=summary, ax=ax)
title = 'Spearman correlation of antibody titer assays'
title += '\n(test v. preds; {} trials; leave 5 out)'.format(trials)
ax.set_title(title, pad=10)
ax.set_xlabel('Model', labelpad=10)
ax.set_ylabel('Spearman CC')
fig.savefig(summary_fn, dpi=300)
