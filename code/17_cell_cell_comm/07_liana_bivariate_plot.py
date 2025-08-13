#   Summarize statistics for ligand-receptor pairs across donors and plot
#   interesting pairs (based on high cosine similarity or high Moran's R)

import pandas as pd
import scanpy as sc
from pyhere import here
import os
import re
import session_info
import matplotlib.pyplot as plt

in_dir = here(
    'processed-data', '17_cell_cell_comm', 'liana', 'bivariate_adatas'
)
plot_dir = here('plots', '17_cell_cell_comm', 'liana')
min_num_donors = 20

#   Read input files as they exist (for one donor, no file was generated)
in_files = [
    os.path.join(in_dir, f)
    for f in os.listdir(in_dir)
    if re.compile(r'.*\.h5ad$').match(f)
]

#   Read in DataFrames of ligand-receptor stats for each donor and concatenate
lr_df_list = []
for f in in_files:
    donor_id = f.split('/')[-1].replace('.h5ad', '')
    lr_df = sc.read(f).var
    lr_df['donor_id'] = donor_id
    lr_df_list.append(lr_df)

lr_df = pd.concat(lr_df_list, axis=0, ignore_index=True)

#   Require a pair to be present in some minimum number of donors
lr_df = lr_df[
    lr_df.groupby(['ligand', 'receptor'])['ligand'].transform('count') >= min_num_donors
]

#   Average stats across samples
lr_mean_df = (
    lr_df
        .groupby(['ligand', 'receptor'], as_index=False)
        .agg({'mean': 'mean', 'morans': 'mean'})
)

lr_top_df = pd.concat(
    [
        lr_mean_df.sort_values("morans", ascending=False).head(),
        lr_mean_df.sort_values("mean", ascending=False).head()
    ]
)
lr_top_df = lr_top_df.drop_duplicates(subset=['ligand', 'receptor'])

top_pairs = list(lr_top_df['ligand'] + '^' + lr_top_df['receptor'])

#   Here we actually just got lucky-- choose the first sample for plotting,
#   since it happens to contain all top pairs
ad_lr = sc.read(in_files[0])
assert all([x in ad_lr.var.index for x in top_pairs])

#   Fix spatial coordinates format
ad_lr.obsm['spatial'] = ad_lr.obsm['spatial'].to_numpy()

#   Plot scores, permutation-based p-values, and local categories
sc.pl.spatial(
    ad_lr, color=top_pairs, cmap='viridis', ncols=2,
    spot_size=200.0
)
plt.savefig(os.path.join(plot_dir, 'bivariate_scores.pdf'))
plt.close('all')

sc.pl.spatial(
    ad_lr, color=top_pairs, cmap='viridis_r', ncols=2,
    spot_size=200.0, layer='pvals'
)
plt.savefig(os.path.join(plot_dir, 'bivariate_p_vals.pdf'))
plt.close('all')

sc.pl.spatial(
    ad_lr, color=top_pairs, cmap='coolwarm', ncols=2,
    spot_size=200.0, layer='cats'
)
plt.savefig(os.path.join(plot_dir, 'bivariate_cats.pdf'))
plt.close('all')

session_info.show()
