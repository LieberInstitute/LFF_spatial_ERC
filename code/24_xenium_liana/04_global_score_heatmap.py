import re
import anndata as ad
import pandas as pd
import plotnine as p9
from pyhere import here
import os
import session_info

data_dir = here('processed-data', '24_xenium_liana', '03_inflow_score')
plot_dir = here('plots', '24_xenium_liana', '04_global_score_heatmap')
out_path = here(
    'processed-data', '24_xenium_liana', '04_global_score_heatmap',
    'global_interactions_summary.csv'
)
n_samples_sig = 4

os.makedirs(plot_dir, exist_ok=True)
os.makedirs(out_path.parent, exist_ok=True)

#################################################################################
#   Gather global scores for all samples
#################################################################################

dfs = []
for path in sorted(data_dir.glob('lrdata_*.h5ad')):
    sample_id = re.search(r'lrdata_(.+)\.h5ad', path.name).group(1)
    adata = ad.read_h5ad(path)
    df = adata.uns['global_interactions'].copy()
    df['sample_id'] = sample_id
    dfs.append(df)

global_interactions = pd.concat(dfs, ignore_index=True)

#   Note that there are only 36 LR pairs being assessed
n_pairs = (
    global_interactions[['ligand_complex', 'receptor_complex']]
    .drop_duplicates()
    .shape[0]
)
print(f'Only {n_pairs} unique ligand-receptor pairs being assessed, ignoring cell type.')

n_samples = global_interactions['sample_id'].nunique()

#################################################################################
#   Aggregate significant interactions across samples
#################################################################################

#   Count how many samples each source-LR-target combo achieves significance
sig_counts = (
    global_interactions[global_interactions['pval'] < 0.05]
    .groupby(['source', 'ligand_complex', 'receptor_complex', 'target'])['sample_id']
    .nunique()
    .reset_index(name='n_sig_samples')
)

#   Keep only combos significant in at least N samples
sig_N = sig_counts[sig_counts['n_sig_samples'] >= n_samples_sig][
    ['source', 'ligand_complex', 'receptor_complex', 'target']
]

#   Average lr_mean across significant samples for each retained combo
avg_lr = (
    global_interactions[global_interactions['pval'] < 0.05]
    .merge(sig_N, on=['source', 'ligand_complex', 'receptor_complex', 'target'])
    .groupby(['source', 'ligand_complex', 'receptor_complex', 'target'], as_index=False)['lr_mean']
    .mean()
)

#   Compute lr_specificity: mean lr_mean per source-LR-target (across all samples)
#   divided by the max mean lr_mean for that LR pair across all source-target combos
lr_means = (
    global_interactions
    .groupby(['source', 'ligand_complex', 'receptor_complex', 'target'], as_index=False)['lr_mean']
    .mean()
    .rename(columns={'lr_mean': 'mean_lr_mean'})
)
lr_max_per_pair = (
    lr_means
    .groupby(['ligand_complex', 'receptor_complex'], as_index=False)['mean_lr_mean']
    .max()
    .rename(columns={'mean_lr_mean': 'max_lr_mean'})
)
lr_specificity = (
    lr_means
    .merge(lr_max_per_pair, on=['ligand_complex', 'receptor_complex'])
    .assign(lr_specificity=lambda x: x['mean_lr_mean'] / x['max_lr_mean'])
    [['source', 'ligand_complex', 'receptor_complex', 'target', 'lr_specificity']]
)
avg_lr = avg_lr.merge(
    lr_specificity, on=['source', 'ligand_complex', 'receptor_complex', 'target']
)

avg_lr.to_csv(out_path, index=False)

#################################################################################
#   Global communication heatmap
#################################################################################

#   Sum averaged lr_mean across LR pairs for each source-target combo,
#   filling missing source-target pairs with 0
cell_types = sorted(global_interactions['source'].unique())
full_grid = pd.MultiIndex.from_product(
    [cell_types, cell_types], names=['source', 'target']
).to_frame(index=False)

heatmap_df = (
    avg_lr
    .groupby(['source', 'target'], as_index=False)['lr_mean']
    .sum()
    .merge(full_grid, on=['source', 'target'], how='right')
    .fillna({'lr_mean': 0})
)

p = (
    p9.ggplot(heatmap_df, p9.aes(x='target', y='source', fill='lr_mean'))
    + p9.geom_tile()
    + p9.scale_fill_gradientn(
        colors=['lightgray', '#440154', '#31688e', '#35b779', '#fde725'],
        values=[0, 1e-6, 0.33, 0.66, 1],
        name='Sum of mean\nlr_mean'
    )
    + p9.labs(
        x='Target',
        y='Source',
        title=f'Sum of mean lr_mean for interactions\nsignificant in ≥{n_samples_sig} samples'
    )
    + p9.theme_bw()
    + p9.theme(
        axis_text_x=p9.element_text(rotation=90, ha='center'),
        figure_size=(10, 8)
    )
)

p.save(
    filename=str(plot_dir / f'source_target_heatmap_sig{n_samples_sig}.pdf'),
    verbose=False
)

session_info.show()
