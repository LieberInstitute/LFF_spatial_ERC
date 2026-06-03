import re
import anndata as ad
import pandas as pd
import plotnine as p9
from pyhere import here
import os
import session_info

data_dir = here('processed-data', '24_xenium_liana', '03_inflow_score')
plot_dir = here('plots', '24_xenium_liana', '04_global_score_heatmap')
out_summary_path = here(
    'processed-data', '24_xenium_liana', '04_global_score_heatmap',
    'global_interactions_summary.csv'
)
out_unique_path = here(
    'processed-data', '24_xenium_liana', '04_global_score_heatmap',
    'unique_ligand_receptor_pairs.csv'
)
n_samples_sig = 4

os.makedirs(plot_dir, exist_ok=True)
os.makedirs(out_summary_path.parent, exist_ok=True)

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
unique_combos = (
    global_interactions[['ligand_complex', 'receptor_complex']]
    .drop_duplicates()
)
print(f'Only {unique_combos.shape[0]} unique ligand-receptor pairs being assessed, ignoring cell type.')
print(f'"APOE" is present as a ligand in {(unique_combos['ligand_complex'] == 'APOE').sum()} unique interaction(s).')

unique_combos.to_csv(out_unique_path, index=False)

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

#   Compute lr_specificity: lr_mean in avg_lr divided by the max lr_mean for
#   that LR pair across all significant source-target combos
lr_max_per_pair = (
    avg_lr
    .groupby(['ligand_complex', 'receptor_complex'], as_index=False)['lr_mean']
    .max()
    .rename(columns={'lr_mean': 'max_lr_mean'})
)
avg_lr = (
    avg_lr
    .merge(lr_max_per_pair, on=['ligand_complex', 'receptor_complex'])
    .assign(lr_specificity=lambda x: x['lr_mean'] / x['max_lr_mean'])
    .drop(columns='max_lr_mean')
)

avg_lr.to_csv(out_summary_path, index=False)

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

#################################################################################
#   SPON1->APP heatmap
#################################################################################

spon1_app_df = (
    avg_lr
    .query("ligand_complex == 'SPON1' and receptor_complex == 'APP'")
    .merge(full_grid, on=['source', 'target'], how='right')
    .fillna({'lr_mean': 0, 'ligand_complex': 'SPON1', 'receptor_complex': 'APP'})
)

p_spon1 = (
    p9.ggplot(spon1_app_df, p9.aes(x='target', y='source', fill='lr_mean'))
    + p9.geom_tile()
    + p9.scale_fill_gradientn(
        colors=['lightgray', '#440154', '#31688e', '#35b779', '#fde725'],
        values=[0, 1e-6, 0.33, 0.66, 1],
        name='Mean\nlr_mean'
    )
    + p9.labs(
        x='Target',
        y='Source',
        title=f'SPON1→APP: mean lr_mean\n(significant in ≥{n_samples_sig} samples)'
    )
    + p9.theme_bw()
    + p9.theme(
        axis_text_x=p9.element_text(rotation=90, ha='center'),
        figure_size=(10, 8)
    )
)

p_spon1.save(
    filename=str(plot_dir / f'spon1_app_heatmap_sig{n_samples_sig}.pdf'),
    verbose=False
)

session_info.show()
