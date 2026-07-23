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
n_samples_sig = 3

os.makedirs(plot_dir, exist_ok=True)
os.makedirs(out_summary_path.parent, exist_ok=True)

#################################################################################
#   Functions
#################################################################################

def save_faceted_heatmap(df, fill_col, fill_name, title, out_path):
    p = (
        p9.ggplot(df, p9.aes(x='target', y='source', fill=fill_col))
        + p9.geom_tile()
        + p9.scale_fill_gradientn(
            colors=['lightgray', '#440154', '#31688e', '#35b779', '#fde725'],
            values=[0, 1e-6, 0.33, 0.66, 1],
            name=fill_name
        )
        + p9.labs(
            x='Target',
            y='Source',
            title=title
        )
        + p9.facet_wrap('~APOE_carrier')
        + p9.theme_bw()
        + p9.theme(
            axis_text_x=p9.element_text(rotation=90, ha='center'),
            figure_size=(16, 8)
        )
    )

    p.save(filename=str(out_path), verbose=False)


def save_difference_heatmap(df, title, out_path):
    carrier_a, carrier_b = carrier_groups

    diff_df = (
        df
        .pivot(index=['source', 'target'], columns='APOE_carrier', values='lr_mean')
        .reset_index()
        .merge(full_grid_diff, on=['source', 'target'], how='right')
        .fillna(0)
    )
    diff_col = f'{carrier_b}_minus_{carrier_a}'
    diff_df[diff_col] = diff_df[carrier_b] - diff_df[carrier_a]
    max_abs_diff = diff_df[diff_col].abs().max()
    max_abs_diff = max_abs_diff if max_abs_diff > 0 else 1

    p = (
        p9.ggplot(diff_df, p9.aes(x='target', y='source', fill=diff_col))
        + p9.geom_tile()
        + p9.scale_fill_gradient2(
            low='#3b4cc0',
            mid='white',
            high='#b40426',
            midpoint=0,
            limits=(-max_abs_diff, max_abs_diff),
            name=f'{carrier_b} -\n{carrier_a}'
        )
        + p9.labs(
            x='Target',
            y='Source',
            title=title
        )
        + p9.theme_bw()
        + p9.theme(
            axis_text_x=p9.element_text(rotation=90, ha='center'),
            figure_size=(10, 8)
        )
    )

    p.save(filename=str(out_path), verbose=False)

#################################################################################
#   Gather global scores for all samples
#################################################################################

dfs = []
for path in sorted(data_dir.glob('lrdata_*.h5ad')):
    sample_id = re.search(r'lrdata_(.+)\.h5ad', path.name).group(1)
    adata = ad.read_h5ad(path)
    df = adata.uns['global_interactions'].copy()
    df['sample_id'] = sample_id
    df['APOE_carrier'] = adata.obs['APOE_carrier'].iloc[0]
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
#   within each APOE_carrier group
sig_counts = (
    global_interactions[global_interactions['pval'] < 0.05]
    .groupby([
        'source', 'ligand_complex', 'receptor_complex', 'target', 'APOE_carrier'
    ])['sample_id']
    .nunique()
    .reset_index(name='n_sig_samples')
)

#   Keep only combos significant in at least N samples within each carrier group
sig_N = sig_counts[sig_counts['n_sig_samples'] >= n_samples_sig][
    ['source', 'ligand_complex', 'receptor_complex', 'target', 'APOE_carrier']
]

#   Average lr_mean across significant samples for each retained combo within
#   each carrier group
avg_lr = (
    global_interactions[global_interactions['pval'] < 0.05]
    .merge(
        sig_N,
        on=['source', 'ligand_complex', 'receptor_complex', 'target', 'APOE_carrier']
    )
    .groupby([
        'source', 'ligand_complex', 'receptor_complex', 'target', 'APOE_carrier'
    ], as_index=False)['lr_mean']
    .mean()
)

#   Add number of significant samples for each retained combo within each
#   carrier group
avg_lr = avg_lr.merge(
    sig_counts,
    on=['source', 'ligand_complex', 'receptor_complex', 'target', 'APOE_carrier'],
    how='left'
)

#   Compute lr_specificity: lr_mean in avg_lr divided by the max lr_mean for
#   that LR pair across all significant source-target combos within each
#   carrier group
lr_max_per_pair = (
    avg_lr
    .groupby(['ligand_complex', 'receptor_complex', 'APOE_carrier'], as_index=False)['lr_mean']
    .max()
    .rename(columns={'lr_mean': 'max_lr_mean'})
)
avg_lr = (
    avg_lr
    .merge(lr_max_per_pair, on=['ligand_complex', 'receptor_complex', 'APOE_carrier'])
    .assign(lr_specificity=lambda x: x['lr_mean'] / x['max_lr_mean'])
    .drop(columns='max_lr_mean')
)

avg_lr.to_csv(out_summary_path, index=False)

cell_types = sorted(global_interactions['source'].unique())
carrier_groups = sorted(global_interactions['APOE_carrier'].dropna().unique())
full_grid = pd.MultiIndex.from_product(
    [cell_types, cell_types, carrier_groups],
    names=['source', 'target', 'APOE_carrier']
).to_frame(index=False)
full_grid_diff = pd.MultiIndex.from_product(
    [cell_types, cell_types],
    names=['source', 'target']
).to_frame(index=False)

#################################################################################
#   Global communication heatmaps
#################################################################################

#   Sum averaged lr_mean across LR pairs for each source-target combo within
#   each carrier group, filling missing source-target pairs with 0
heatmap_df = (
    avg_lr
    .groupby(['source', 'target', 'APOE_carrier'], as_index=False)['lr_mean']
    .sum()
    .merge(full_grid, on=['source', 'target', 'APOE_carrier'], how='right')
    .fillna({'lr_mean': 0})
)

save_faceted_heatmap(
    heatmap_df,
    fill_col='lr_mean',
    fill_name='Sum of mean\nlr_mean',
    title=f'Sum of mean lr_mean for interactions\nsignificant in ≥{n_samples_sig} samples',
    out_path=plot_dir / f'source_target_heatmap_sig{n_samples_sig}_by_APOE_carrier.pdf'
)

save_difference_heatmap(
    heatmap_df,
    title=f'Difference in summed mean lr_mean\n({carrier_groups[1]} - {carrier_groups[0]})',
    out_path=plot_dir / f'source_target_heatmap_sig{n_samples_sig}_APOE_difference.pdf'
)

#################################################################################
#   SPON1->APP heatmap
#################################################################################

spon1_app_df = (
    avg_lr
    .query("ligand_complex == 'SPON1' and receptor_complex == 'APP'")
    .merge(full_grid, on=['source', 'target', 'APOE_carrier'], how='right')
    .fillna({
        'lr_mean': 0,
        'ligand_complex': 'SPON1',
        'receptor_complex': 'APP'
    })
)

save_faceted_heatmap(
    spon1_app_df,
    fill_col='lr_mean',
    fill_name='Mean\nlr_mean',
    title=f'SPON1→APP: mean lr_mean\n(significant in ≥{n_samples_sig} samples)',
    out_path=plot_dir / f'spon1_app_heatmap_sig{n_samples_sig}_by_APOE_carrier.pdf'
)

save_difference_heatmap(
    spon1_app_df,
    title=f'SPON1→APP mean lr_mean difference\n({carrier_groups[1]} - {carrier_groups[0]})',
    out_path=plot_dir / f'spon1_app_heatmap_sig{n_samples_sig}_APOE_difference.pdf'
)

session_info.show()
