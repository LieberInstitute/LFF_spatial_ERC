import numpy as np
import pandas as pd
import scanpy as sc
import plotnine as p9
import liana as li
import muon as mu
import mofax as mofa
import decoupler as dc
from pyhere import here
import matplotlib.pyplot as plt
import os

ad_sc_path = here('processed-data', '17_cell_cell_comm', 'sce.h5ad')
mofa_model_path = here(
    'processed-data', '17_cell_cell_comm', 'liana', 'mofa_model.h5ad'
)
plot_dir = here('plots', '17_cell_cell_comm', 'liana')
mofa_meta_vars = ['APOE_carrier', 'Ancestry', 'Sex', 'Age']
mofa_num_factors = 10
samples_per_view_prop = 0.8

os.makedirs(mofa_model_path.parent, exist_ok=True)
os.makedirs(plot_dir, exist_ok=True)

ad_sc = sc.read_h5ad(ad_sc_path)

#   Must use gene symbols for LIANA+
ad_sc.var.index = ad_sc.var['gene_name'].astype(str)
ad_sc.var.index.name = None

# ad_sc = ad_sc[ad_sc.obs['sample_id'].isin(ad_sc.obs['sample_id'].unique()[:2]), :].copy()

li.mt.rank_aggregate.by_sample(
    ad_sc,
    groupby='cell_type_anno',
    resource_name='consensus',
    sample_key='sample_id',
    expr_prop = 0.1,
    use_raw=False,
    n_perms=100,
    return_all_lrs=False,
    verbose=True
)

num_samples_per_view = int(
    round(samples_per_view_prop * ad_sc.obs['sample_id'].nunique())
)
mdata = li.multi.lrs_to_views(
    ad_sc,
    sample_key='sample_id',
    score_key='magnitude_rank',
    obs_keys=mofa_meta_vars,
    lr_prop = 0.3, # minimum required proportion of samples to keep an LR
    lrs_per_sample = 20, # min num of interactions to keep a sample in a view
    lrs_per_view = 20, # minimum number of interactions to keep a view
    samples_per_view = num_samples_per_view, # min num samples to keep a view
    min_variance = 0, # minimum variance to keep an interaction
    lr_fill = 0, # fill missing LR values across samples with this
    verbose=True
)

mu.tl.mofa(
    mdata, 
    use_obs='union',
    convergence_mode='medium',
    outfile=mofa_model_path,
    n_factors=mofa_num_factors
)

adata = sc.AnnData(obs = mdata.obs, obsm = mdata.obsm)

dc.tl.rankby_obsm(adata, key='X_mofa', uns_key='rank_obsm')
dc.pl.obsm(
    adata,
    key='rank_obsm',
    names = mofa_meta_vars
    titles=['Principle component scores', 'Adjusted p-values from ANOVA'],
    figsize=(7, 5),
    nvar=10
)
plt.savefig(os.path.join(plot_dir, 'mofa_heatmap.pdf'))
plt.close('all')
