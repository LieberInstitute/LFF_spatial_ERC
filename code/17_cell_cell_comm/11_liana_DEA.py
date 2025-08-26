#   Integrate fine snRNA-seq DEA results with LR statistics, ultimately
#   producing a tileplot. Largely based on:
#   https://liana-py.readthedocs.io/en/latest/notebooks/targeted.html#dea-to-ligand-receptor-interactions

import pandas as pd
import numpy as np
import scanpy as sc
import liana as li
from pyhere import here
import os
import json
import session_info

ad_sc_path = here('processed-data', '17_cell_cell_comm', 'sce.h5ad')
dea_path = here(
    'processed-data', '13_compile_DGE', '01_compile_DGE', 'sn_fine',
    'DGE_results_carrier_sn_fine.csv.gz'
)
plot_dir = here('plots', '17_cell_cell_comm', 'liana')

#   Read in and drop the incorrectly labeled donor
ad_sc = sc.read_h5ad(ad_sc_path)
ad_sc = ad_sc[ad_sc.obs['sample_id'] != 'Br1289', :].copy()

#   Force use of logcounts (assay seems to default to '.X')
ad_sc.X = ad_sc.layers['logcounts']

#   Must use gene symbols for LIANA+
ad_sc.var.index = ad_sc.var['gene_name'].astype(str)

dea_df = (
    pd.read_csv(dea_path, index_col='gene_name')
        .rename(columns={'cluster': 'cell_type_anno'})
)

lr_df = li.multi.df_to_lr(
    ad_sc,
    dea_df=dea_df,
    resource_name='consensus',
    expr_prop=0.1,
    groupby='cell_type_anno',
    stat_keys=['vlmf_t', 'vlmf_P.Value', 'vlmf_adj.P.Val'],
    use_raw=False,
    complex_col='vlmf_t',
    verbose=True,
    return_all_lrs=False
)

p = li.pl.tileplot(
    liana_res=lr_df,
    fill = 'expr',
    label='interaction_vlmf_adj.P.Val',
    label_fun = lambda x: '*' if x < 0.05 else np.nan,
    top_n=15,
    orderby = 'interaction_vlmf_t',
    orderby_ascending = False,
    orderby_absolute = False,
    source_title='Ligand',
    target_title='Receptor',
)
p.save(
    filename=os.path.join(plot_dir, f'dea_tileplot_fine.pdf'),
    width=12, height=7
)

session_info.show()
