#   Run MOFA on LIANA+ identified LR pairs, the second of two steps largely
#   based on:
#   https://liana-py.readthedocs.io/en/latest/notebooks/mofatalk.html

import os
import scanpy as sc
import liana as li
import muon as mu
from pyhere import here
import session_info
import mudata
import numpy as np
import pandas as pd

task_id = int(os.getenv('SLURM_ARRAY_TASK_ID')) - 1
data_description = ['broad', 'fine', 'visium'][task_id]

ad_in_path = here(
    'processed-data', '17_cell_cell_comm', 'liana', 'ranked_adatas',
    f'{data_description}.h5ad'
)
mofa_model_path = here(
    'processed-data', '17_cell_cell_comm', 'liana', 'mofa_results',
    f'model_{data_description}.h5ad'
)
mdata_out_path = here(
    'processed-data', '17_cell_cell_comm', 'liana', 'mofa_results',
    f'mdata_{data_description}.h5mu'
)
mofa_meta_vars = ['APOE_carrier', 'Ancestry', 'Sex', 'Age']
mofa_num_factors = 10
samples_per_view_prop = 0.5

if data_description == 'broad':
    group_var = 'cell_type_broad'
elif data_description == 'fine':
    group_var = 'cell_type_anno'
else:
    group_var = 'SpD'

os.makedirs(mofa_model_path.parent, exist_ok=True)

adata = sc.read_h5ad(ad_in_path)

num_samples_per_view = int(
    round(samples_per_view_prop * adata.obs['sample_id'].nunique())
)
mdata = li.multi.lrs_to_views(
    adata,
    sample_key='sample_id',
    score_key='magnitude_rank',
    obs_keys=mofa_meta_vars,
    lr_prop = 0.2, # minimum required proportion of samples to keep an LR
    lrs_per_sample = 20, # min num of interactions to keep a sample in a view
    lrs_per_view = 20, # minimum number of interactions to keep a view
    samples_per_view = num_samples_per_view, # min num samples to keep a view
    min_variance = 0, # minimum variance to keep an interaction
    lr_fill = 0, # fill missing LR values across samples with this
    verbose=True
)

#   Print some info on number of views, samples, and pairs retained
m_df = pd.DataFrame(
    np.array([x.shape for x in dict(mdata.mod).values()]),
    columns = ['num_samples', 'num_pairs']
)

print(f'{m_df.shape[0]} views kept of {adata.obs[group_var].nunique() ** 2}.')
print('Distribution of number of samples by view:')
print(m_df['num_samples'].describe())
print('Distribution of number of pairs by view:')
print(m_df['num_pairs'].describe())

mu.tl.mofa(
    mdata, 
    use_obs='union',
    convergence_mode='medium',
    outfile=mofa_model_path,
    n_factors=mofa_num_factors
)

mudata.write(str(mdata_out_path), mdata)

session_info.show()
