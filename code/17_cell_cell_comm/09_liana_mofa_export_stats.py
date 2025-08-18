#   Really should've been part of 04_liana_mofa_run.py-- export initial LR stats
#   as CSV for analysis outside of python

import scanpy as sc
from pyhere import here
import os
import session_info

task_id = int(os.getenv('SLURM_ARRAY_TASK_ID')) - 1
data_description = ['broad', 'fine', 'visium'][task_id]

ad_in_path = here(
    'processed-data', '17_cell_cell_comm', 'liana', 'ranked_adatas',
    f'{data_description}.h5ad'
)
out_path = here(
    'processed-data', '17_cell_cell_comm', 'liana', 'ranked_results',
    f'{data_description}.csv.gz'
)

os.makedirs(out_path.parent, exist_ok=True)

adata = sc.read_h5ad(ad_in_path)
adata.uns["liana_res"].to_csv(out_path, index=False)

session_info.show()
