#   Run LIANA by sample on the given data: either broad/fine snRNA-seq data or
#   Visiums data. This is one of two steps largely based on:
#   https://liana-py.readthedocs.io/en/latest/notebooks/mofatalk.html

import scanpy as sc
import liana as li
from pyhere import here
import session_info
import os

task_id = int(os.getenv('SLURM_ARRAY_TASK_ID'))

if task_id == 1:
    ad_in_path = here('processed-data', '17_cell_cell_comm', 'sce.h5ad')
    ad_out_path = here(
        'processed-data', '17_cell_cell_comm', 'liana', 'ranked_adatas',
        'broad.h5ad'
    )
    group_var = 'cell_type_broad'
elif task_id == 2:
    ad_in_path = here('processed-data', '17_cell_cell_comm', 'sce.h5ad')
    ad_out_path = here(
        'processed-data', '17_cell_cell_comm', 'liana', 'ranked_adatas',
        'fine.h5ad'
    )
    group_var = 'cell_type_anno'
else:
    ad_in_path = here('processed-data', '17_cell_cell_comm', 'spe.h5ad')
    ad_out_path = here(
        'processed-data', '17_cell_cell_comm', 'liana', 'ranked_adatas',
        'visium.h5ad'
    )
    group_var = 'SpD'

os.makedirs(ad_out_path.parent, exist_ok=True)
num_cores = int(os.getenv('SLURM_CPUS_ON_NODE'))

#   Read in and drop the incorrectly labeled donor
adata = sc.read_h5ad(ad_in_path)
adata = adata[adata.obs['sample_id'] != 'Br1289', :]

#   Must use gene symbols for LIANA+
adata.var.index = adata.var['gene_name'].astype(str)
adata.var.index.name = None

li.mt.rank_aggregate.by_sample(
    adata,
    groupby=group_var,
    resource_name='consensus',
    sample_key='sample_id',
    expr_prop = 0.1,
    use_raw=False,
    n_perms=1000,
    return_all_lrs=False,
    verbose=True,
    n_jobs=num_cores
)

sc.write(ad_out_path, adata)

session_info.show()
