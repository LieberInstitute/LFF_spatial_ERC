#   Run LIANA by sample on the fine snRNA-seq data, one of two steps largely
#   based on:
#   https://liana-py.readthedocs.io/en/latest/notebooks/mofatalk.html

import scanpy as sc
import liana as li
from pyhere import here
import session_info
import os

ad_in_path = here('processed-data', '17_cell_cell_comm', 'sce.h5ad')
ad_out_path = here(
    'processed-data', '17_cell_cell_comm', 'liana', 'ad_sc_ranked.h5ad'
)

os.makedirs(ad_out_path.parent, exist_ok=True)
num_cores = int(os.getenv('SLURM_CPUS_ON_NODE'))

#   Read in and drop the incorrectly labeled donor
ad_sc = sc.read_h5ad(ad_in_path)
ad_sc = ad_sc[ad_sc.obs['sample_id'] != 'Br1289', :]

#   Must use gene symbols for LIANA+
ad_sc.var.index = ad_sc.var['gene_name'].astype(str)
ad_sc.var.index.name = None

li.mt.rank_aggregate.by_sample(
    ad_sc,
    groupby='cell_type_anno',
    resource_name='consensus',
    sample_key='sample_id',
    expr_prop = 0.1,
    use_raw=False,
    n_perms=1000,
    return_all_lrs=False,
    verbose=True,
    n_jobs=num_cores
)

sc.write(ad_out_path, ad_sc)

session_info.show()
