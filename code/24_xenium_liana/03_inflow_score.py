import scanpy as sc
import liana as li
import os
from pyhere import here
import session_info

task_id = int(os.environ.get('SLURM_ARRAY_TASK_ID', '1')) - 1

ad_path = here(
    'processed-data', '24_xenium_liana', '01_to_python',
    'cell_rna_converted_gobject.h5ad'
)
bandwidth = 32

#   Read in just one sample depending on the array task ID
adata = sc.read_h5ad(ad_path)
sample_id = adata.obs['sample_id'].cat.categories[task_id]
adata = adata[adata.obs['sample_id'] == sample_id]

out_path = here(
    'processed-data', '24_xenium_liana', '03_inflow_score',
    f'lrdata_{sample_id}.h5ad'
)

os.makedirs(out_path.parent, exist_ok=True)

print(f"Processing sample {sample_id} with bandwidth {bandwidth}...")

li.ut.spatial_neighbors(
    adata=adata, bandwidth=bandwidth, spatial_key="spatial"
)

lrdata = li.mt.inflow(
    adata, groupby='cell_type_anno',
    resource=li.rs.select_resource("consensus"), use_raw=False
)

li.mt.compute_global_specificity(
    lrdata, groupby='cell_type_anno', use_raw=False, verbose=False
)

lrdata.write_h5ad(out_path)

session_info.show()
