import anndata as ad
import scanpy as sc
from pyhere import here
import os
import session_info
import scipy.sparse as sp
import numpy as np

donor = 'Br1039'
interaction = 'Astro.5^SPON1^APP'
in_path = here(
    'processed-data', '24_xenium_liana', '03_inflow_score',
    f'lrdata_{donor}.h5ad'
)
plot_dir = here('plots', '24_xenium_liana', '06_spatial_plots')

os.makedirs(plot_dir, exist_ok=True)

lr_data = ad.read_h5ad(in_path)

#   Cap at 99th percentile to improve color scale
vals = lr_data[:, interaction].X
if sp.issparse(vals):
    vals = vals.toarray().flatten()
else:
    vals = vals.flatten()
pct_99 = np.quantile(vals, 0.99)

p = sc.pl.embedding(
    lr_data, basis='spatial', color=interaction, vmax=pct_99, s=10,
    cmap='magma', frameon=False, return_fig=True
)
p.savefig(plot_dir / f"inflow_{donor}.pdf", dpi=300, bbox_inches='tight')

session_info.show()
