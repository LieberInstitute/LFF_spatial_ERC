import anndata as ad
import scanpy as sc
from pyhere import here
import os
import session_info
import scipy.sparse as sp
import numpy as np

donor = 'Br1039'
interaction = 'Astro.5^SPON1^APP'
lr_in_path = here(
    'processed-data', '24_xenium_liana', '03_inflow_score',
    f'lrdata_{donor}.h5ad'
)
ad_in_path = here(
    'processed-data', '24_xenium_liana', '01_to_python',
    'cell_rna_converted_gobject.h5ad'
)
plot_dir = here('plots', '24_xenium_liana', '06_spatial_plots')

os.makedirs(plot_dir, exist_ok=True)

lr_data = ad.read_h5ad(lr_in_path)
adata = ad.read_h5ad(ad_in_path)
adata = adata[adata.obs['sample_id'] == donor, :]

#   Cap at 99th percentile to improve color scale
vals = lr_data[:, interaction].X
if sp.issparse(vals):
    vals = vals.toarray().flatten()
else:
    vals = vals.flatten()
pct_99 = np.quantile(vals, 0.99)

#   Inflow score spatial plot
p = sc.pl.embedding(
    lr_data, basis='spatial', color=interaction, vmax=pct_99, s=10,
    cmap='magma', frameon=False, return_fig=True
)
p.savefig(plot_dir / f"inflow_{donor}.pdf", dpi=300, bbox_inches='tight')

#   Ligand and receptor expression individually
comp = interaction.split("^")
color_max = [
    float(np.quantile(adata[:, comp[1]].X.toarray().flatten(), 0.99)),
    float(np.quantile(adata[:, comp[2]].X.toarray().flatten(), 0.99))
]
p = sc.pl.embedding(
    adata, basis='spatial', color=[comp[1], comp[2]], s=10, use_raw=False,
    ncols=2, frameon=False, return_fig=True, vmax=color_max
)
p.savefig(plot_dir / f"LR_{donor}.pdf", dpi=300, bbox_inches='tight')

session_info.show()
