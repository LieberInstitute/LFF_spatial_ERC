import pandas as pd
import scanpy as sc
import decoupler as dc
import liana as li
from matplotlib import pyplot as plt
from mudata import MuData
from pyhere import here
import os

ad_sp_path = here('processed-data', '17_cell_cell_comm', 'spe.h5ad')
plot_dir = here('plots', '17_cell_cell_comm', 'liana')

ad_sp = sc.read_h5ad(ad_sp_path)

#   Must use gene symbols for LIANA+
ad_sp.var.index = ad_sp.var['gene_name'].astype(str)
ad_sp.var.index.name = None

plot, _ = li.ut.query_bandwidth(
    coordinates=ad_sp.obsm['spatial'], start=0, end=500, interval_n=20
)
plot.show()
plt.savefig(os.path.join(plot_dir, 'neighbors_vs_radius.pdf'))
plt.close('all')

lr_df = li.mt.bivariate(
    ad_sp,
    resource_name='consensus', # NOTE: uses HUMAN gene symbols!
    local_name='cosine', # Name of the function
    global_name='morans', # Name global function
    n_perms=100, # Number of permutations to calculate a p-value
    mask_negatives=False, # Whether to mask LowLow/NegativeNegative interactions
    add_categories=True, # Whether to add local categories to the results
    nz_prop=0.2, # Minimum expr. proportion for ligands/receptors and their subunits
    use_raw=False,
    verbose=True
)
