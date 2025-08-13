import pandas as pd
import scanpy as sc
import decoupler as dc
import liana as li
from matplotlib import pyplot as plt
from mudata import MuData
from pyhere import here
import os
import json

ad_sp_path = here('processed-data', '17_cell_cell_comm', 'spe.h5ad')
sample_info_path = here(
    'processed-data', '00_project_prep', '02_get_online_metadata',
    'metadata_visium_plan.csv'
)
scale_path = here(
    'processed-data', '01_spaceranger', '{}', 'outs', 'spatial',
    'scalefactors_json.json'
)
plot_dir = here('plots', '17_cell_cell_comm', 'liana')
spot_diameter_um = 65
inter_spot_dist_um = 100

sample_info = pd.read_csv(sample_info_path)
sample_info['sample_id'] = sample_info['Visium Slide #'] + '_' + sample_info['Visium_subslide']

ad_sp = sc.read_h5ad(ad_sp_path)

#   Must use gene symbols for LIANA+
ad_sp.var.index = ad_sp.var['gene_name'].astype(str)
ad_sp.var.index.name = None

#   Subset to just this one sample
donor_id = ad_sp.obs['sample_id'].cat.categories[0]
sample_id = sample_info.loc[sample_info['BrNum'] == donor_id, 'sample_id'].values[0]
ad_sp_small = ad_sp[ad_sp.obs['sample_id'] == donor_id, :]
scale_path = str(scale_path).format(sample_id)

with open(scale_path, 'r') as f:
    scale_factors = json.load(f)

#   Calculate the bandwidth (radius from the center of a spot), in number of
#   pixels, that just barely encompasses the 6 spot neighbors
good_bandwidth = scale_factors['spot_diameter_fullres'] * (0.5 + inter_spot_dist_um / spot_diameter_um)

plot, _ = li.ut.query_bandwidth(coordinates=ad_sp_small.obsm['spatial'], start=0, end=500, interval_n=20)
plot.show()
plt.savefig(os.path.join(plot_dir, 'aaa.pdf'))
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
