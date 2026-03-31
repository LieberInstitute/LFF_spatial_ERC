import scanpy as sc
import decoupler as dc
from pyhere import here
import matplotlib.pyplot as plt
import os
import mudata
import session_info
import pandas as pd

task_id = int(os.getenv('SLURM_ARRAY_TASK_ID')) - 1
data_description = ['broad', 'fine', 'visium'][task_id]

mdata_in_path = here(
    'processed-data', '17_cell_cell_comm', 'liana', 'mofa_results',
    f'mdata_{data_description}.h5mu'
)
path_path = here(
    "processed-data", "00_project_prep", "05_pathology", "sample_taupathy.csv"
)
plot_dir = here('plots', '17_cell_cell_comm', 'liana')

os.makedirs(plot_dir, exist_ok=True)

#   Read in results as an AnnData and append some pathology covariates
path_df = pd.read_csv(path_path, index_col='BrNum')
mdata = mudata.read(mdata_in_path)
mdata.obs['taupathy'] = pd.Series(
    ['t+' if x else 't-' for x in path_df['taupathy']], index = path_df.index
)
adata = sc.AnnData(obs = mdata.obs, obsm = mdata.obsm)

dc.tl.rankby_obsm(adata, key='X_mofa', uns_key='rank_obsm')
dc.pl.obsm(
    adata,
    key='rank_obsm',
    names = list(mdata.obs.columns),
    titles=['Principle component scores', 'Adjusted p-values from ANOVA'],
    figsize=(10, 7),
    nvar=10
)
plt.savefig(os.path.join(plot_dir, f'mofa_heatmap_{data_description}.pdf'))
plt.close('all')

session_info.show()
