import numpy as np
import pandas as pd
import scanpy as sc
import plotnine as p9
import liana as li
import muon as mu
import mofax as mofa
import decoupler as dc
from pyhere import here
import matplotlib.pyplot as plt
import os
import mudata

task_id = int(os.getenv('SLURM_ARRAY_TASK_ID')) - 1
data_description = ['broad', 'fine', 'visium'][task_id]

mdata_in_path = here(
    'processed-data', '17_cell_cell_comm', 'liana', 'mofa_results',
    f'mdata_{data_description}.h5mu'
)
plot_dir = here('plots', '17_cell_cell_comm', 'liana')

os.makedirs(plot_dir, exist_ok=True)

mdata = mudata.read(mdata_in_path)
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
