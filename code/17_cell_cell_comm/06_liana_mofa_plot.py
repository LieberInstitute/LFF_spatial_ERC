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

mdata_in_path = here(
    'processed-data', '17_cell_cell_comm', 'liana', 'mdata.h5mu'
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
plt.savefig(os.path.join(plot_dir, 'mofa_heatmap.pdf'))
plt.close('all')
