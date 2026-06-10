#   Grab the inflow scores for a select interaction and export for later use
#   in R, where we have info about spatial domains

import re
import anndata as ad
import pandas as pd
from pyhere import here
import os
import session_info

in_dir = here('processed-data', '24_xenium_liana', '03_inflow_score')
out_path = here(
    'processed-data', '24_xenium_liana', '08_prep_inflow', 'inflow_scores.csv'
)
interaction = 'Astro.5^SPON1^APP'

os.makedirs(out_path.parent, exist_ok=True)

inflow_df_list = []
for path in sorted(in_dir.glob('lrdata_*.h5ad')):
    sample_id = re.search(r'lrdata_(.+)\.h5ad', path.name).group(1)
    adata = ad.read_h5ad(path)

    inflow_df = pd.DataFrame({
        'inflow_score': adata[:, interaction].X.toarray().flatten(),
        'cell_id': adata.obs_names
    })
    inflow_df_list.append(inflow_df)

inflow_df = pd.concat(inflow_df_list, axis=0, ignore_index=True)
inflow_df.to_csv(out_path, index=False)

session_info.show()
