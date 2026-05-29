import re
import anndata as ad
import pandas as pd
from pyhere import here

data_dir = here('processed-data', '24_xenium_liana', '03_inflow_score')

dfs = []
for path in sorted(data_dir.glob('lrdata_*.h5ad')):
    sample_id = re.search(r'lrdata_(.+)\.h5ad', path.name).group(1)
    adata = ad.read_h5ad(path)
    df = adata.uns['global_interactions'].copy()
    df['sample_id'] = sample_id
    dfs.append(df)

global_interactions = pd.concat(dfs, ignore_index=True)
