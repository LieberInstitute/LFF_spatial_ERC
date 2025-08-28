#   Ideally would've been part of 08_liana_bivariate_plot.py-- export spot-level
#   LR scores as a large CSV for analysis outside of python

import pandas as pd
import scanpy as sc
from pyhere import here
import os
import re
import session_info

in_dir = here(
    'processed-data', '17_cell_cell_comm', 'liana', 'bivariate_adatas'
)
out_path = here(
    'processed-data', '17_cell_cell_comm', 'liana', 'bivariate_stats_spot.csv.gz'
)

#   Read input files
in_files = [
    os.path.join(in_dir, f)
    for f in os.listdir(in_dir)
    if re.compile(r'.*\.h5ad$').match(f)
]

pair_df_list = []
for f in in_files:
    adata = sc.read(f)
    
    pair_df = pd.DataFrame(
        adata.X.toarray(), columns = adata.var.index,
        index = adata.obs.index
    )
    pair_df[['key', 'sample_id', 'SpD']] = adata.obs[['key', 'sample_id', 'SpD']]
    pair_df = pair_df.melt(
        id_vars = ['key', 'sample_id', 'SpD'], value_vars = list(adata.var.index),
        var_name = 'pair_id', value_name = 'local_score'
    )
    pair_df_list.append(pair_df)

pair_df = pd.concat(pair_df_list, axis=0, ignore_index=True)
pair_df.to_csv(out_path, index = False)

session_info.show()
