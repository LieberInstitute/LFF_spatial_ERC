import scanpy as sc
import liana as li
import plotnine as p9
import os
from pyhere import here
import session_info

ad_path = here(
    'processed-data', '24_xenium_liana', '01_to_python',
    'cell_rna_converted_gobject.h5ad'
)
plot_dir = here('plots', '24_xenium_liana', '02_bandwidth')

os.makedirs(plot_dir, exist_ok=True)

#   Grab just the first sample. Spatial coordinates are in microns for Xenium,
#   so units that apply to one sample should generally apply decently to the
#   full dataset
adata = sc.read_h5ad(ad_path)
adata = adata[
    adata.obs['sample_id'] == adata.obs['sample_id'].unique()[0], :
]

plot, df = li.ut.query_bandwidth(
    coordinates=adata.obsm["spatial"],
    start=10, end=50,           
    interval_n=40    
)

#   Show that a bandwidth of 42 microns yields about 5 neighbors
bandwidth = 42
plot = (
    plot + 
        p9.geom_vline(xintercept=bandwidth, linetype="dashed", color="red") + 
        p9.scale_y_continuous(
            breaks=range(
                int(df.neighbours.min()), int(df.neighbours.max()) + 1
            )
        )
)
plot.save(filename=str(plot_dir / "bandwidth_neighbors.pdf"), verbose=False)

#   Compute neighbors (just for visualization purposes)
li.ut.spatial_neighbors(
    adata=adata, bandwidth=bandwidth, spatial_key="spatial"
)

#   Plot what the spatial kernel looks like on the example sample
p = (
    li.pl.connectivity(
        adata, idx=5500, size=0.01, figure_size=(6, 5), spatial_key="spatial"
    ) +
    p9.coord_fixed()
)
p.save(filename=str(plot_dir / "spatial_kernel.pdf"), verbose=False)

session_info.show()
