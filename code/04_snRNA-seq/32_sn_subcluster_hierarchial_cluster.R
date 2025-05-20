## Louise Huuki-Myers, May 2025
## Run hierchial clustering on sn subclustering to examine different resolutions

library("SingleCellExperiment")
library("dendextend")
library("dynamicTreeCut")
library("here")

## Prep directories
plot_dir <- here("plots", "04_snRNA-seq", "32_sn_subcluster_hierarchical_cluster")
if(!dir.exists(plot_dir)) dir.create(plot_dir)

data_dir <- here("processed-data", "04_snRNA-seq", "32_sn_subcluster_hierarchical_cluster")
if(!dir.exists(data_dir)) dir.create(data_dir)


#### load pseudobulk data ####
load(here("processed-data", "04_snRNA-seq", "29_sn_subcluster_model_pseudobulk"))