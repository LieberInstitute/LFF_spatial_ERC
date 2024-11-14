## November 2024, Louise Huuki-Myers
## Run BayesSpace as outlined in the tutorial https://www.bioconductor.org/packages/release/bioc/vignettes/BayesSpace/inst/doc/BayesSpace.html

## Required libraries

library("SpatialExperiment")
library("spatialLIBD")
library("BayesSpace")
library("Polychrome")
library("tidyverse")
library("HDF5Array")

library("getopt")
library("here")
library("sessioninfo")

plot_dir <- here("plots", "05_spe_correct_cluster", "05_BayesSpace_by_the_book")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

data_dir <- here("processed-data", "05_spe_correct_cluster", "05_BayesSpace_by_the_book")
if(!dir.exists(data_dir)) dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

#### Load data ####
message(Sys.time(), " - Load HDF5 SPE")
spe_hdf5_path <- here("processed-data", "spe_objects", "spe_postQC")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here(spe_hdf5_path))

reducedDimNames(spe)

#### PreProcess ####
message(Sys.time(), " - Preprocess")
spe <- spatialPreprocess(spe, 
                         platform="Visium",
                         n.PCs=15,
                         n.HVGs=2000, 
                         log.normalize=TRUE)

reducedDimNames(spe)

#### Clustering ####
## Set Seed
set.seed(20241114)

## Fix colnames for row and col
spe$row <- spe$array_row
spe$col <- spe$array_col

## pick number of clusters
message(Sys.time(), " - qTune")
spe <- qTune(spe, qs=seq(2, 28), platform="Visium")

pdf(here(plot_dir, "BayesSpace_erc_qplot.pdf"))
qPlot(spe)
dev.off()

## choose k/q with plot
k = 9
k_nice <- sprintf("k%02d", k)

## Set the BayesSpace metadata using code from
## https://github.com/edward130603/BayesSpace/blob/master/R/spatialPreprocess.R#L43-L46
metadata(spe)$BayesSpace.data <- list(platform = "Visium", is.enhanced = FALSE)

## Run BayesSpace
message(Sys.time(), " - Running spatialCluster: k=", k, ", dimred = ", dimred)
spe <- BayesSpace::spatialCluster(spe, 
                                  use.dimred = dimred, 
                                  q = k, 
                                  nrep = 20000)

message(Sys.time(), " - Format and Export")

table(spe$spatial.cluster)

spe$bayesSpace_temp <- as.factor(spe$spatial.cluster)
bayesSpace_name <- paste0("BayesSpace_by_the_books_", k_nice)
colnames(colData(spe))[ncol(colData(spe))] <- bayesSpace_name

cluster_export(
    spe,
    bayesSpace_name,
    cluster_dir = here(data_dir, "clusters_BayesSpace_by_the_book")
)

## Visualize BayesSpace results
message(Sys.time(), " - Visualize clusters")
sample_ids <- unique(spe$sample_id)

## create color pallet
cols <- Polychrome::palette36.colors(k)
if(k ==2) cols <- cols[seq(2)] ## fix return 3 colors bug

names(cols) <- sort(unique(spe[[bayesSpace_name]]))

#   Use 'vis_grid_clus' to preserve all spots (including overlaps)
p_list = vis_grid_clus(
    spe = spe,
    clustervar = bayesSpace_name,
    pdf_file = here(dir_plots, paste0(bayesSpace_name, "-ALL.pdf")),
    sort_clust = FALSE,
    colors = cols,
    spatial = FALSE,
    point_size = 1
)

# #  Vis_clus for each sample
# walk(sample_ids, function(samp){
#     spot_plot <- vis_clus(
#         spe = spe,
#         sampleid = samp,
#         clustervar = bayesSpace_name,
#         # sort_clust = FALSE,
#         colors = cols,
#         point_size = 2
#     )
#     ggsave(spot_plot, filename = here(dir_plots, paste0(bayesSpace_name, "-", samp, ".pdf")))
# })

# slurmjobs::job_single('05_BayesSpace_by_the_book', create_shell = TRUE, memory = '25G', command = "Rscript 05_BayesSpace_by_the_book.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
