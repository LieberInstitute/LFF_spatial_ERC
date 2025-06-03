## Louise Huuki-Myers, May 2025
## Run  dreamlet on sn data

library("SpatialExperiment")
library("dreamlet")
library("here")
library("sessioninfo")

#### Set up dirs ####
data_dir <- here("processed-data", "10_dreamlet_sn", "01_prep_dreamlet_sn")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "10_dreamlet_sn", "01_prep_dreamlet_sn")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC_subcluster"))

colnames(sce) <- paste0(sce$sample_id, "_", sce$Barcode)
rownames(sce) <- rowData(sce)$gene_name

## duplicated gene_names
sce <- sce[!duplicated(rownames(sce)),]

AD_risk <- read.csv(here("processed-data", "00_project_prep", "07_OpenTargets_AD_data", "clin_var_genes.csv")) |>
    dplyr::filter(symbol %in% rowData(sce)$gene_name) 

#### Aggregate to pseudobulk ####
message(Sys.time(), ' - Aggregate to pseudobulk')
pb <- aggregateToPseudoBulk(sce,
                            assay = "counts",
                            cluster_id = "cell_type_broad",
                            sample_id = "sample_id",
                            verbose = FALSE
)

## one assay per cell type
assayNames(pb)

dim(pb)

# Normalize and apply voom/voomWithDreamWeights
res.proc <- processAssays(pb, ~APOE, min.count = 5)

# the resulting object of class dreamletProcessedData stores
# normalized data and other information
res.proc

details(res.proc)

# view details of dropping samples
details(res.proc)

## save 
saveRDS(res.proc, file = here(data_dir, "sn_res_proc.rds"))

# show voom plot for each cell clusters
pdf(here(plot_dir, "sn_dreamlet_voom.pdf"))
plotVoom(res.proc)
dev.off()

# slurmjobs::job_single('01_prep_dreamlet_sn', create_shell = TRUE, memory = '50G', command = "Rscript 01_prep_dreamlet_sn.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

