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

## add syntactic APOE vars
sce$APOE_syn <- factor(gsub("/", ".", sce$APOE))
levels(sce$APOE_syn)

sce$APOE_carrier_syn <- factor(gsub("\\+", "", sce$APOE_carrier))
levels(sce$APOE_carrier_syn)

## add E4/E4 variable
sce$APOE_E4E4 <- sce$APOE == "E4/E4"
table(sce$APOE_E4E4)

## duplicated gene_names
sce <- sce[!duplicated(rownames(sce)),]

#### Aggregate to pseudobulk ####
message(Sys.time(), ' - Aggregate to pseudobulk')
pb <- aggregateToPseudoBulk(sce,
                            assay = "counts",
                            cluster_id = "cell_type_broad",
                            sample_id = "sample_id",
                            verbose = FALSE
)


## save 
saveRDS(pb, file = here(data_dir, "sn_dreamlet_pb.rds"))

# pb <- readRDS(here("processed-data", "10_dreamlet_sn", "01_prep_dreamlet_sn", "sn_dreamlet_pb.rds"))

# slurmjobs::job_single('01_prep_dreamlet_sn', create_shell = TRUE, memory = '50G', command = "Rscript 01_prep_dreamlet_sn.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

