## Louise Huuki-Myers, May 2025
## Run  dreamlet on sn data

library("SpatialExperiment")
library("dreamlet")
library("here")
library("sessioninfo")
library("HDF5Array")

#### Set up dirs ####
data_dir <- here("processed-data", "11_dreamlet_Visium", "01_prep_dreamlet_Visium")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "11_dreamlet_Visium", "01_prep_dreamlet_Visium")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC_annotated"))
spe

## Use gene names
rownames(spe) <- rowData(spe)$gene_name
## deal with duplicated gene_names
rownames(spe)[duplicated(rownames(spe))] <- rowData(spe)$gene_id[duplicated(rownames(spe))]
stopifnot(!any(duplicated(rownames(spe))))

## add syntactic APOE vars
spe$APOE_syn <- factor(gsub("/", ".", spe$APOE))
levels(spe$APOE_syn)

spe$APOE_carrier_syn <- factor(gsub("\\+", "", spe$APOE_carrier))
levels(spe$APOE_carrier_syn)

## add E4/E4 variable
spe$APOE_E4E4 <- spe$APOE == "E4/E4"
table(spe$APOE_E4E4)

## syntatic SpD
spe$SpD_syn <- gsub("~", "_", spe$SpD)

## drop uneeded 10x clusters -> lead to error
colData(spe)[,grep("^10x", colnames(colData(spe)))] <- NULL

#### Aggregate to pseudobulk ####
message(Sys.time(), ' - Aggregate to pseudobulk')
pb <- aggregateToPseudoBulk(spe,
                            assay = "counts",
                            cluster_id = "SpD_syn",
                            sample_id = "sample_id",
                            verbose = FALSE
)

## one assay per SpD
assayNames(pb)
dim(pb)


## save 
saveRDS(pb, file = here(data_dir, "Visium_dreamlet_pb.rds"))

# pb <- readRDS(here("processed-data", "11_dreamlet_Visium", "01_prep_dreamlet_Visium", "Visium_dreamlet_pb.rds"))



# slurmjobs::job_single('01_prep_dreamlet_Visium', create_shell = TRUE, memory = '25G', command = "Rscript 01_prep_dreamlet_Visium.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

