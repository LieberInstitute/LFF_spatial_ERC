## April 2025, Louise Huuki-Myers
## Check sn data against external data 

library("spatialLIBD")
library("tidyverse")
# library("ComplexHeatmap")
# library("ggrepel")
library("here")
library("sessioninfo")

#### Set up dirs ####
data_dir <- here("processed-data", "04_snRNA-seq", "24_external_data_check")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "04_snRNA-seq", "24_external_data_check")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## load colors
load(here("processed-data", "04_snRNA-seq", "cell_type_colors.Rdata"), verbose = TRUE)
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

#### Franjic et al. ####

load("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/snRNAseq_hpc/processed-data/sce/sce_sestan_all.rda", verbose = TRUE)
# sce_sestan

table(sce_sestan$region)

## subset to EC
sce_sestan <- sce_sestan[,sce_sestan$region == "EC"]
sce_sestan$cluster <- droplevels(sce_sestan$cluster)
table(sce_sestan$cluster)

table(sce_sestan$cluster, sce_sestan$samplename)

rowData(sce_sestan)
colData(sce_sestan)

sestan_modeling <- registration_wrapper(
    sce = sce_sestan,
    var_registration = "cluster",
    var_sample_id = 'samplename',
    covars = NULL,
    gene_ensembl = "gene_id",
    gene_name = "gene_name",
    min_ncells = 10,
    pseudobulk_rds_file = NULL
)

colnames(sestan_modeling$enrichment)

saveRDS(sestan_modeling, file = here(data_dir, "sestan_EC_modeling.rds"))

