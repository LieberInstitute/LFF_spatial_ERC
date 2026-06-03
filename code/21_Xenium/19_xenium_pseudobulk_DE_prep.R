## Louise Huuki-Myers, May 2026
## Run spatial registration on 

#### Set Up ####

library("qs2")
library("spatialLIBD")
library("dplyr")
library("here")
library("sessioninfo")
library("getopt")
library("scater")


data_dir <- here("processed-data", "21_Xenium", "19_xenium_pseudobulk_DE_prep")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)


# Import command-line parameters
scec <- matrix(
    c(  "cluster", "c", "1", "character", "Name of cluster"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)
print(opt)

# opt$cluster <- "cell_type_anno"

spe_fn <- here("processed-data", "21_Xenium", "10_xenium_cell_types","spe_xenium_cell_types.qs2")
message(Sys.time(), " - load data from: ", basename(spe_fn))
(sce <- qs_read(spe_fn))

message("ncells: ", ncol(sce))

## make APOE syntatic
sce$APOE_syn <- gsub("/", ".", sce$APOE)


#### run pseudobulk - cell_type_anno ####
message(Sys.time(), " - pseudobulk ", opt$cluster)
sce_pseudo <- registration_pseudobulk(
    sce,
    var_registration = opt$cluster,
    var_sample_id = "sample_id",
    covars = NULL,
    min_ncells = 10,
    pseudobulk_rds_file = NULL,
    filter_expr = FALSE
)

message(Sys.time(), " - Done pseudobulk")

message(sprintf("nrow: %d, ncol: %d", nrow(sce_pseudo), ncol(sce_pseudo)))

#### Check n samples for each cell type ####
table(sce_pseudo$APOE_carrier, sce_pseudo$cell_type_broad)
table(sce_pseudo$APOE_carrier, sce_pseudo$registration_variable)

cell_type_count <- colData(sce_pseudo) |> as.data.frame() |> dplyr::count(APOE_carrier, registration_variable)

## cell type must have two or more samples on either side of DEG split (APOE carrier)
enough_samples <- cell_type_count |>
    filter(n >=2) |>
    dplyr::count(registration_variable) |>
    filter(n >=2) |>
    dplyr::pull(registration_variable)

message("Too few samples in: ", 
        paste(levels(sce_pseudo$registration_variable)[!levels(sce_pseudo$registration_variable) %in% enough_samples], collapse = ", ")
)

cell_type_count |>
    mutate(enough_samples = registration_variable %in% enough_samples)|>
    write.csv(here(data_dir, sprintf("spe_xenium_psuedobulk_sample_count-%s.csv", opt$cluster)))

## drop too few sample cell types
sce_pseudo <- sce_pseudo[, sce_pseudo$registration_variable %in% enough_samples]
sce_pseudo$registration_variable <- droplevels(sce_pseudo$registration_variable)

#### Add PCAs ####
sce_pseudo <- scater::runPCA(sce_pseudo, 
                             ncomponents = 50,
                             name = "PCA")


#### Additional edits + Save ####
## drop all NA cols
all_na <- sapply(colData(sce_pseudo), function(x)all(is.na(x)))
colData(sce_pseudo) <- colData(sce_pseudo)[, names(all_na)[!all_na]]

## save 
message(Sys.time(), " - Save")
saveRDS(sce_pseudo, file = here(data_dir, sprintf("spe_xenium_pseudo_DGE-%s.RDS", opt$cluster)))

# slurmjobs::job_single('19_xenium_pseudobulk_DE_prep', create_shell = TRUE, memory = '10G', command = "Rscript 19_xenium_pseudobulk_DE_prep.R --var cell_type_anno")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()


