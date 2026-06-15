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
    c("cluster", "c", "1", "character", "Name of cluster"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)
print(opt)

# opt$cluster <- "cell_type_anno_SpX"

spe_fn <- here("processed-data", "21_Xenium", "13_xenium_bansky_embedding","spe_xenium_bansky.qs2")
message(Sys.time(), " - load data from: ", basename(spe_fn))

spe <- qs_read(spe_fn)

spe <- spe[,spe$spot_class == "singlet"]
(spe <- qs_read(spe_fn))

message("ncells: ", ncol(spe))

## make APOE syntatic
spe$APOE_syn <- gsub("/", ".", spe$APOE)

## Simplify SpX

spe$SpX_simple <- factor(gsub("~SpX[0-9]", "", spe$SpX), levels = c("Vasc", "L1a", "L1b", "L2.3", "Inhib", "L5", "L6", "WMtz", "WM"))
table(spe$SpX, spe$SpX_simple)

if(opt$cluster == "cell_type_anno_SpX"){
    spe$cell_type_anno_SpX <- paste0(spe$cell_type_anno, "_", spe$SpX_simple)
    
    message("cell type x SpX combindations: ", length(unique(spe$cell_type_anno_SpX)))
}



#### run pseudobulk - cell_type_anno ####
message(Sys.time(), " - pseudobulk ", opt$cluster)
spe_pseudo <- registration_pseudobulk(
    spe,
    var_registration = opt$cluster,
    var_sample_id = "sample_id",
    covars = NULL,
    min_ncells = 10,
    pseudobulk_rds_file = NULL,
    filter_expr = FALSE
)

message(Sys.time(), " - Done pseudobulk")

message(sprintf("nrow: %d, ncol: %d", nrow(spe_pseudo), ncol(spe_pseudo)))

#### Check n samples for each cell type ####
table(spe_pseudo$APOE_carrier, spe_pseudo$cell_type_broad)
table(spe_pseudo$APOE_carrier, spe_pseudo$registration_variable)

cell_type_count <- colData(spe_pseudo) |> as.data.frame() |> dplyr::count(APOE_carrier, registration_variable)

cell_type_count_wide <- cell_type_count |> tidyr::pivot_wider(values_from = "n", names_from = "APOE_carrier")

## cell type must have two or more samples on either side of DEG split (APOE carrier)
enough_samples <- cell_type_count |>
    filter(n >=2) |>
    dplyr::count(registration_variable) |>
    filter(n >=2) |>
    dplyr::pull(registration_variable)

message("Too few samples in: ", 
        paste(unique(spe_pseudo$registration_variable)[!unique(spe_pseudo$registration_variable) %in% enough_samples], collapse = ", ")
)

cell_type_count |>
    mutate(enough_samples = registration_variable %in% enough_samples)|>
    write.csv(here(data_dir, sprintf("spe_xenium_psuedobulk_sample_count-%s.csv", opt$cluster)))

## drop too few sample cell types
spe_pseudo <- spe_pseudo[, spe_pseudo$registration_variable %in% enough_samples]
if(is.factor(spe_pseudo$registration_variable)){
    spe_pseudo$registration_variable <- droplevels(spe_pseudo$registration_variable)
}

#### Add PCAs ####
spe_pseudo <- scater::runPCA(spe_pseudo, 
                             ncomponents = 50,
                             name = "PCA")


#### Additional edits + Save ####
## drop all NA cols
all_na <- sapply(colData(spe_pseudo), function(x)all(is.na(x)))
colData(spe_pseudo) <- colData(spe_pseudo)[, names(all_na)[!all_na]]

## save 
message(Sys.time(), " - Save")
saveRDS(spe_pseudo, file = here(data_dir, sprintf("spe_xenium_pseudo_DGE-%s.RDS", opt$cluster)))

# slurmjobs::job_single('19_xenium_pseudobulk_DE_prep', create_shell = TRUE, memory = '10G', command = "Rscript 19_xenium_pseudobulk_DE_prep.R --cluster cell_type_anno")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()


