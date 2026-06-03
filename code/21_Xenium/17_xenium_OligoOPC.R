## Louise Huuki-Myers, May 2026
## Explore Oligo + OPC in Xenium data

#### Set Up ####
library("SpatialExperiment")
library("qs2")
library("here")
library("sessioninfo")

library("tidyverse")
library("crumblr")
library("variancePartition")
library("spatialLIBD")

data_dir <- here("processed-data", "21_Xenium", "17_xenium_OligoOPC")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "21_Xenium", "17_xenium_OligoOPC")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE) 
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

#### Load data ####
message(Sys.time(), "- Load xenium data")
spe <- qs_read(here("processed-data", "21_Xenium", "13_xenium_bansky_embedding","spe_xenium_bansky.qs2"))

spe <- spe[,spe$spot_class == "singlet"]

#### make Oligo.3 only SpD ####
if(FALSE){
    
    spe_O3 <- spe[,spe$cell_type_anno == 'Oligo.3']
    table(spe_O3$SpX, spe_O3$BrNum)
    
    ## make APOE syntatic
    spe_O3$APOE_syn <- gsub("/", ".", spe_O3$APOE)
    
    ## add syntacticly valid version of SpX
    spe_O3$SpX_syn <- gsub("~", "_", spe_O3$SpX)
    table(spe_O3$SpX_syn)
    
    #### run pseudobulk - SpX ####
    message(Sys.time(), " - pseudobulk ")
    sce_pseudo <- registration_pseudobulk(
        spe_O3,
        var_registration = "SpX_syn",
        var_sample_id = "sample_id",
        covars = NULL,
        min_ncells = 10,
        pseudobulk_rds_file = NULL,
        filter_expr = FALSE
    )
    
    message(Sys.time(), " - Done pseudobulk")
    
    message(sprintf("nrow: %d, ncol: %d", nrow(sce_pseudo), ncol(sce_pseudo)))
    
    #### Check n samples for each cell type ####
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
        write.csv(here(data_dir, "spe_xenium_psuedobulk_sample_count-O3_SpX.csv"))
    
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
    saveRDS(sce_pseudo, file = here(data_dir, "spe_xenium_pseudo_DGE-O3_SpX.RDS"))
    
}




