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

plot_dir <- here("plots", "21_Xenium", "19_xenium_pseudobulk_DE_prep")
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

## make APOE syntatic
spe$APOE_syn <- gsub("/", ".", spe$APOE)

## Simplify SpX

spe$SpX_simple <- factor(gsub("~SpX[0-9]", "", spe$SpX), levels = c("Vasc", "L1a", "L1b", "L2.3", "Inhib", "L5", "L6", "WMtz", "WM"))
table(spe$SpX, spe$SpX_simple)

if(opt$cluster == "cell_type_anno_SpX"){
    
    ## filter to Singlets
    spe <- spe[,spe$spot_class == "singlet"]
    message("filter to singlets ncells: ", ncol(spe))
    
    spe$cell_type_anno_SpX <- paste0(spe$cell_type_anno, "_", spe$SpX_simple)
    
    message("cell type x SpX combindations: ", length(unique(spe$cell_type_anno_SpX)))
    
} else if(opt$cluster == "doublet_Astro_Oligo.3"){
    
    ## filter to doublets
    spe <- spe[,spe$spot_class == "doublet_certain"]
    message("filter to doublet_certain, ncells: ", ncol(spe))
    
    ## filter to Astro-Oligo.3 doublets
    
    spe <- spe[,grepl("Oligo.3|Astro", spe$first_type)]
    spe <- spe[,grepl("Oligo.3|Astro", spe$second_type)]
    
    pd_temp <- colData(spe) |>
        as.data.frame() |>
        mutate(first_type = as.character(first_type),
               second_type = as.character(second_type),
               doublet_fine = paste(pmin(first_type, second_type), 
                             pmax(first_type, second_type), 
                             sep="_"),
               doublet = gsub("Astro\\.[1-5]", "Astro", doublet_fine)
               )
    
    
    spe$doublet_fine <- pd_temp$doublet_fine
    spe$doublet <- pd_temp$doublet
    spe$doublet_Astro_Oligo.3 <- pd_temp$doublet
    
    table(spe$doublet)
}

message("check input: ")
table(spe$spot_class)

table(spe$cell_type_anno)


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

cell_type_count <- colData(spe_pseudo) |> 
    as.data.frame() |> 
    dplyr::group_by(APOE_carrier, registration_variable) |>
    summarise(n = n(),
              min_ncells = min(ncells),
              median_ncells = median(ncells),
              max_ncells = max(ncells),
              )

cell_type_count_wide <- cell_type_count |> 
    select(-ends_with("ncells")) |> 
    tidyr::pivot_wider(values_from = "n", names_from = "APOE_carrier")

## cell type must have two or more samples on either side of DEG split (APOE carrier)
enough_samples <- cell_type_count |>
    ungroup() |>
    filter(n >=2) |>
    dplyr::count(registration_variable) |>
    filter(n >=2) |>
    dplyr::pull(registration_variable)

message("Too few samples in: ", 
        paste(unique(spe_pseudo$registration_variable)[!unique(spe_pseudo$registration_variable) %in% enough_samples], collapse = ", ")
)

cell_type_count |>
    mutate(enough_samples = registration_variable %in% enough_samples)|>
    write.csv(here(data_dir, sprintf("spe_xenium_psuedobulk_sample_count-%s.csv", opt$cluster)), row.names = FALSE)

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


#### Explore number of cells + donors ####

if(opt$cluster == "cell_type_anno_SpX"){
    # spe_pseudo <- readRDS(here("processed-data", "21_Xenium", "19_xenium_pseudobulk_DE_prep", "spe_xenium_pseudo_DGE-cell_type_anno_SpX.RDS"))
    
    n_cells_tb <- colData(spe_pseudo) |>
        as.data.frame() |>
        select(cell_type_anno_SpX, cell_type_anno, SpX, ncells) |>
        as_tibble()
    
    n_cells_tb |> 
        group_by(cell_type_anno) |>
        count(ncells < 50)
    
    n_cell_boxplot <- n_cells_tb |>
        ggplot(aes(x = SpX, y = ncells)) +
        geom_boxplot() +
        facet_wrap(~cell_type_anno) +
        ylim(0, 100)
    
    ggsave(n_cell_boxplot, filename = here(plot_dir, "xenium_pseudobulk_n_cell_boxplot.png"))
    
    
    min_cell_tile <- n_cells_tb |>
        group_by(SpX, cell_type_anno) |>
        slice_min(ncells) |>
        ggplot(aes(x = SpX, y = cell_type_anno, fill = ncells)) +
        geom_tile() +
        geom_text(aes(label = ncells, color = ncells < 50)) +
        theme_bw()  +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    
    ggsave(min_cell_tile, filename = here(plot_dir, "xenium_pseudobulk_min_cell_tile.png"))    
    
    median_cell_tile <- n_cells_tb |>
        group_by(SpX, cell_type_anno) |>
        summarise(median_cells = median(ncells)) |>
        ggplot(aes(x = SpX, y = cell_type_anno, fill = median_cells)) +
        geom_tile() +
        geom_text(aes(label = median_cells, color = median_cells < 50), size = 2) +
        theme_bw()  +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    
    ggsave(median_cell_tile, filename = here(plot_dir, "xenium_pseudobulk_median_cell_tile.png"))
    
    
    n_cell_histo <- n_cells_tb |>
        ggplot(aes(x =  ncells)) +
        geom_histogram(binwidth = 25) +
        facet_grid(cell_type_anno~SpX) +
        ylim(0, 100)
    
    ggsave(n_cell_histo, filename = here(plot_dir, "xenium_pseudobulk_n_cell_histo_grid.png"), height = 25, width = 9)
    
}


# slurmjobs::job_single('19_xenium_pseudobulk_DE_prep', create_shell = TRUE, memory = '10G', command = "Rscript 19_xenium_pseudobulk_DE_prep.R --cluster cell_type_anno")
# slurmjobs::job_single('19_xenium_pseudobulk_DE_prep_cell_type_SpX', create_shell = TRUE, memory = '10G', command = "Rscript 19_xenium_pseudobulk_DE_prep.R --cluster cell_type_anno_SpX")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()


