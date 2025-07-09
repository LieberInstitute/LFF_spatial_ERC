## Louise Huuki-Myers, July 2025
## Detirmine optimal resolution for DE

#### set-up ####
library("edgeR")
library("limma")
library("SingleCellExperiment")
library("here")
library("data.table")
library("getopt")
library("tidyverse")
library("sessioninfo")
library("spatialLIBD")
library("dendextend")
library("dynamicTreeCut")

# Import command-line parameters
scec <- matrix(
    c("celltype", "c", "1", "character", "Data type"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)

## test
# opt$celltype <- "Astro"

## set up output dirs
data_dir <- here("processed-data", "12_voomLmFit", "02_optimal_cluster_voomLmFit", opt$celltype)
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "12_voomLmFit", "02_optimal_cluster_voomLmFit", opt$celltype)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC_subcluster"))

sce <- sce[,sce$cell_type_broad == opt$celltype]

sce$cell_type_anno <- droplevels(sce$cell_type_anno)
table(sce$cell_type_anno)


## hierarchical cluster data
load(here("processed-data", "04_snRNA-seq", "32_sn_subcluster_hierarchical_cluster", "sn_subcluster_hierarchical_cluster.Rdata"), verbose = TRUE)
# dend
# tree.clusCollapsed
# dist.clusCollapsed

## prune tree to cell types
prune_ct <- labels(dend)[!labels(dend) %in% levels(sce$cell_type_anno)]

dend_ct <- prune(dend, prune_ct)

str(dend_ct)
## get internal node heights
node_heights <- get_nodes_attr(dend_ct, "height")[is.na(get_nodes_attr(dend_ct, "label"))]
cut_heights <- sort(node_heights)[-length(node_heights)]

pdf(here(plot_dir, sprintf("dend_%s.pdf", opt$celltype)))
plot(dend_ct)
abline(h = cut_heights, lty = 2)
dev.off()

hc_ct <- map_dfr(cut_heights, ~cutree(dend_ct, h = .x))
hc_ct$hc.n <- rowMax(as.matrix(hc_ct))

hc_ct_long <- hc_ct |>
    pivot_longer(!hc.n, names_to = "cell_type_anno", values_to = "subgroup") |>
    rowwise() |>
    mutate(hc.ct = paste0("cell_type_hc", hc.n),
           cell_type_hc = paste0(c(opt$celltype, hc.n, subgroup), collapse = ".")
           )
## create new labels
hc_ct_wide <- hc_ct_long |> 
    select(cell_type_anno, cell_type_hc, hc.ct) |>
    pivot_wider(names_from = "hc.ct", values_from = "cell_type_hc")

## add to sce object
hc_ct_wide_sce <- hc_ct_wide |>
    column_to_rownames("cell_type_anno") |>
    mutate_all(factor)

hc_ct_wide_sce <- hc_ct_wide_sce[sce$cell_type_anno, ]

stopifnot(nrow(hc_ct_wide_sce) == ncol(sce))

colData(sce) <- cbind(colData(sce), hc_ct_wide_sce)

## check tables
hc_levels <- colnames(hc_ct_wide_sce)
names(hc_levels) <- hc_levels
map(hc_levels, ~table(sce$cell_type_anno, sce[[.x]]))

#### pseudobulk + run DGE ####

source(here("code", "utils", "sce_pseudobulk.R"))
source(here("code", "utils", "clusterwise_voomLmFit.R"))

FDR_summary <- map(hc_levels, function(hc_level){
    
    ## pseudobulk
    sce_pb <- sce_pseudobulk(sce, hc_level)
    
    ## run DGE by subtype
    hc_cell_types <- levels(sce_pb$registration_variable)
    names(hc_cell_types) <- hc_cell_types
    
    vlmf_tt <- map(hc_cell_types, ~clusterwise_voomLmFit(sce_pb, clus = .x))
    
    saveRDS(vlmf_tt, file = here(data_dir, sprintf("voomLmFit_%s.rds", hc_level)))
    
    return(map_int(vlmf_tt, ~sum(.x$adj.P.Val < 0.05)))
})

FDR_summary_tb <- tibble(cell_type_hc = names(unlist(FDR_summary)),
                      n_FDR05 = unlist(FDR_summary)) |>
    separate(cell_type_hc, sep = "\\.", into = c("hc", "cell_type_hc"), extra = "merge") |>
    mutate(cell_type_broad = opt$celltype, .before = 1)

write_csv(FDR_summary_tb, file = here("processed-data", "12_voomLmFit", "02_optimal_cluster_voomLmFit", 
                                   sprintf("optimal_cluster_FDR05_summary-%s.csv", opt$celltype)))

# slurmjobs::job_loop(loops = list(celltype = c("Astro", "Micro", "OPC", "Oligo", "Vasc",  "Excit", "Inhib")),
#                     name = "02_optimal_cluster_voomLmFit",
#                     create_shell = TRUE,
#                     memory = "25G",
#                     create_script = FALSE)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


