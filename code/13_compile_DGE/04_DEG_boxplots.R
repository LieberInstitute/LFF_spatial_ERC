## Louise Huuki-Myers, July 2025
## Plot boxplots for DGE data

#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")
library("getopt")
library("DeconvoBuddies")
library("SingleCellExperiment")

# Import command-line parameters
scec <- matrix(
    c("datatype", "d", "1", "character", "Data type"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)

## test
# opt$datatype <- "sn_broad"

print(opt)

#### set up dirs ####
# data_dir <- here("processed-data", "13_compile_DGE", "04_DEG_boxplots")
# if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "13_compile_DGE", "04_DEG_boxplots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load the data ####
if(opt$datatype == "Visium"){
    pb_fn <- here("processed-data", "09_pseudoBulkDGE_Visium", "01_pseudobulk_data_Visium", "spe_pseudo_DGE.RDS")
    batch <- "Visium_slide"
} else if(opt$datatype == "sn_broad"){
    pb_fn <- here("processed-data", "08_pseudoBulkDGE_sn", "01_pseudobulk_data_sn","sce_pseudo_DGE-cell_type_broad.RDS")
    batch <- "exp_round"
}else if(opt$datatype == "sn_fine"){
    pb_fn <- here("processed-data", "08_pseudoBulkDGE_sn", "01_pseudobulk_data_sn","sce_pseudo_DGE-cell_type_anno.RDS")
    batch <- "exp_round"
} else {
    stop("non-valid datatype")
}

message(Sys.time(), sprintf(" - Datatype = %s, loading '%s'", opt$datatype, basename(pb_fn)))
sce_pb <- readRDS(pb_fn)
dim(sce_pb)
table(sce_pb$registration_variable)

rownames(sce_pb) <- rowData(sce_pb)$gene_name

load(here("processed-data", "project_colors.Rdata"))

## load DE data ##

DE_data_fn <- here("processed-data", "13_compile_DGE", "01_compile_DGE", opt$datatype, sprintf("DGE_results_carrier_%s.Rds", opt$datatype))
DE_data <- readRDS(DE_data_fn)

DEGs_signif <- DE_data |> 
    filter(vlmf_adj.P.Val < 0.05) |>
    group_by(cluster) |>
    arrange(vlmf_adj.P.Val) |>
    mutate(DE_class = case_when(vlmf_logFC > 0 ~ "up",
                                vlmf_logFC < 0 ~ "down",
                                TRUE ~ "None"),
           DE_class_cluster = paste0(gsub("\\.", "-", cluster), "_",DE_class),
           rank = row_number())

DEGs_signif |> filter(rank <= 5)


"COL4A5" %in% rownames(sce_pb)

test <- plot_gene_express(sce = sce_pb[sce_pb$registration_variable == "Inhib",],
                          genes = "ST18",
                          # genes = "CAMTA1",
                          category = "APOE",
                          plot_points = TRUE)

ggsave(test, filename = here(plot_dir, "test.png"))



