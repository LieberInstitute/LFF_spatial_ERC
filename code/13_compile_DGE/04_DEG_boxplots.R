## Louise Huuki-Myers, July 2025
## Plot boxplots for DGE data

#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")
library("getopt")
library("DeconvoBuddies")
library("SingleCellExperiment")

source(here("code","utils","plot_DEG_express.R"))
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

# Import command-line parameters
scec <- matrix(
    c("datatype", "d", "1", "character", "Data type"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)

## test
# opt$datatype <- "sn_fine"
# opt$datatype <- "sn_broad"
# opt$datatype <- "Visium"

print(opt)

#### set up dirs ####
# data_dir <- here("processed-data", "13_compile_DGE", "04_DEG_boxplots")
# if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "13_compile_DGE", "04_DEG_boxplots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load the data ####
if(opt$datatype == "Visium"){
    pb_fn <- here("processed-data", "09_pseudoBulkDGE_Visium", "01_pseudobulk_data_Visium", "spe_pseudo_DGE.RDS")
    cluster_var <- "SpD"
} else if(opt$datatype == "sn_broad"){
    pb_fn <- here("processed-data", "08_pseudoBulkDGE_sn", "01_pseudobulk_data_sn","sce_pseudo_DGE-cell_type_broad.RDS")
    cluster_var <- "cell_type_broad"
}else if(opt$datatype == "sn_fine"){
    pb_fn <- here("processed-data", "08_pseudoBulkDGE_sn", "01_pseudobulk_data_sn","sce_pseudo_DGE-cell_type_anno.RDS")
    cluster_var <- "cell_type_anno"
} else {
    stop("non-valid datatype")
}

message(Sys.time(), sprintf(" - Datatype = %s, loading '%s'", opt$datatype, basename(pb_fn)))
sce_pb <- readRDS(pb_fn)
dim(sce_pb)
table(sce_pb$registration_variable)

rownames(sce_pb) <- rowData(sce_pb)$gene_name

load(here("processed-data", "project_colors.Rdata"))

#### Carrier DE data ####

DE_data_fn <- here("processed-data", "13_compile_DGE", "01_compile_DGE", opt$datatype, sprintf("DGE_results_carrier_%s.Rds", opt$datatype))
DE_data <- readRDS(DE_data_fn)

DEGs_signif <- DE_data |> 
    group_by(cluster) |>
    arrange(vlmf_adj.P.Val) |>
    mutate(DE_class = case_when(vlmf_logFC > 0 ~ "up",
                                vlmf_logFC < 0 ~ "down",
                                TRUE ~ "None"),
           DE_class_cluster = paste0(gsub("\\.", "-", cluster), "_",DE_class),
           rank = row_number())

DEGs_signif |> filter(vlmf_adj.P.Val < 0.05) |> dplyr::count(cluster)

head(model.matrix(~0 + APOE_syn + Sex + Age + Anc_Afr + pseudo_expr_chrM_ratio, colData(sce_pb)))

cluster_levels <- levels(sce_pb[[cluster_var]])
cluster_levels2 <- as.character(unique(DE_data$cluster))

cluster_levels <- intersect(cluster_levels, cluster_levels2)

all(cluster_levels %in% sce_pb[[cluster_var]])

# plot_DEG_express(sce = sce_pb,
#                     stats = DE_data,
#                     clus = "Astro",
#                     n_genes = 10,
#                     pval_col = "vlmf_adj.P.Val",
#                     fc_col = "vlmf_logFC",
#                     gene_col = "gene_name",
#                     cluster_col = cluster_var,
#                     category_col = "APOE_carrier",
#                     mod = ~0 + APOE_syn + Sex + Age + Anc_Afr + pseudo_expr_chrM_ratio,
#                     color_pal = APOE_carrier_colors,
#                     plot_points = FALSE,
#                     ncol = 2,
#                     cleanY_P = 4)

table(sce_pb[[cluster_var]])

pdf(here(plot_dir, sprintf("DEG_boxplots_carrier_%s.pdf", opt$datatype)))
map(cluster_levels, ~plot_DEG_express(sce = sce_pb,
                                          stats = DE_data,
                                          clus = .x,
                                          n_genes = 10,
                                          pval_col = "vlmf_adj.P.Val",
                                          fc_col = "vlmf_logFC",
                                          gene_col = "gene_name",
                                          cluster_col = cluster_var,
                                          category_col = "APOE_carrier",
                                          mod = ~0 + APOE_syn + Sex + Age + Anc_Afr + pseudo_expr_chrM_ratio,
                                          color_pal = APOE_carrier_colors,
                                          plot_points = TRUE,
                                          ncol = 2,
                                          cleanY_P = 4)
         )
dev.off()

#### plot taupath data ####

## load tau pathology data
tau_tb <- read.csv(here("processed-data", "00_project_prep","05_pathology","sample_taupathy.csv")) |>
    column_to_rownames("BrNum")

sce_pb$taupathy <- ifelse(tau_tb[sce_pb$BrNum,]$taupathy, "t+", "t-")
table(sce_pb$taupathy)

sce_pb$carrier_tau <- paste(sce_pb$APOE_carrier, sce_pb$taupathy)

table(sce_pb$registration_variable, sce_pb$carrier_tau)

carrier_tau_colors <- c(`E2+ t-` = "#398A84",
                        `E2+ t+` = "#60BEB8",
                        `E4+ t-` = "#D46B43",
                        `E4+ t+` = "#DD8A69")

pdf(here(plot_dir, sprintf("DEG_boxplots_carrier_taupathy_%s.pdf", opt$datatype)))
map(cluster_levels, ~plot_DEG_express(sce = sce_pb,
                                         stats = DE_data,
                                         clus = .x,
                                         n_genes = 10,
                                         pval_col = "vlmf_adj.P.Val",
                                         fc_col = "vlmf_logFC",
                                         gene_col = "gene_name",
                                         cluster_col = cluster_var,
                                         category_col = "carrier_tau",
                                         mod = ~0 + APOE_syn + Sex + Age + Anc_Afr + pseudo_expr_chrM_ratio,
                                         color_pal = carrier_tau_colors,
                                         plot_points = TRUE,
                                         ncol = 2,
                                         cleanY_P = 4)
)
dev.off()

#### Carrier DE data ####

DE_interaction_data_fn <- here("processed-data", "13_compile_DGE",  "05_compile_DGE_interaction", opt$datatype, sprintf("DGE_results_interaction_%s.Rds", opt$datatype))
DE_interaction_data <- readRDS(DE_interaction_data_fn)

DE_interaction_data |> filter(vlmf_adj.P.Val < 0.2) |> dplyr::count(cluster)
DE_interaction_data |> filter(vlmf_adj.P.Val < 0.05) 

sce_pb$carrier_anc <- paste(sce_pb$APOE_carrier, sce_pb$Ancestry)

table(sce_pb$carrier_anc)

interaction_mod <- model.matrix(~APOE_carrier_syn*Ancestry + Sex + Age + pseudo_expr_chrM_ratio, colData(sce_pb))
interaction_mod <-interaction_mod[,c("(Intercept)", "APOE_carrier_synE4:AncestryEA", "APOE_carrier_synE4","AncestryEA", "SexM","Age","pseudo_expr_chrM_ratio")]
head(interaction_mod)

## resubset clusters
cluster_levels2 <- as.character(unique(DE_interaction_data$cluster))
cluster_levels <- intersect(cluster_levels, cluster_levels2)


pdf(here(plot_dir, sprintf("DEG_boxplots_interaction_%s.pdf", opt$datatype)))
map(cluster_levels, ~plot_DEG_express(sce = sce_pb,
                                      stats = DE_interaction_data,
                                      clus = .x,
                                      n_genes = 10,
                                      pval_col = "vlmf_adj.P.Val",
                                      fc_col = "vlmf_logFC",
                                      gene_col = "gene_name",
                                      cluster_col = cluster_var,
                                      category_col = "carrier_anc",
                                      mod = interaction_mod,
                                      # color_pal = carrier_tau_colors,
                                      plot_points = TRUE,
                                      ncol = 2,
                                      cleanY_P = 2)
)
dev.off()

# slurmjobs::job_single('04_DEG_boxplots', create_shell = TRUE, memory = '5G', command = "Rscript 04_DEG_boxplots.R --datatype sn_broad")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

