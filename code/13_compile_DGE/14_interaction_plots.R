## Louise Huuki-Myers, April 2026
## Create summary plots for Interaction data

#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")
library("ggrepel")
library("getopt")
library("ComplexHeatmap")
library("SingleCellExperiment")
library("DeconvoBuddies")

## set plot dir
plot_dir <- here("plots", "13_compile_DGE", "14_interaction_plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

load(here("processed-data", "project_colors.Rdata"))
AD_risk <- read.csv(here("processed-data", "00_project_prep", "07_OpenTargets_AD_data", "clin_var_genes.csv")) 

#### colors and factors ####

datatype <- "sn_fine"

cluster_colors <- c()
cluster_levels <- c()

if(datatype == "sn_broad"){
  load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE)
  cluster_colors <- cell_type_colors$broad
  cluster_levels <- names(cell_type_colors$broad)
} else if(datatype == "sn_fine"){
  load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE)
  cluster_colors <- cell_type_colors$anno
  cluster_levels <- names(cell_type_colors$anno)
}else if(datatype == "Visium"){
  load(here("processed-data", "SpD_colors.Rdata"), verbose = TRUE)
  cluster_colors <- SpD_colors
  cluster_levels <- names(SpD_colors)
}

cluster_levels <- cluster_levels[cluster_levels != "Other"]

#### Read summary tables ####

summary_table_fn <- list.files(here("processed-data", "13_compile_DGE", "05_compile_DGE_interaction"), pattern = 'summary', recursive = TRUE, full.names = TRUE)
names(summary_table_fn) <- gsub("vlmf_model_summary_interaction_|.csv", "", basename(summary_table_fn))

dge_count_combined <- map2_dfr(summary_table_fn, names(summary_table_fn), ~read.csv(.x, row.names = 1) |> mutate(datatype = gsub("Age_|Anc_","",.y), .before = 1))

dge_count_combined |> filter(n_FDR05 > 0) |> select(1:4)

# datatype             mod          cluster n_FDR05 nUP nDown n_FDR10 n_FDR20
# 1  sn_fine interaction Age     Excit.L2_5.2     947 702   245    1727    2976
# 2  sn_fine interaction Age       Excit.L5.2       7   7     0       7      13
# 3  sn_fine interaction Age            OPC.3       1   0     1       1       1
# 4  sn_fine interaction Anc          Astro.5       1   1     0       1       1
# 5  sn_fine interaction Anc          Micro.2       1   1     0       1       1
# 6  sn_fine interaction Anc          Micro.4       1   1     0       1       1
# 7  sn_fine interaction Anc Inhib.Chandelier       1   1     0       1       1
# 8   Visium interaction Age    Inhib~Sp09D09       1   1     0       1       1
# 9   Visium interaction Anc       L5~Sp09D03       1   0     1       1       2


#### combined summary barplot ####

dge_summary_bar_reg <- dge_count_combined |>
  filter(n_FDR05 > 0) |>
  select(datatype, cluster, Up = nUP, Down = nDown) |>
  mutate(Down = -1*Down) |>
  pivot_longer(!c(cluster, datatype), names_to = "reg", values_to = "n_genes") |>
  ggplot(aes(x = cluster, y = n_genes, fill = reg)) +
  geom_col() +
  geom_text(aes(y = n_genes + (n_genes/abs(n_genes)*40), label = ifelse(n_genes != 0, abs(n_genes), ""), color = reg)) +
  # geom_text(aes(label = ifelse(n_genes != 0, abs(n_genes), ""), color = reg), 
  #                 size = 3) +
  facet_wrap(~datatype, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = c(Up = APOE_carrier_colors[["E4+"]],
                               Down = APOE_carrier_colors[["E2+"]])) +    
  scale_color_manual(values = c(Up = APOE_carrier_colors_dark[["E4+"]],
                                Down = APOE_carrier_colors_dark[["E2+"]])) +
  theme_bw() +
  labs(title = sprintf("DGE interaction - %s", interaction), y = "n DE genes (FDR < 0.05)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position = "right")

ggsave(dge_summary_bar_reg, filename = here(plot_dir, sprintf("DGE_interaction_%s_summary_bar_reg_combined.png", interaction)), height = 5, width = 6)

#### load DEG data ####
interaction <- "Age"

DE_data <- readRDS(here("processed-data", "13_compile_DGE", "05_compile_DGE_interaction", "sn_fine", sprintf("DGE_results_interaction_%s_sn_fine.Rds", interaction))) |>
  mutate(cluster = factor(cluster, levels = cluster_levels))

#### v gene vs. n signif ####
dge_n_v_fdr <- DE_data |>
  group_by(cluster) |>
  summarize(n_gene = n(),
            n_signif = sum(vlmf_adj.P.Val < 0.05))

dge_n_v_fdr_scatter <- dge_n_v_fdr |>
  ggplot(aes(x = n_gene, y = n_signif, color = cluster)) +
  geom_point() +
  geom_text_repel(aes(label = cluster)) +
  scale_color_manual(values = cluster_colors) +
  # facet_wrap(~contrast) +
  scale_y_log10() +
  theme_bw() +
  theme(legend.position = "None") +
  labs(title = sprintf("DGE interaction - %s", interaction))

ggsave(dge_n_v_fdr_scatter, filename = here(plot_dir, sprintf("DGE_interaction_%s_gene_v_signif_scatter.png", interaction)), width = 10)

#### logFC heatmaps ####
source(here("code", "13_compile_DGE", "logFC_heatmap.R"))

## top DGEs
topDEGs <- DE_data |>
  group_by(cluster) |>
  filter(vlmf_adj.P.Val < 0.05) |>
  arrange(vlmf_adj.P.Val) |>
  dplyr::slice(1:5) |>
  pull(gene_name) |>
  unique()

length(topDEGs)

logFC_Heatmap(data = DE_data, gene_list = topDEGs, title = sprintf("topDEGs_interaction_%s", interaction), datatype = datatype, h = 6)

## Risk gene heatmap
logFC_Heatmap(data = DE_data, gene_list = AD_risk$symbol, title = sprintf("topDEGs_interaction_%s_ADrisk", interaction), datatype = datatype, h = 6)

if(datatype == "sn_fine"){
  
  ## Age interaction signal in Excit.L2_5.2
  top_ct_DEGs <- DE_data |>
    filter(vlmf_adj.P.Val < 0.05, cluster == "Excit.L2_5.2") |>
    group_by(cluster, vlmf_logFC > 0) |>
    arrange(vlmf_adj.P.Val) |>
    slice(1:25) |>
    pull(gene_name) |>
    unique()
  
  logFC_Heatmap(data = DE_data |>
                  filter(grepl("Excit", cluster)), 
                gene_list = top_ct_DEGs, 
                title = sprintf("topDEGs_interaction_%s_Excit.L2_5.2", interaction), 
                cluster_col = TRUE,
                datatype = datatype)
  
  logFC_Heatmap(data = DE_data |>
                  filter(grepl("Excit", cluster)), 
                gene_list = AD_risk$symbol, 
                title = "ADrisk_Excit")
}


#### Interaction Plots ####
source(here("code","utils","plot_DEG_express.R"))
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)


#### Load pseudobulk data ####
if(datatype == "Visium"){
  pb_fn <- here("processed-data", "09_pseudoBulkDGE_Visium", "01_pseudobulk_data_Visium", "spe_pseudo_DGE.RDS")
  cluster_var <- "SpD"
} else if(datatype == "sn_broad"){
  pb_fn <- here("processed-data", "08_pseudoBulkDGE_sn", "01_pseudobulk_data_sn","sce_pseudo_DGE-cell_type_broad.RDS")
  cluster_var <- "cell_type_broad"
}else if(datatype == "sn_fine"){
  pb_fn <- here("processed-data", "08_pseudoBulkDGE_sn", "01_pseudobulk_data_sn","sce_pseudo_DGE-cell_type_anno.RDS")
  cluster_var <- "cell_type_anno"
} else {
  stop("non-valid datatype")
}

message(Sys.time(), sprintf(" - Datatype = %s, loading '%s'", datatype, basename(pb_fn)))
sce_pb <- readRDS(pb_fn)
dim(sce_pb)

## Drop sample Br1289
sce_pb <- sce_pb[,sce_pb$BrNum != 'Br1289']

table(sce_pb$registration_variable)

rownames(sce_pb) <- rowData(sce_pb)$gene_name

## check model
head(model.matrix(~0 + APOE_carrier + Age + Sex + Anc_Afr + pseudo_expr_chrM_ratio, colData(sce_pb)))

cluster_levels <- levels(sce_pb[[cluster_var]])
cluster_levels2 <- as.character(unique(DE_data$cluster))

cluster_levels <- intersect(cluster_levels, cluster_levels2)

all(cluster_levels %in% sce_pb[[cluster_var]])

dge_n_v_fdr |> filter(n_signif > 0)

table(sce_pb[[cluster_var]])

table(sce_pb[[cluster_var]], sce_pb$APOE_carrier)

summary(sce_pb[, sce_pb$cell_type_anno == "Excit.L2_5.2"]$ncells)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 11.00   15.00   30.00   36.64   50.00   95.00 

colData(sce_pb)[sce_pb$cell_type_anno == "Excit.L2_5.2",c("BrNum","ncells")]

pdf(here(plot_dir, sprintf("DEG_interaction_%s_boxplots_carrier_%s.pdf", interaction, datatype)))
plot_DEG_express_top(sce = sce_pb,
                 stats = DE_data,
                 clus = "Excit.L2_5.2",
                 # n_genes = 10,
                 pval_col = "vlmf_adj.P.Val",
                 fc_col = "vlmf_logFC",
                 gene_col = "gene_name",
                 cluster_col = cluster_var,
                 category_col = "APOE_carrier",
                 mod = ~0 + APOE_carrier + Age + Sex + Ancestry + pseudo_expr_chrM_ratio,
                 color_pal = APOE_carrier_colors,
                 plot_points = TRUE,
                 ncol = 2,
                 cleanY_P = 3)
dev.off()

rps_genes <- DE_interaction_data |> 
  filter(vlmf_adj.P.Val < 0.05) |> 
  select(cluster, gene_name, vlmf_P.Value, vlmf_adj.P.Val, vlmf_logFC) |> 
  filter(grepl("^RPS" , gene_name)) |> 
  dplyr::slice(1:6) |> 
  pull(gene_name)

interaction_plot_RPS_genes <- plot_DEG_express_interaction(
  sce = sce_pb,
  stats = DE_data,
  gene = rps_genes,
  cluster_col = "cell_type_anno",
  clus = "Excit.L2_5.2",
  gene_col = "gene_name"
)

ggsave(interaction_plot_RPS_genes, filename = here(plot_dir, "DEG_scatter_interaction_Age_sn_fine_Excit.L2_5.2_RPSgenes.png"))

## risk gene "MT-ND1"
interaction_plot_MTND1 <- plot_DEG_express_interaction(
  sce = sce_pb,
  stats = DE_data,
  gene = "MT-ND1",
  cluster_col = "cell_type_anno",
  clus = "Excit.L2_5.2",
  gene_col = "gene_name"
)

ggsave(interaction_plot_MTND1, filename = here(plot_dir, "DEG_scatter_interaction_Age_sn_fine_Excit.L2_5.2_MT-ND1.png"), height = 5, width = 5)

top50_degs <- DE_data |> filter(clus == "Excit.L2_5.2", vlmf_adj.P.Val < 0.05) |> arrange(vlmf_adj.P.Val) |> dplyr::slice(1:50)

top50_interaction_plots <- map(top50_degs$gene_name, 
                               ~plot_DEG_express_interaction(
                                 sce = sce_pb,
                                 stats = DE_data,
                                 gene = .x,
                                 cluster_col = "cell_type_anno",
                                 clus = "Excit.L2_5.2",
                                 gene_col = "gene_name"
                               ))

pdf(here(plot_dir, "DEG_scatter_interaction_Age_sn_fine_Excit.L2_5.2_top50.pdf"), height = 5, width = 5)
print(top50_interaction_plots)
dev.off()

# 
# slurmjobs::job_loop(loops = list(datatype = c("sn_broad","sn_fine","Visium")),
#                     create_shell = TRUE,
#                     name = "13_interaction_plots",
#                     create_script = FALSE)


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

