## Compile Xenium ancestry-stratified voomLmFit results and validate against snRNA-seq ancestry findings
## Merges the ancestry-contrast compiling/plotting pattern of 05_compile_DGE_ancestry.R
## with the cross-platform sn-vs-Xenium validation pattern of 01_2_compile_DGE_Xenium.R

#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")
library("ggrepel")
library("GGally")
library("broom")
library("patchwork")

opt <- list()
opt$datatype <- "Xenium"

data_dir <- here("processed-data", "13_compile_DGE", "05_compile_DGE_ancestry", opt$datatype)
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "13_compile_DGE", "05_compile_DGE_ancestry", opt$datatype)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)
AD_risk <- read.csv(here("processed-data", "00_project_prep", "07_OpenTargets_AD_data", "clin_var_genes.csv"))

#### colors and factors ####
## Xenium ancestry clusters are cell_type_anno (fine) labels - same granularity as sn_fine,
## which is what lets us join directly on `cluster` below
load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE)
cluster_colors <- cell_type_colors$anno
cluster_levels <- names(cell_type_colors$anno)

#### voomLmFit data - Xenium ancestry ####
message(Sys.time(), " - load Xenium ancestry VLMF data")

vlmf_fn <- list.files(here("processed-data", "12_voomLmFit", "04_Clusterwise_voomLmFit_ancestry", "vlmf_Xenium"),
                      full.names = TRUE, pattern = ".rds")

names(vlmf_fn) <- map_chr(vlmf_fn, ~gsub("voomLmFit_ancestry_Xenium_|.rds", "", basename(.x)))

## read data - each element is a named list of contrast-level topTables for that cluster
vlmf_data <- map(vlmf_fn, readRDS)

## what contrasts actually came out of the voomLmFit step (should be carrier_AA/carrier_EA)
# message("Contrasts found in Xenium ancestry results: ", paste(contrast_names, collapse = ", "))

vlmf_data_tb <- map_dfr(vlmf_data, ~do.call("rbind", .x)) |>
    dplyr::rename(vlmf_logFC = logFC,
                  vlmf_AveExpr = AveExpr,
                  vlmf_t = t,
                  vlmf_P.Value = P.Value,
                  vlmf_adj.P.Val = adj.P.Val,
                  vlmf_B = B) |>
    as_tibble()

## flag + drop any cluster names that don't match the known color/level scheme, rather than
## silently turning them to NA (or hard-stopping the whole compile over one bad cluster)
unmatched_clusters <- setdiff(unique(vlmf_data_tb$cluster), cluster_levels)
if (length(unmatched_clusters) > 0) {
    message("Note: dropping cluster(s) not found in cell_type_colors$anno: ", paste(unmatched_clusters, collapse = ", "))
}

vlmf_data_tb <- vlmf_data_tb |>
    filter(cluster %in% cluster_levels) |>
    mutate(cluster = factor(cluster, levels = cluster_levels))

levels(vlmf_data_tb$cluster)

vlmf_data_tb |> count(cluster, contrast)
vlmf_data_tb |> filter(vlmf_adj.P.Val < 0.05) |> count(cluster, contrast)

#### summary + bar plots, one per contrast (as in 05_compile_DGE_ancestry.R) ####
vlmf_model_summary <- vlmf_data_tb |>
    group_by(cluster, contrast) |>
    summarize(n_genes = n(),
              n_FDR05 = sum(vlmf_adj.P.Val < 0.05),
              nUP = sum(vlmf_adj.P.Val < 0.05 & vlmf_logFC > 0),
              nDown = sum(vlmf_adj.P.Val < 0.05 & vlmf_logFC < 0),
              .groups = "drop")

vlmf_model_summary

vlmf_model_summary_bar <- vlmf_model_summary |>
    ggplot(aes(x = cluster, y = n_FDR05, fill = cluster)) +
    geom_col() +
    geom_text(aes(label = n_FDR05), vjust = -.5) +
    scale_fill_manual(values = cluster_colors) +
    facet_wrap(~contrast, ncol = 1) +
    theme_bw() +
    labs(title = "voomLmFit - Xenium - carrier by ancestry", subtitle = "FDR < 0.05") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = "None")

ggsave(vlmf_model_summary_bar, filename = here(plot_dir, "ancestry_Xenium_vlmf_model_summary_bar.png"), height = 10, width = 10)

vlmf_model_summary_bar_reg <- vlmf_model_summary |>
    select(-n_FDR05) |>
    mutate(nDown = -1 * nDown) |>
    pivot_longer(!c(cluster, contrast), names_to = "reg", values_to = "n_genes") |>
    ggplot(aes(x = cluster, y = n_genes, fill = reg)) +
    geom_col() +
    geom_text(aes(label = abs(n_genes))) +
    facet_wrap(~contrast, ncol = 1) +
    theme_bw() +
    labs(title = "voomLmFit - Xenium - carrier by ancestry", subtitle = "FDR < 0.05") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(vlmf_model_summary_bar_reg, filename = here(plot_dir, "ancestry_reg_Xenium_vlmf_model_summary_bar.png"), width = 10)

write.csv(vlmf_model_summary, file = here(data_dir, "vlmf_ancestry_model_summary_Xenium.csv"), row.names = FALSE)

#### vlmf volcano plots ####
custom_volcano <- function(data, FDR_cut = 0.05, model_name){

    signif_colors <- c("purple", "blue", "red")
    names(signif_colors) <- c("both", paste("FDR<", FDR_cut), "abs(logFC)>1")

    volcano <- data |>
        mutate(DE_class = case_when(vlmf_adj.P.Val < FDR_cut & abs(vlmf_logFC) > 1 ~ "both",
                                     vlmf_adj.P.Val < FDR_cut ~ paste("FDR<", FDR_cut),
                                     abs(vlmf_logFC) > 1 ~ "abs(logFC)>1",
                                     TRUE ~ "None")) |>
        ggplot(aes(x = vlmf_logFC, y = -log10(vlmf_P.Value), color = DE_class)) +
        geom_point(alpha = 0.5, size = 0.5) +
        scale_color_manual(values = signif_colors) +
        facet_grid(contrast ~ cluster) +
        theme_bw() +
        theme(legend.position = "bottom") +
        labs(title = model_name)

    if (nrow(data) > 1000) {
        volcano <- volcano + geom_text_repel(aes(label = ifelse(DE_class != "None", gene_name, "")), size = 1.5)
    } else {
        volcano <- volcano + geom_text_repel(aes(label = gene_name), size = 1.5)
    }
    ggsave(volcano, filename = here(plot_dir, sprintf("Volcano_plot_%s.png", model_name)), height = 10, width = 14)
}

custom_volcano(data = vlmf_data_tb, model_name = "Xenium-ancestry")
# custom_volcano(data = vlmf_data_tb |> filter(gene_name %in% AD_risk$symbol), model_name = "Xenium-ancestry-risk")

#### save compiled Xenium ancestry data ####
saveRDS(vlmf_data_tb, file = here(data_dir, "DGE_results_ancestry_Xenium.Rds"))
write.csv(vlmf_data_tb, file = here(data_dir, "DGE_results_ancestry_Xenium.csv"), row.names = FALSE)

# vlmf_data_tb <- readRDS(here(data_dir, "DGE_results_ancestry_Xenium.Rds"))

#### Add snRNA-seq ancestry results for validation ####
message(Sys.time(), " - load sn_fine ancestry results for validation")

## requires 05_compile_DGE_ancestry.R --datatype sn_fine to have been run already
sn_ancestry_fn <- here("processed-data", "13_compile_DGE", "05_compile_DGE_ancestry", "sn_fine", "DGE_results_ancestry_sn_fine.Rds")
if (!file.exists(sn_ancestry_fn)) {
    stop(sprintf("sn_fine ancestry results not found at '%s' - run 05_compile_DGE_ancestry.R --datatype sn_fine first", sn_ancestry_fn))
}

sn_ancestry_data <- readRDS(sn_ancestry_fn) |>
    select(cluster, gene_id, gene_name, contrast, starts_with("vlmf")) |>
    rename_with(~str_replace(.x, "vlmf_", "vlmf_sn_"), starts_with("vlmf_")) |>
    mutate(cluster = as.character(cluster))

## 05_compile_DGE_ancestry.R only carries carrier_AA/carrier_EA 
shared_contrasts <- intersect(unique(as.character(vlmf_data_tb$contrast)), unique(sn_ancestry_data$contrast))

message("Shared contrasts between Xenium and sn_fine ancestry results: ", paste(shared_contrasts, collapse = ", "))

vlmf_data_tb_xenium <- vlmf_data_tb |>
    filter(contrast %in% shared_contrasts) |>
    mutate(cluster = as.character(cluster), contrast = as.character(contrast)) |>
    rename_with(~str_replace(.x, "vlmf_", "vlmf_xenium_"), starts_with("vlmf_")) |>
    left_join(sn_ancestry_data, by = c("cluster", "gene_id", "gene_name", "contrast")) |>
    mutate(cluster = factor(cluster, levels = cluster_levels),
           signif_xenium = vlmf_xenium_P.Value < 0.1,
           signif_sn = vlmf_sn_adj.P.Val < 0.05,
           signif_both = signif_xenium & signif_sn,
           dir_match = (vlmf_xenium_t > 0) == (vlmf_sn_t > 0),
           validate = signif_both & dir_match)

## save data
saveRDS(vlmf_data_tb_xenium, file = here(data_dir, "DGE_results_ancestry_Xenium_wSN.Rds"))
write.csv(vlmf_data_tb_xenium, file = here(data_dir, "DGE_results_ancestry_Xenium_wSN.csv"), row.names = FALSE)

vlmf_data_tb_xenium |>
    filter(validate) |>
    select(cluster, contrast, gene_id, gene_name, vlmf_sn_logFC, vlmf_sn_adj.P.Val, vlmf_xenium_logFC, vlmf_xenium_P.Value, dir_match) |>
    print(n = 5)

vlmf_data_tb_xenium |>
    filter(gene_name == "NTRK3", cluster == 'Oligo.3') |>
    select(cluster, contrast, gene_id, gene_name, signif_sn, vlmf_sn_logFC, vlmf_sn_adj.P.Val, vlmf_xenium_logFC, vlmf_xenium_P.Value, dir_match)

#### Summarize validation experiments (per cluster x contrast) ####
message(Sys.time(), " - summarize validation experiments")

validation_summary <- vlmf_data_tb_xenium |>
    filter(!is.na(vlmf_sn_P.Value)) |>
    group_by(cluster, contrast) |>
    summarise(n_signif_xenium = sum(signif_xenium),
              n_signif_sn = sum(signif_sn),
              n_signif_both = sum(signif_both),
              n_validate = sum(validate),
              percent_valid = ifelse(n_signif_sn > 0, 100 * n_validate / n_signif_sn, NA_real_),
              .groups = "drop")

validation_summary |> filter(n_signif_sn > 0) |> arrange(-n_validate) |> print(n = 10)

write_csv(validation_summary, file = here(data_dir, "ancestry_Xenium_DEG_validation_summary.csv"))

## validation bar plot
validation_bar <- validation_summary |>
    filter(n_signif_sn > 0) |>
    ggplot(aes(x = cluster, y = n_validate, fill = cluster)) +
    geom_col() +
    geom_text(aes(label = n_validate), vjust = -0.5) +
    scale_fill_manual(values = cluster_colors) +
    facet_wrap(~contrast, ncol = 1) +
    theme_bw() +
    labs(title = "Xenium validation of snRNA-seq ancestry DEGs",
         subtitle = "n genes: signif in sn (FDR<0.05) AND signif + direction-matched in Xenium (p<0.10)") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = "None")

ggsave(validation_bar, filename = here(plot_dir, "ancestry_Xenium_validation_bar.png"), height = 8, width = 10)

#### t-stat correlation, sn vs Xenium, per cluster x contrast ####
message(Sys.time(), " - Xenium vs. snRNA-seq ancestry t-stat correlation")

sn_xenium_cor <- vlmf_data_tb_xenium |>
    filter(!is.na(vlmf_sn_t), !is.na(vlmf_xenium_t)) |>
    group_by(cluster, contrast) |>
    filter(n() >= 10) |>
    summarise(tidy(cor.test(vlmf_sn_t, vlmf_xenium_t, method = "spearman")),
              n_genes = n(),
              .groups = "drop") |>
    select(cluster, contrast, n_genes, cor = estimate, p_value = p.value) |>
    mutate(FDR = p.adjust(p_value, method = "BH"))

write_csv(sn_xenium_cor, file = here(data_dir, "ancestry_Xenium_v_sn_tstat_cor.csv"))

sn_xenium_cor |> arrange(-cor) |> print(n = 50)

#### t-stat scatter plots (Xenium vs sn, colored by validation status) ####
compare_stats_scatter_ancestry <- function(dge_tb, cor_tb, model_name, clus){

    dge_tb <- dge_tb |>
        filter(cluster == clus, !is.na(vlmf_sn_t)) |>
        mutate(cluster = droplevels(cluster),
               DE_class = case_when(validate ~ "validated",
                                     signif_sn ~ "sn only (FDR<0.05)",
                                     signif_xenium ~ "Xenium only (p<0.10)",
                                     TRUE ~ "None"))

    signif_colors <- c(validated = "purple",
                       `sn only (FDR<0.05)` = "blue",
                        `Xenium only (p<0.10)` = "red", None = "grey70")

    scatter <- dge_tb |>
        ggplot(aes(x = vlmf_sn_t, y = vlmf_xenium_t, color = DE_class)) +
        geom_point(alpha = 0.5, size = 0.5) +
        geom_text_repel(aes(label = ifelse(DE_class != "None", gene_name, "")), size = 2) +
        geom_abline(linetype = "dashed") +
        scale_color_manual(values = signif_colors) +
        facet_grid(contrast ~ cluster) +
        theme_bw() +
        theme(legend.position = "bottom") +
        labs(x = "sn t-stat", y = "Xenium t-stat")

    ggsave(scatter, filename = here(plot_dir, sprintf("ancestry_Xenium_v_sn_tstat_scatter_%s.png", model_name)),
           height = 7,
           width = 7)
}

select_clusters <- validation_summary |> filter(n_signif_sn > 0) |> pull(cluster) |> unique()

walk(select_clusters, ~compare_stats_scatter_ancestry(vlmf_data_tb_xenium, clus = .x, model_name = paste0("all_genes_", .x)))

#### correlation heatmap ####
sn_xenium_cor_heatmap <- sn_xenium_cor |>
    ggplot(aes(x = contrast, y = cluster)) +
    geom_tile(aes(fill = cor)) +
    geom_text(aes(label = ifelse(FDR < 0.05, "*", ""))) +
    theme_bw() +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#D6604D", midpoint = 0) +
    labs(title = "Xenium vs sn_fine ancestry t-stat correlation (Spearman)")

ggsave(sn_xenium_cor_heatmap, filename = here(plot_dir, "ancestry_Xenium_v_sn_tstat_cor_heatmap.png"), width = 6, height = 8)


#### contrast scatter ####

validation_genes <- vlmf_data_tb_xenium |> filter(validate) |> pull(gene_name)

valid_contrast <- vlmf_data_tb_xenium |> 
    filter(validate)  |>
    select(gene_name, contrast) 

contrast_t_stats <- vlmf_data_tb_xenium |> 
    filter(gene_name %in% validation_genes, cluster == "Oligo.3") |>
    select(gene_name, contrast, vlmf_xenium_t) |>
    arrange(gene_name) |>
    pivot_wider(names_from = "contrast", values_from = "vlmf_xenium_t") |>
    mutate(datatype = "Xenium") |> rbind(
        vlmf_data_tb_xenium |> 
            filter(gene_name %in% validation_genes, cluster == "Oligo.3") |>
            select(gene_name, contrast, vlmf_sn_t) |>
            arrange(gene_name) |>
            pivot_wider(names_from = "contrast", values_from = "vlmf_sn_t") |>
            mutate(datatype = "sn")
    ) |>
    left_join(valid_contrast) 


names(ancestry_colors) <- paste0("carrier_", names(ancestry_colors))

contrast_scatter <- contrast_t_stats |>
    ggplot(aes(x = carrier_AA, y = carrier_EA, color = contrast)) +
    geom_point(alpha = 0.5, size = 0.5) +
    geom_text_repel(aes(label = gene_name)) +
    theme_bw() +
    geom_vline(xintercept = 0)+
    geom_hline(yintercept = 0) +
    labs(x = "t-stat: AA", y = "t-stat: EA") +
    scale_color_manual(values = ancestry_colors) +
    facet_wrap(~datatype)

ggsave(contrast_scatter, filename = here(plot_dir, "ancestry_valid_scatter.png"), width = 8, height = 4)


# slurmjobs::job_single('05.1_compile_DEG_ancestry_Xenium', create_shell = TRUE, memory = '10G', command = "Rscript 05.1_compile_DEG_ancestry_Xenium.R")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
