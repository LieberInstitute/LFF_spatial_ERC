## November 2024, Louise Huuki-Myers
## Compare preprint SpD annotation ("SpD_preprint") to the updated SVGm k11
## clustering used in 18_compare_clusters.R ("SpD_update")

library("spatialLIBD")
library("tidyverse")
library("bluster")
library("ComplexHeatmap")
library("patchwork")
library("here")
library("sessioninfo")
library("readxl")

#### define dirs ####
data_dir <- here("processed-data", "05_spe_correct_cluster", "18_compare_clusters")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "05_spe_correct_cluster", "18_compare_clusters")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### source plotting function ####
source(here("code", "utils", "my_plot_reduced_dim.R"))

#### specify and load data ####
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC_annotated"))
spe

reducedDimNames(spe)
# [1] "10x_pca"      "10x_tsne"     "10x_umap"     "PCA_p1"       "PCA_p2"       "TSNE"         "UMAP"        
# [8] "HARMONY"      "UMAP.HARMONY" "TSNE.HARMONY"

## preprint annotation, kept under an explicit name for clarity
##
## 19_SpD_update_spe.R is about to be re-run and will overwrite
## spe_ERC_annotated with the k11-based updated annotation (SpDv/SpDv_k11),
## dropping the old BayesSpace_SVGm_k09 + SpD (preprint) columns loaded here.
## Stash those columns to a csv now, while they're still around, keyed on
## `key` so they can be matched back onto the post-overwrite spe. Once
## 19_SpD_update_spe.R has been re-run, flip this to `if (FALSE)` to reload
## the stashed preprint labels instead of expecting spe$SpD to still exist.
if (FALSE) { # Overwriting spe object on Aug 20 , 2026
    preprint_stash <- colData(spe) |>
        as.data.frame() |>
        select(key, BayesSpace_SVGm_k09, SpD)

    write_csv(preprint_stash, file = here(data_dir, "preprint_SpD_stash.csv"))

    spe$SpD_preprint <- spe$SpD
} else {
    preprint_stash <- read_csv(here(data_dir, "preprint_SpD_stash.csv"))

    spe$SpD_preprint <- preprint_stash$SpD[match(spe$key, preprint_stash$key)]
}

table(spe$SpD_preprint)

#### Load updated (SVGm k11) clustering + annotation ####

bayes_cluster_fn <- list.files(here("processed-data", "05_spe_correct_cluster", "05_BayesSpace", "clusters_BayesSpace-SVGm", "clusters_BayesSpace"),
                               recursive = TRUE,
                               full.names = TRUE)

map(bayes_cluster_fn, ~as.Date(file.info(.x)$mtime))

bayes_clusters <- map2(bayes_cluster_fn, map(str_split(bayes_cluster_fn, "/"), 12), function(file, name){
    sp = sprintf("Sp%02d",parse_number(name))
    clus <- read.csv(file) |>
        mutate(cluster = factor(paste0(sp, "D", str_pad(cluster, 2,  pad = "0"))))
    
    colnames(clus)[2] <- name
    return(clus)
})


bayes_clusters_tab <- purrr::reduce(bayes_clusters, dplyr::left_join, by = c("key"))
## fix double brain number in key
bayes_clusters_tab$key <- gsub("(^.*?_Br[0-9]+)_Br[0-9+]+$", "\\1", bayes_clusters_tab$key)
identical(colnames(spe), bayes_clusters_tab$key)
head(bayes_clusters_tab)

pd <- colData(spe) |> as.data.frame()
pd <- pd[, !grepl("BayesSpace_SVGm_", colnames(pd))]

pd <- left_join(pd, bayes_clusters_tab)

colData(spe) <- DataFrame(pd)

## same convention as 18_compare_clusters.R: BayesSpace_SVGm_k11 clusters,
## annotated via the k11 spatial registration summary
anno_table_k11 <- read_xlsx(here(
    "processed-data", "05_spe_correct_cluster", "10_spatial_registration_DLPFC",
    "ERC_SpD_spatial_registration_anno_summary_k11.xlsx"
)) |>
    mutate(SpD = fct_reorder(paste0(anno, "~", cluster), order))

spe$SpD_update <- anno_table_k11$SpD[match(spe$BayesSpace_SVGm_k11, anno_table_k11$cluster)]

table(spe$SpD_update, spe$BayesSpace_SVGm_k11)

#### Match colors ####
## preprint colors, stashed here directly (not in a file) since
## metadata(spe)$SpD_colors will be overwritten with the new V4 palette once
## 19_SpD_update_spe.R re-runs
SpD_colors_simple <- c("Vasc~Sp09D08" = "#E05AD2", #Orchid
                       "L1~Sp09D05" = "#0220DE", #Chrystler Blue
                       "L2.3~Sp09D01" = "#FEAF16", #light orange
                       "LD~Sp09D02" = "#00BCF9", #dark sky blue
                       "Inhib~Sp09D09" = "#C82100", #Engineering red
                       "L5~Sp09D03" = "#16FF32", #lime
                       "L6~Sp09D04" = "#178C6D", #forest green
                       "WM.uf~Sp09D07" = "#E4E1E3", # purpe white
                       "WM~Sp09D06" = "#581009") #brown

## updated (k11) uses the new V4 domain palette, mapped onto the k11 SpD labels
SpD_colors_V4 <- c(
    "Vasc"      = "#E05AD2",
    "L1"        = "#47C281",
    "L2"        = "#41D4EB",
    "L3"        = "#889DF0",
    "L3_Inhib"  = "#B6686F",
    "LD"        = "grey80",
    "L5_LD"     = "#B1C2CE",
    "L5"        = "#0072CE",
    "L6"        = "#0A2E5C",
    "L6a"       = "#24487A",
    "L6b"       = "#05193B",
    "WMuf"      = "#F4A460",
    "WMtz"      = "#E8720C",
    "WM"        = "#F57A00",
    "WMd"       = "#581009",
    "Inhib"     = "#C82100"
)

#' Build a color vector for one resolution's annotation table, keyed by SpD
make_SpD_colors <- function(anno_table, domain_colors = SpD_colors_V4) {
    missing_domains <- setdiff(anno_table$anno, names(domain_colors))
    if (length(missing_domains) > 0) {
        stop("No color defined for domain(s): ", paste(missing_domains, collapse = ", "))
    }
    colors_out <- domain_colors[anno_table$anno]
    names(colors_out) <- anno_table$SpD
    colors_out
}

SpD_colors_update <- make_SpD_colors(anno_table_k11)

#### Spot plots ####
sample_order <- sort(unique(spe$BrNum))

vis_grid_clus(
    spe = spe,
    clustervar = "SpD_update",
    pdf = here(plot_dir, "spe_erc-BayesSpace_SVGm_k11_update.pdf"),
    sort_clust = FALSE,
    point_size = 1.2,
    colors = SpD_colors_update,
    sample_order = sample_order
)

#### Compare preprint vs. update with a Jaccard-style correspondence matrix ####
cluster_tab <- table(spe$SpD_preprint, spe$SpD_update)

cluster_tab_long <- cluster_tab |>
    reshape2::melt() |>
    rename(SpD_preprint = Var1, SpD_update = Var2, n = value)

jacc_mat <- linkClustersMatrix(spe$SpD_preprint, spe$SpD_update)

jacc_mat_long <- jacc_mat |>
    reshape2::melt() |>
    rename(SpD_preprint = Var1, SpD_update = Var2, Jacc = value) |>
    left_join(cluster_tab_long)

write_csv(jacc_mat_long, file = here(data_dir, "preprint_v_update_jacc_mat_long_k11.csv"))

## for each updated (k11) cluster, its single best-matching preprint domain
cluster_match <- jacc_mat_long |>
    group_by(SpD_update) |>
    arrange(-Jacc) |>
    slice(1) |>
    ungroup() |>
    arrange(-Jacc)

sum(cluster_match$n) / ncol(spe)

#### Jacc heatmap, annotated with each clustering's own colors ####
row_ha <- rowAnnotation(
    df = data.frame(SpD_preprint = rownames(jacc_mat)),
    col = list(SpD_preprint = SpD_colors_simple),
    show_legend = FALSE
)

col_ha <- HeatmapAnnotation(
    df = data.frame(SpD_update = colnames(jacc_mat)),
    col = list(SpD_update = SpD_colors_update),
    show_legend = FALSE
)

pdf(here(plot_dir, "jacc_mat_preprint_v_update_k11.pdf"))
print(Heatmap(
    jacc_mat,
    name = "Correspondence",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    right_annotation = row_ha,
    bottom_annotation = col_ha,
    col = c("black", viridisLite::plasma(100)),
    na_col = "black"
))
dev.off()

#### Representative section grid: SpD_preprint (top) vs SpD_update (bottom) ####
rep_sections_tb <- read.csv(here("processed-data", "05_spe_correct_cluster", "22_SpD_clean_plots", "rep_section.csv"))

rep_samples_4 <- rep_sections_tb |>
    filter(rep_section_4) |>
    pull(sample_id)

row_clustervars <- c(SpD_preprint = "SpD_preprint", SpD_update = "SpD_update")
row_colors <- list(SpD_preprint = SpD_colors_simple, SpD_update = SpD_colors_update)

cluster_plots_by_row <- map(names(row_clustervars), function(row_label) {
    clustervar <- row_clustervars[[row_label]]
    row_pal <- row_colors[[row_label]]

    cluster_row_plots <- imap(rep_samples_4, function(s, i) {
        is_last <- i == length(rep_samples_4)

        vis_clus_plot <- vis_clus(
            spe = spe,
            point_size = 1.3,
            colors = row_pal,
            sampleid = s,
            clustervar = clustervar
        ) +
            labs(title = s, y = if (i == 1) row_label else NULL) +
            theme(
                legend.position = if (is_last) "right" else "none", ## only right-most sample keeps its legend
                axis.title.x = element_blank(),
                text = element_text(size = 12),
                plot.title = element_text(hjust = 0.5)
            )

        vis_clus_plot
    })

    Reduce("+", cluster_row_plots) + plot_layout(nrow = 1)
})

cluster_grid <- Reduce("/", cluster_plots_by_row)

# ggsave(cluster_grid, filename = here(plot_dir, "vis_SpD_preprint_v_update_rep_sections.pdf"), width = 18, height = 6)
ggsave(cluster_grid, filename = here(plot_dir, "vis_SpD_preprint_v_update_rep_sections.png"), width = 18, height = 9)

# slurmjobs::job_single('18.1_compare_clusters_offset', create_shell = TRUE, memory = '10G', command = "Rscript 18.1_compare_clusters_offset")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
