## Louise Huuki-Myers July 2026
## Check total sets of 

#### set up ####
library("here")
library("tidyverse")
library("qs2")
library("sessioninfo")
library("patchwork")
library("ComplexHeatmap")
library("circlize")

data_dir <- here("processed-data", "13_compile_DGE", "19_validate_summary")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "13_compile_DGE", "19_validate_summary")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## colors
load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE)
load(here("processed-data", "SpX_colors.Rdata"), verbose = TRUE)
load(here("processed-data", "project_colors.Rdata"))
# APOE_carrier_colors
# E2+       E4+
#"#51B8B1" "#D97D59"

## CHECK: logFC_Heatmap() looks up `cluster_levels` from the global env for
## row ordering (it's not a function argument) - loaded per-datatype below,
## right before each call, so it's always the current resolution's levels.

## fill_stat-adaptable logFC/t-stat heatmap helper (logFC_Heatmap / logFC_Heatmap_contrast)
source(here("code", "13_compile_DGE", "logFC_heatmap.R")) 


#### Validation data ####
validation_data_types <- c("Xenium_cell_type_anno", "Xenium_cell_type_anno_SpX", "Xenium_Oligo.3_Astro")

load_DE_valid_data <- function(datatype){
    DE_data_fn <- here("processed-data", "13_compile_DGE", "01_compile_DGE", datatype, sprintf("DGE_results_carrier_%s_wSN.Rds", datatype))
    stopifnot("DE validation data file not found" = file.exists(DE_data_fn))
    message("loading ", datatype, ": ", DE_data_fn)
    readRDS(DE_data_fn)
}

data_type_lookup <- tribble(
    ~data_type,                    ~data_type_short,
    "Xenium_cell_type_anno",       "cell_type",
    "Xenium_cell_type_anno_SpX",   "SpX",
    "Xenium_Oligo.3_Astro",        "Oligo.3_Nbr"
)


DE_valid_data <- validation_data_types |>
    map(load_DE_valid_data) |>
    list_rbind() |>
    mutate(cluster = ifelse(data_type == "Xenium_cell_type_anno_SpX", cluster_SpX, as.character(cluster))) |>
    left_join(data_type_lookup)

all(anchor_genes %in% DE_valid_data$gene_name)

message(nrow(DE_valid_data), " rows across ", n_distinct(DE_valid_data$data_type), " validation data types")
count(DE_valid_data, data_type)

# quick look at which clusters actually validate, per resolution
DE_valid_data |>
    filter(validate) |> 
    count(data_type_short, cluster)

#### Validation check ####

DE_valid_data |> filter(cell_type_anno == "Oligo.3", gene_name == "OPALIN") |> count(validate)

DE_valid_data |> 
    filter(cell_type_anno == "Oligo.3", gene_name == "OPALIN", validate) |> 
    select(data_type, cluster, cell_type_anno, gene_name)

#### Validation summary ####    
valid_genes_summary <- DE_valid_data |>
    filter(validate) |>
    group_by(data_type_short, gene_name, cell_type_anno) |>
    summarise(n = n(), 
              enviro = paste0(unique(cluster), collapse = ", "),
              reg = case_when(all(vlmf_sn_t < 0) ~ "down",
                              all(vlmf_sn_t > 0) ~ "up",
                              TRUE ~ "ERROR")) |>
    arrange(-n)  |>
    mutate(enviro2=ifelse(n>1, "multi", enviro))

valid_genes_summary |> 
    filter(cell_type_anno == "Oligo.3") |> 
    arrange(gene_name) |> 
    print(n = 53) 

DE_valid_data |> filter(data_type == "Xenium_SpX") |> distinct(cluster)


DE_validated <- DE_valid_data |>
    filter(validate)

write_csv(DE_valid_data |> filter(cell_type_anno == "Oligo.3"), file = here(data_dir, "DGE_Xenium_validate_check_Oligo.3.csv"))

write_csv(DE_validated, file = here(data_dir, "DGE_Xenium_validated_All.csv"))
write_csv(valid_genes_summary, file = here(data_dir, "DGE_Xenium_validated_summary.csv"))


#### bar plot ####


DE_valid_data |> 
    filter(cell_type_anno == "Oligo.3", validate) |>
    ggplot(aes(x = data_type_short)) +
    geom_bar(aes(fill = cluster))

valid_genes_summary |>
    filter(cell_type_anno == "Oligo.3") |>
    mutate(enviro2 = ifelse(n > 1, "multi", enviro)) |>
    ggplot(aes(x = data_type_short)) +
    geom_bar(aes(fill = enviro2))

# valid_genes_summary |>
#     filter(cell_type_anno=="Oligo.3") |> 
#     write_csv(here(data_dir, "valid_genes_summary_test.csv"))


valid_genes_summary_filter <- valid_genes_summary |>
    filter(cell_type_anno=="Oligo.3") |>
    ungroup() |>
    mutate(enviro2=fct_relevel(factor(enviro2), "multi"),
           gene_name_reg = paste(gene_name, ifelse(reg == "up", '^', "")))    

# extract ordered SpX suffixes from color names
spx_order <- paste0("Oligo.3_", str_remove(names(SpX_colors) ,"~.*")) 
SpX_colors2 <- SpX_colors
names(SpX_colors2) <- spx_order

spx_levels <- spx_order[spx_order %in% levels(valid_genes_summary_filter$enviro2)]

# full level order: multi first, then Astro enviros, then SpX order
astro_levels <- levels(valid_genes_summary_filter$enviro2) |>
    keep(~str_detect(.x, "APOE|nn_Astro|nnA"))

full_levels <- c("multi", astro_levels, spx_levels, "Oligo.3")

valid_genes_summary_filter <- valid_genes_summary_filter |>
    mutate(enviro2=fct_relevel(enviro2, full_levels))

nnA_colors <- c(
    "APOE_high_nnA_Astro.2" = "#32CD32",
    "APOE_low_nnA_Astro.1"  = "#228B22",
    "APOE_low_nnA_Astro.3"  = "#808000"
)

enviro2_colors <- c(multi = "#2D2D2D", nnA_colors, cell_type_colors$anno, SpX_colors2)[levels(valid_genes_summary_filter$enviro2)]

valid_genes_summary_bar <- valid_genes_summary_filter |>
    arrange(data_type_short, enviro2) |>
    group_by(data_type_short) |>
    mutate(y_pos=row_number() - 0.5) |>
    ungroup() |>
    ggplot(aes(x=data_type_short)) +
    geom_bar(aes(fill=enviro2)) +
    # geom_text(aes(y=y_pos, label=gene_name_reg), size=3, hjust=0.5, color = "white") +
    scale_fill_manual(values = enviro2_colors) +
    labs(x = "Xenium validation context", y = "n validated DEGs") +
    theme_bw() 

ggsave(valid_genes_summary_bar, filename = here(plot_dir, "valid_genes_summary_bar_Oligo.3.png"), width = 5, height = 6)


valid_genes_summary_bar_text <- valid_genes_summary_filter |>
    arrange(data_type_short, enviro2) |>
    group_by(data_type_short) |>
    mutate(y_pos=row_number() - 0.5) |>
    ungroup() |>
    ggplot(aes(x=data_type_short)) +
    geom_bar(aes(fill=enviro2)) +
    geom_text(aes(y=y_pos, label=gene_name_reg), size=3, hjust=0.5, color = "white") +
    scale_fill_manual(values = enviro2_colors) +
    labs(x = "Xenium validation context", y = "n validated DEGs") +
    theme_bw() 

ggsave(valid_genes_summary_bar_text, filename = here(plot_dir, "valid_genes_summary_bar_text_Oligo.3.png"), width = 5, height = 5)
ggsave(valid_genes_summary_bar_text, filename = here(plot_dir, "valid_genes_summary_bar_text_Oligo.3.pdf"), width = 5, height = 5)


valid_genes_summary_bar_text_gap <- valid_genes_summary_filter |>
    arrange(data_type_short, enviro2) |>
    group_by(data_type_short) |>
    mutate(y_pos=row_number() - 0.5) |>
    ungroup() |>
    ggplot(aes(x=data_type_short)) +
    geom_bar(aes(fill=enviro2), width=0.5) +
    geom_text(aes(y=y_pos, label=gene_name_reg), size=3, hjust=0.5, color = "white") +
    scale_fill_manual(values = enviro2_colors) +
    scale_x_discrete(expand=expansion(mult=0.2)) +
    labs(x = "Xenium validation context", y = "n validated DEGs") +
    theme_bw() +
    theme(legend.position = "None")

ggsave(valid_genes_summary_bar_text_gap, filename = here(plot_dir, "valid_genes_summary_bar_text_Oligo.3_gap.png"), width = 5, height = 5)
ggsave(valid_genes_summary_bar_text_gap, filename = here(plot_dir, "valid_genes_summary_bar_text_Oligo.3_gap.pdf"), width = 5, height = 5)


label_data <- valid_genes_summary_filter |>
    arrange(data_type_short, desc(enviro2)) |>
    group_by(data_type_short) |>
    mutate(y_pos=row_number() - 0.5) |>
    ungroup() |>
    group_by(data_type_short, enviro2) |>
    summarise(y_pos=mean(y_pos), .groups="drop")


valid_genes_summary_bar_text_enviro <- valid_genes_summary_filter |>
    arrange(data_type_short, enviro2) |>
    group_by(data_type_short) |>
    mutate(y_pos=row_number() - 0.5) |>
    ungroup() |>
    ggplot(aes(x=data_type_short)) +
    geom_bar(aes(fill=enviro2)) +
    geom_text(data=label_data, aes(y=y_pos, label=enviro2),
              size=2.5, color="white", fontface="bold",) +
    scale_fill_manual(values=enviro2_colors) +
    theme_bw() +
    theme(legend.position = "None")

ggsave(valid_genes_summary_bar_text_enviro, filename = here(plot_dir, "valid_genes_summary_bar_text_enviro_Oligo.3.png"), width = 3, height = 4)


## facet by datatype

DE_valid_data |> 
    filter(cell_type_anno == "Oligo.3", validate) |>
    ggplot(aes(x = cluster)) +
    geom_bar(aes(fill=cluster)) +
    facet_wrap(~data_type, scales = "free_x", space = "free_x")  +
    theme_bw() +
    theme(axis.text.x=element_text(angle=45, hjust=1))

#### Upset plot ####

library(ComplexUpset)

valid_ct <- DE_validated |>
    distinct(cell_type_anno, data_type_short) |>
    count(cell_type_anno) |>
    pull(cell_type_anno)


upset_data <- DE_validated |>
    filter(cell_type_anno == "Oligo.3", validate==TRUE) |>
    distinct(gene_name, data_type_short) |>
    mutate(present=TRUE) |>
    pivot_wider(names_from=data_type_short, values_from=present, values_fill=FALSE)

data_type_short_cols <- upset_data |> select(-gene_name) |> colnames()

datatype_upset <- upset(upset_data, intersect=data_type_short_cols,
                        base_annotations=list(
                            "Intersection size"=intersection_size()
                        ),
                        min_size=1
) + 
    labs(title = "Oligo.3")

ggsave(datatype_upset, filename = here(plot_dir, "valid_datatype_upset_Oligo.3.png"))


