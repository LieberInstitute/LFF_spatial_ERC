## Louise Huuki-Myers, June 2025
## Compile and plot all DGE data

#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")
library("ggrepel")

data_dir <- here("processed-data", "13_compile_DGE", "01_compile_DGE")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "13_compile_DGE", "01_compile_DGE")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

load(here("processed-data", "project_colors.Rdata"))
AD_risk <- read.csv(here("processed-data", "00_project_prep", "07_OpenTargets_AD_data", "clin_var_genes.csv")) 
main_mods <- c("carrier", "apoe", "e4e4", "carrier_i", "apoe_i", "e4e4_i")

#### pseudobulkDGE data ####

pseudobulkDGE_fn <- list(sn = here("processed-data", "08_pseudoBulkDGE_sn", "05_pseudoBulkDGE_sn_broad", "sn_pseudoBulkDGE.rds"),
                          Visium = here("processed-data", "09_pseudoBulkDGE_Visium", "03_pseudoBulkDGE_Visium", "Visium_pseudoBulkDGE.rds"))

map_lgl(pseudobulkDGE_fn, file.exists)

## load main mod data
pseudobulkDGE_data <- map(pseudobulkDGE_fn, ~readRDS(.x)[c("carrier", "APOE", "E4E4", "carrier_i", "APOE_i", "E4E4_i")])
## fix names
names(pseudobulkDGE_data$sn) <- main_mods
names(pseudobulkDGE_data$Visium) <- main_mods


pseudobulkDGE_data$sn$carrier$Astro

edit_pseudobulkDEG_data <- function(df, cluster){
    
    df2 <- df |>
        as.data.frame() |>
        mutate(cluster = cluster) |>
        select(gene_id, gene_name, logFC, t, P.Value, adj.P.Val)
    
    names_keep <- c()
    
}
    
    
#### dreamlet data ####
dreamlet_fn <- list(sn = here("processed-data", "10_dreamlet_sn", "05_compile_dreamlet_sn", "dreamlet_sn_TopTables.Rdata"),
                         Visium = here("processed-data", "11_dreamlet_Visium", "05_compile_dreamlet_Visium", "dreamlet_Visium_TopTables.Rdata"))

map_lgl(dreamlet_fn, file.exists)

dreamlet_data <- map(dreamlet_fn, ~get(load(.x)))

names(dreamlet_data$sn)

head(dreamlet_data$sn$carrier)

#### voomLmFit data ####

vlmf_fn <- map(list(sn = here("processed-data", "12_voomLmFit", "01_Clusterwise_voomLmFit", "vlmf_sn_broad"),
                Visium = here("processed-data", "12_voomLmFit", "01_Clusterwise_voomLmFit", "vlmf_Visium")),
               ~list.files(.x, full.names = TRUE, pattern = ".rds"))
names <- map_depth(vlmf_fn, 1, ~gsub("voomLmFit_sn_broad_|voomLmFit_Visium_|.rds", "", basename(.x)))

vlmf_fn <- map2(vlmf_fn, names, ~setNames(.x, .y))

##bug ?
# vlmf_data <- map_depth(vlmf_fn, 1, readRDS)

vlmf_data <- list()
vlmf_data$sn <- map(vlmf_fn[["sn"]], readRDS)
vlmf_data$Visium <- map(vlmf_fn[["Visium"]], readRDS)

vlmf_data <- map_depth(vlmf_data, 1, list_transpose)
names(vlmf_data[["sn"]])
head(vlmf_data$sn$"E2E2_E4E4"$Astro)

vlmf_data_tb <- map_depth(vlmf_data, 2, 
                          ~as_tibble(do.call("rbind", .x)) |>
                              dplyr::rename(vlmf_logFC = logFC,
                                            vlmf_AveExpr = AveExpr,
                                            vlmf_t = t,
                                            vlmf_P.Value = P.Value,
                                            vlmf_adj.P.Val = adj.P.Val,
                                            vlmf_B = B
                                            )
                          )


#### Add other data to vlmc data ####

test <- vlmf_data_tb$sn$anyE2_anyE4 |>
    left_join(dreamlet_data$sn$carrier |>
                  as.data.frame() |> 
                  select(gene_name = ID,
                         cluster = assay,
                         dream_logFC = logFC,
                         dream_AveExpr = AveExpr,
                         dream_t = t,
                         dream_P.Value = P.Value,
                         dream_adj.P.Val = adj.P.Val.cell_type,
                         dream_B = B)) |>
    mutate(DE_class = case_when(vlmf_adj.P.Val < 0.05 & dream_adj.P.Val < 0.2 ~ "signif_both",
                                vlmf_adj.P.Val < 0.05  ~ "signif_vlmf",
                                dream_adj.P.Val < 0.2 ~ "signif_dreamlet",
                                TRUE ~"None"
    ))

test |> count(DE_class)

vlmf_dreamlet_scatter <- test |>
    ggplot(aes(x = vlmf_t, y = dream_t, color = DE_class)) +
    geom_point(alpha = 0.5, size = 0.5) +
    geom_text_repel(aes(label = ifelse(DE_class != "None", gene_name, "")), size = 1.5) +
    geom_abline(linetype = "dashed") +
    scale_color_manual(values = c(signif_both = "purple", signif_vlmf = "blue", signif_dreamlet = "red")) +
    facet_wrap(~cluster) +
    theme_bw()

ggsave(vlmf_dreamlet_scatter, filename = here(plot_dir, "vlmf_dreamlet_scatter_carrier.png"))
