## Louise Huuki-Myers, July 2025
## Compile and plot all DGE data

#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")
library("UpSetR")

# data_dir <- here("processed-data", "13_compile_DGE", "03_DEG_upset")
# if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "13_compile_DGE", "03_DEG_upset")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

load(here("processed-data", "project_colors.Rdata"))

## load DE data ##
data_types <- c("sn_broad", "sn_fine", "Visium")

DE_data_fn <- map(data_types, ~here("processed-data", "13_compile_DGE", "01_compile_DGE", .x, sprintf("DGE_results_carrier_%s.Rds", .x)))
names(DE_data_fn) <- data_types
map_lgl(DE_data_fn, file.exists)

DE_data <- map(DE_data_fn, readRDS)

DEGs_signif <- map(DE_data, ~.x |> 
                       filter(vlmf_adj.P.Val < 0.05) |>
                       mutate(DE_class = case_when(vlmf_logFC > 0 ~ "up",
                                                   vlmf_logFC < 0 ~ "down",
                                                   TRUE ~ "None"),
                              DE_class_cluster = paste0(gsub("\\.", "-", cluster), "_",DE_class)))

map(DEGs_signif, ~.x |> count(DE_class_cluster))

## with direction
get_gene_list_dir <- function(DEG_tb){
    map(rafalib::splitit(DEG_tb$DE_class_cluster), ~DEG_tb$gene_id[.x])
}

## w/o direction
get_gene_list <- function(DEG_tb){
    map(rafalib::splitit(DEG_tb$cluster), ~DEG_tb$gene_id[.x])
}

# get_gene_list(DEGs_signif$sn_broad)

DEGs_signif_list <- map(DEGs_signif, get_gene_list)
DEGs_signif_list_dir <- map(DEGs_signif, get_gene_list_dir)

map2(DEGs_signif_list, names(DEGs_signif_list_dir), function(de_list, data_type){
    pdf(here(plot_dir, sprintf("DEG_Upset_%s.pdf", data_type)))
    print(upset(fromList(de_list), nsets = length(de_list)))
    dev.off()
})


## visium vs. single cell broad

pdf(here(plot_dir, "DEG_Upset_Visium-sn_broad.pdf"))
print(upset(fromList(c(DEGs_signif_list$sn_broad, DEGs_signif_list$Visium)),
            nsets = 10))
dev.off()

## single cell broad vd. fine

pdf(here(plot_dir, "DEG_Upset_sn_broad_v_fine.pdf"))
print(upset(fromList(c(DEGs_signif_list$sn_broad, DEGs_signif_list$sn_fine)),
            nsets = 10))
dev.off()


## All gene sets
DEGs_signif_list_all <- do.call("c", DEGs_signif_list)

pdf(here(plot_dir, "DEG_Upset_ALL.pdf"))
print(upset(fromList(DEGs_signif_list_all), nsets = 10))
dev.off()



common_WMuf_Vasc <- intersect(DEGs_signif_list$Visium$`WM.uf~Sp09D07`, DEGs_signif_list$Visium$`Vasc~Sp09D08`)

DEGs_signif$Visium |> 
    filter(gene_id %in% common_WMuf_Vasc) |> 
    select(cluster, gene_name, vlmf_logFC, vlmf_adj.P.Val)

# cluster       gene_name vlmf_logFC vlmf_adj.P.Val
# <chr>         <chr>          <dbl>          <dbl>
# 1 Vasc~Sp09D08  KLK6          -1.89          0.0152
# 2 Vasc~Sp09D08  CFL2          -0.656         0.0337
# 3 Vasc~Sp09D08  ELOVL1        -1.01          0.0337
# 4 WM.uf~Sp09D07 KLK6          -1.82          0.0279
# 5 WM.uf~Sp09D07 CFL2          -0.952         0.0363
# 6 WM.uf~Sp09D07 ELOVL1        -1.39          0.0438

DEGs_signif$Visium |> 
    filter(gene_name == "KLK6")


