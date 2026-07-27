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
load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE)
load(here("processed-data", "SpD_colors.Rdata"), verbose = TRUE)

cell_type_colors$broad <- cell_type_colors$broad[names(cell_type_colors$broad) != "Other"]

## load DE data ##
data_types <- c("sn_broad", "sn_fine", "Visium")
names(data_types) <- data_types

DE_data_fn <- map(data_types, ~here("processed-data", "13_compile_DGE", "01_compile_DGE", .x, sprintf("DGE_results_carrier_%s.Rds", .x)))
map_lgl(DE_data_fn, file.exists)

DE_data <- map(DE_data_fn, readRDS)

DE_data[["Visium"]]$cluster <- factor(DE_data[["Visium"]]$cluster, levels = names(SpD_colors))
DE_data[["sn_broad"]]$cluster <- factor(DE_data[["sn_broad"]]$cluster, levels = names(cell_type_colors$broad))

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

## single cell broad vs. fine

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

#### combined bar plots ####

DEG_count <- map_dfr(DEGs_signif, ~.x |>
        count(data_type, cluster)) |>
    mutate(data_type = factor(data_type, levels = c("Visium", "sn_broad", "sn_fine")))

DEG_count_data_type_bar <- DEG_count |> 
    ggplot(aes(x = cluster, y = n, fill = cluster)) +
    geom_col() +
    facet_wrap(~data_type, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = c(SpD_colors, cell_type_colors$broad, cell_type_colors$anno)) +
    labs( y = "n DEGs") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = "None")

ggsave(DEG_count_data_type_bar, filename = here(plot_dir, "DEG_count_data_type_bar.png"))

ggsave(DEG_count_data_type_bar, filename = here(plot_dir, "DEG_count_data_type_bar_small.png"), height = 3, width = 6)

library(ggbreak)

DEG_count_data_type_bar_vertical <- DEG_count |> 
    ggplot(aes(y = cluster, x = n, fill = cluster)) +
    geom_col() +
    geom_text(aes(label = ifelse(cluster == "Oligo.3", "", n)), hjust=-0.07) +
    facet_wrap(~data_type, scales = "free_y", space = "free_y", nrow = 1) +
    scale_fill_manual(values = c(SpD_colors, cell_type_colors$broad, cell_type_colors$anno)) +
    scale_x_continuous(
        limits=c(0, 305),
        # expand=expansion(mult=c(0, 0.15)),
        oob=scales::squish
    ) +
    labs(x = "n DEGs") +
    theme_bw()+
    theme(legend.position = "None")

ggsave(DEG_count_data_type_bar_vertical, filename = here(plot_dir, "DEG_count_data_type_bar_vertical.png"), width = 4, height =  5.5)



#### contrast data ####

datatype <- "sn_fine"
contrast <- "ancestry"

## contrast data locations
if(contrast == "ancestry"){
    dge_dir = "05_compile_DGE_ancestry"
    
    contrast_1 <- "carrier_AA"
    contrast_2 <- "carrier_EA"
    
} else if(contrast == "Sex"){
    dge_dir = "09_compile_DGE_Sex"
    
    contrast_1 <- "carrier_F"
    contrast_2 <- "carrier_M"
} 

contrast_levels <- c(contrast_1, contrast_2)
names(contrast_levels) <- contrast_levels
 
DE_data_contrast <- map(data_types, function(datatype){
    
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
    
    DE_data_contrast <- readRDS(here("processed-data", "13_compile_DGE", dge_dir, datatype,
                                     sprintf("DGE_results_%s_%s.Rds", contrast, datatype))) |>
        mutate(cluster = factor(cluster, levels = cluster_levels))
})

DEGs_contrast_signif <- map(DE_data_contrast, function(data){
    map(contrast_levels, ~data |> 
                       filter(vlmf_adj.P.Val < 0.05,
                              contrast == .x) |>
                       mutate(DE_class = case_when(vlmf_logFC > 0 ~ "up",
                                                   vlmf_logFC < 0 ~ "down",
                                                   TRUE ~ "None"),
                              DE_class_cluster = paste0(gsub("\\.", "-", cluster), "_",DE_class))
    )
    })



DEGs_contrast_signif_list <- map_depth(DEGs_contrast_signif, 2, get_gene_list)
DEGs_contrast_signif_list_dir <- map_depth(DEGs_contrast_signif, 2, get_gene_list_dir)


pdf(here(plot_dir, "DEG_Upset_contrast_ancestry_Oligo.3.pdf"))
print(upset(fromList(list(All_up = DEGs_signif_list_dir$sn_fine$`Oligo-3_up`,
                          All_down = DEGs_signif_list_dir$sn_fine$`Oligo-3_down`,
                          AA_up = DEGs_contrast_signif_list_dir$sn_fine$carrier_AA$`Oligo-3_up`,
                          AA_down = DEGs_contrast_signif_list_dir$sn_fine$carrier_AA$`Oligo-3_down`,
                          EA_up = DEGs_contrast_signif_list_dir$sn_fine$carrier_EA$`Oligo-3_up`,
                          EA_down = DEGs_contrast_signif_list_dir$sn_fine$carrier_EA$`Oligo-3_down`)
),
order.by = "freq",
nsets = 6
)
)
dev.off()


#### Venn Diagram ####
library(VennDiagram)

names(DEGs_contrast_signif_list$sn_fine)

DEGs_signif_list$sn_fine$`Oligo.3`
DEGs_contrast_signif_list$sn_fine$carrier_AA$`Oligo.3`
DEGs_contrast_signif_list$sn_fine$carrier_EA$`Oligo.3`


contrast_cols <- c(ALL = "white", ancestry_colors)
# contrast_cols <- contrast_cols[order(names(contrast_cols))]

venn.diagram(
    x = list(DEGs_contrast_signif_list$sn_fine$carrier_AA$`Oligo.3`,
             DEGs_contrast_signif_list$sn_fine$carrier_EA$`Oligo.3`),
    category.names = c("AA " , "EA"),
    fill = ancestry_colors,
    filename = here(plot_dir, "contrast_venn_ancestry_only_Oligo.3.tiff"),
    output=TRUE
)

venn.diagram(
    x = list(DEGs_signif_list$sn_fine$`Oligo.3`,
             DEGs_contrast_signif_list$sn_fine$carrier_AA$`Oligo.3`,
             DEGs_contrast_signif_list$sn_fine$carrier_EA$`Oligo.3`),
    category.names = c("ALL" , "AA " , "EA"),
    fill = contrast_cols,
    filename = here(plot_dir, "contrast_venn_ancestry_Oligo.3.tiff"),
    output=TRUE
)

venn.diagram(
    x = list(DEGs_signif_list_dir$sn_fine$`Oligo-3_up`,
             DEGs_contrast_signif_list_dir$sn_fine$carrier_AA$`Oligo-3_up`,
             DEGs_contrast_signif_list_dir$sn_fine$carrier_EA$`Oligo-3_up`),
    category.names = c("ALL" , "AA " , "EA"),
    fill = contrast_cols,
    filename = here(plot_dir, "contrast_venn_ancestry_Oligo.3_UP.tiff"),
    output=TRUE
)

venn.diagram(
    x = list(DEGs_signif_list_dir$sn_fine$`Oligo-3_down`,
             DEGs_contrast_signif_list_dir$sn_fine$carrier_AA$`Oligo-3_down`,
             DEGs_contrast_signif_list_dir$sn_fine$carrier_EA$`Oligo-3_down`),
    category.names = c("ALL" , "AA " , "EA"),
    fill = contrast_cols,
    filename = here(plot_dir, "contrast_venn_ancestry_Oligo.3_DOWN.tiff"),
    output=TRUE
)


