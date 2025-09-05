## Louise Huuki-Myers, Aug 2025
## MOFA gene weight heatmap

#### set up ####
library("tidyverse")
library("ComplexHeatmap")
library("here")
library("sessioninfo")
library("getopt")

# Import command-line parameters
scec <- matrix(
    c("datatype", "d", "1", "character", "Data type"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)

## test 
# opt$datatype = "sn_broad"
# opt$datatype = "sn_fine"

data_dir <- here("processed-data", "14_MOFA", "02_MOFA_gw_heatmaps", opt$datatype)
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "14_MOFA", "02_MOFA_gw_heatmaps", opt$datatype)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## colors
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)
load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE)
load(here("processed-data", "SpD_colors.Rdata"), verbose = TRUE)

## factors 
spd_factor_levels <- names(SpD_colors)

if(opt$datatype == "sn_broad"){
    cell_type_colors$broad <- cell_type_colors$broad[names(cell_type_colors$broad) != "Other"]
    sn_factor_levels <- names(cell_type_colors$broad)
}else if(opt$datatype == "sn_fine"){
    sn_factor_levels <- names(cell_type_colors$anno)
}

cluster_levels <- c(sn_factor_levels, spd_factor_levels)

#### load gene weight data ####

gene_weights <- readRDS(here("processed-data", "14_MOFA", "01_MOFA", opt$datatype, sprintf('MOFA_gene_weights_%s.rds', opt$datatype)))
                        
gene_weights <- map(gene_weights, ~.x |>
                            mutate(ctype = factor(gsub("_Sp", "~Sp", ctype),
                                   levels = cluster_levels))
                    )


map(gene_weights, ~.x |> group_by(ctype) |> slice_max(value))

#### Factor 3 heatmap ####

summary(gene_weights$Factor3$value)

if(opt$datatype == "sn_broad"){
    my_views <- c("Oligo","WM~Sp09D06")
} else if(opt$datatype == "sn_fine"){
    my_views <- c("Oligo.3","Oligo.4","Oligo.5","WM~Sp09D06") #Excit.L5.2 - ugly when added
}

top_gene_weights <- gene_weights$Factor3 |>
    mutate(abs_value = abs(value),
           weight_pos = value > 0,
           datatype = ifelse(grepl("_Sp", ctype), "Visium", "snRNA-seq")) |>
    filter(ctype %in% my_views) |>
    group_by(ctype,weight_pos) |>
    arrange(-abs_value) |>
    dplyr::slice(1:5)


## prep heatmap 
top_gw_value_matrix <- gene_weights$Factor3 |>
    filter(feature %in% top_gene_weights$feature) |>
    select(gene_name, ctype, value) |>
    pivot_wider(names_from = ctype, values_from = value) |>
    column_to_rownames("gene_name") |>
    as.matrix()

col_order <- cluster_levels[cluster_levels %in% colnames(top_gw_value_matrix)]

# row_order <- order(top_gw_value_matrix[,"Oligo"])

row_order <- top_gene_weights |>
    group_by(gene_name) |> 
    summarise(max_value = max(value)) |>
    arrange(-max_value) |>
    pull(gene_name)

top_gw_value_matrix <- top_gw_value_matrix[row_order,col_order]

dim(top_gw_value_matrix)
top_gw_value_matrix[1:5, 1:5]


pdf_width_all <- (nrow(top_gw_value_matrix)/8) + 3
pdf_width_select <- (length(my_views)/8) + 3

pdf_height <- length(row_order)/5

## weights across all clusters
pdf(here(plot_dir, sprintf("MOFA_gene_weight_heatmap_%s_all.pdf", opt$datatype)), width = pdf_width_all, height = pdf_height)
Heatmap(top_gw_value_matrix,
        name = "feature\nweights",
        cluster_rows = FALSE,
        cluster_columns = FALSE)
dev.off()

## select clusters
pdf(here(plot_dir, sprintf("MOFA_gene_weight_heatmap_%s_select.pdf", opt$datatype)), width = pdf_width_select, height = pdf_height)
Heatmap(top_gw_value_matrix[,my_views],
        name = "feature\nweights",
        cluster_rows = FALSE,
        cluster_columns = FALSE)
dev.off()

#### DEG data ####

spd_dge_data <- readRDS(here("processed-data", "13_compile_DGE", "01_compile_DGE", "Visium",
                         sprintf("DGE_results_carrier_%s.Rds", "Visium"))) 

sn_dge_data <- readRDS(here("processed-data", "13_compile_DGE", "01_compile_DGE", opt$datatype,
                         sprintf("DGE_results_carrier_%s.Rds", opt$datatype))) 


dge_data <- bind_rows(spd_dge_data, sn_dge_data)



source(here("code", "13_compile_DGE", "logFC_heatmap.R"))

## all clusters
logFC_Heatmap(dge_data, 
              gene_list = row_order,
              title = "Factor3", 
              h = pdf_height, 
              w = pdf_width_all, 
              cluster_col = FALSE,
              datatype = "MOFA",
              flip = TRUE,
              save = TRUE,
              order_genes = FALSE)

## select clusters
logFC_Heatmap(dge_data |> filter(cluster %in% my_views), 
              gene_list = row_order,
              title = "Factor3_select", 
              h = pdf_height, 
              w = pdf_width_select, 
              cluster_col = FALSE,
              datatype = "MOFA",
              flip = TRUE,
              save = TRUE,
              order_genes = FALSE)

#### GO on top genes ####
library("org.Hs.eg.db")
library("clusterProfiler")

entrez_search <- bitr(unique(gene_weights$Factor3$feature), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
entrez_search |> count(ENSEMBL) |> count(n)

weight_density <- gene_weights$Factor3 |>
    filter(ctype %in% my_views) |>
    ggplot(aes(x = value, color = ctype)) +
    geom_density()

ggsave(weight_density, filename = here(plot_dir, "weight_density.png"))
    

gene_weights_GO <- gene_weights$Factor3 |>
    left_join(entrez_search, by = c("feature" = "ENSEMBL"), relationship = "many-to-many") |>
    mutate(abs_value = abs(value),
           weight_pos = value > 0,
           datatype = ifelse(grepl("_Sp", ctype), "Visium", "snRNA-seq")) |>
    filter(ctype %in% my_views) |>
    group_by(ctype,weight_pos) |>
    arrange(-abs_value) |>
    mutate(rank = row_number(),
           GO_group = ifelse(rank <= 25, 
                         ifelse(weight_pos, 
                                paste0(gsub("\\.", "", ctype), "_F3+"), 
                                paste0(gsub("\\.", "", ctype), "_F3-")),
                                "None")
    ) |>
    ungroup()

gene_weights_GO |> count(GO_group)

universe <- unique(gene_weights_GO$ENTREZID)
length(universe)

ont_list <- c("CC","BP","MF")
names(ont_list) <- ont_list

## run GO
go_result <- map(ont_list, ~compareCluster(ENTREZID ~ GO_group,
                                           data = gene_weights_GO |> filter(GO_group != "None"), 
                                           OrgDb = org.Hs.eg.db,
                                           fun = enrichGO,
                                           universe = universe,
                                           ont = .x, ##ALL,CC,BP,MF
                                           pAdjustMethod = "BH",
                                           pvalueCutoff = 0.05,
                                           # qvalueCutoff = 0.05,
                                           readable = TRUE))

saveRDS(go_result, file = here(data_dir, sprintf("MOFA_GO_result_%s.rds", opt$datatype)))

## convert to table
compare_clus <- map2_dfr(go_result, names(go_result), ~.x@compareClusterResult |> mutate(ONTOLOGY = .y))

compare_clus |> count(ONTOLOGY, GO_group)
compare_clus |> filter(Count >=5)|> count(ONTOLOGY, GO_group)

## dotplots
pdf(file = here(plot_dir, sprintf("MOFA_GO_dotplot_%s.pdf", opt$datatype)), width = 10, height = 10)
walk2(go_result, names(go_result), 
      ~print(
          dotplot(.x, 
                  x = "GO_group", 
                  showCategory = 5, 
                  label_format = 60)  +
              ggtitle(paste("GO Enrichment:", .y)) +
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
      )
)
dev.off()


#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


