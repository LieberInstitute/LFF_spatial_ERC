## Louise Huuki-Myers, June 2025
## Compile and plot all DGE data

#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")
library("getopt")
library("spatialLIBD")
library("jaffelab")

# Import command-line parameters
scec <- matrix(
    c("datatype", "d", "1", "character", "Data type"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)

## test
# opt$datatype = "sn_broad"
# opt$datatype = "sn_fine"
# opt$datatype = "Visium"

data_dir <- here("processed-data", "13_compile_DGE", "12_DE_enrichment")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "13_compile_DGE", "12_DE_enrichment")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

load(here("processed-data", "project_colors.Rdata"))

if(opt$datatype == "sn_broad"){
    load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE)
    cluster_colors <- cell_type_colors$broad
    cluster_levels <- names(cell_type_colors$broad)
    
    modeling_fn <- here("processed-data", "04_snRNA-seq", "29_sn_subcluster_model_pseudobulk", "sce_subcluster_modeling_results-cell_type_broad.rds")
    
}else if(opt$datatype == "sn_fine"){
    load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE)
    cluster_colors <- cell_type_colors$anno
    cluster_levels <- names(cell_type_colors$anno)
    
    modeling_fn <- here("processed-data", "04_snRNA-seq", "29_sn_subcluster_model_pseudobulk", "sce_subcluster_modeling_results-cell_type_anno.rds")
    
}else if(opt$datatype == "Visium"){
    load(here("processed-data", "SpD_colors.Rdata"), verbose = TRUE)
    cluster_colors <- SpD_colors
    cluster_levels <- names(SpD_colors)
    
    modeling_fn <- here("processed-data", "05_spe_correct_cluster", "20_model_pseudobulk_anno", "modeling_results-SpD.rds")
}

cluster_levels <- cluster_levels[cluster_levels != "Other"]

## load DE data ##
DE_data_fn <- here("processed-data", "13_compile_DGE", "01_compile_DGE", opt$datatype, sprintf("DGE_results_carrier_%s.Rds", opt$datatype))
file.exists(DE_data_fn)

DE_data <- readRDS(DE_data_fn)

cluster_levels <- cluster_levels[cluster_levels %in% DE_data$cluster]

DE_data <- DE_data |> mutate(cluster = factor(cluster, levels = cluster_levels))
DE_data |> filter(vlmf_adj.P.Val < 0.05) |> dplyr::count(cluster)


## load modeling data
modeling_data <- readRDS(modeling_fn)
modeling_data$enrichment$ensembl <- modeling_data$enrichment$gene


#### Gene Sets ####

## load Mathys DE data ##
suppressMessages(source(here("external-data", "Mathys2019", "get_Mathys_DE_data.R")))

## just no vs. path data

Mathys_DE_data_path <- Mathys_DE_data |> 
    filter(DE_type == "no_v_path") |>
    mutate(reg = ifelse(IndModel.FC > 0, "up", "down"),
           DE_ct = ifelse(DEGs.Ind.Mix.models, cell_type, NA),
           DE_ct_reg = ifelse(DEGs.Ind.Mix.models, paste0(DE_ct, "_", reg), NA))

Mathys_DE_data_path |> dplyr::count(cell_type)

## match counts in Fig 1b
# Mathys_DE_data_path |> filter(DEGs.Ind.Mix.models) |> dplyr::count(cell_type)
# Mathys_DE_data_path |> filter(DEGs.Ind.Mix.models) |> dplyr::count(cell_type, IndModel.FC > 0)
Mathys_DE_data_path |> dplyr::count(DE_ct_reg)

# DE_ct_reg     n
# <chr>     <int>
# 1 Ast_down     32
# 2 Ast_up       37
# 3 Ex_down     565
# 4 Ex_up       191
# 5 In_down      51
# 6 In_up         3
# 7 Mic_down     13
# 8 Mic_up       22
# 9 Oli_down     71
# 10 Oli_up      102
# 11 Opc_down     18
# 12 Opc_up       10

mathys_ct_DEGs <- map(splitit(Mathys_DE_data_path$DE_ct), ~Mathys_DE_data_path$gene_name[.x])
mathys_ct_reg_DEGs <- map(splitit(Mathys_DE_data_path$DE_ct_reg), ~Mathys_DE_data_path$gene_name[.x])

## mathys no vs. early data

Mathys_DE_data_early <- Mathys_DE_data |> 
    filter(DE_type == "no_v_early") |>
    mutate(reg = ifelse(IndModel.FC > 0, "up", "down"),
           DE_ct = ifelse(DEGs.Ind.Mix.models, cell_type, NA),
           DE_ct_reg = ifelse(DEGs.Ind.Mix.models, paste0(DE_ct, "_", reg), NA))

Mathys_DE_data_early |> dplyr::count(DE_ct_reg)

mathys_early_ct_DEGs <- map(splitit(Mathys_DE_data_early$DE_ct), ~Mathys_DE_data_early$gene_name[.x])
mathys_early_ct_reg_DEGs <- map(splitit(Mathys_DE_data_early$DE_ct_reg), ~Mathys_DE_data_early$gene_name[.x])


gene_set_list <- list(Mathys = mathys_ct_DEGs,
                      Mathys_reg = mathys_ct_reg_DEGs,
                      Mathys_early = mathys_ct_DEGs,
                      Mathys_early_reg = mathys_ct_reg_DEGs
                      )

map_int(gene_set_list[[1]], length)

## source function
source(here("code", "13_compile_DGE", "DE_gene_set_enrichment.R"))

map2(gene_set_list, names(gene_set_list), function(gene_set, name){
    
    message("gene set enrichment for: ", name)
    
    map_int(gene_set, length) 

    #### cluster modeling ####
    
    cluster_gse <- gene_set_enrichment(
        gene_list = gene_set,
        modeling_results = modeling_data,
        model_type = "enrichment"
    ) |>
        mutate(gene_set = name)
    
    write_csv(cluster_gse, file = here(data_dir, sprintf("Cluster_gene_enrichment_%s_%s.csv", opt$datatype, name)))
    cluster_gse |> filter(Pval < 0.1)
    
    ## plot enrichment
    pdf(here(plot_dir, sprintf("Cluster_gene_enrichment_%s_%s.pdf", opt$datatype, name)))
    print(gene_set_enrichment_plot(cluster_gse))
    dev.off()
    
    #### DE enrichment ####
    
    ## combine reg
    DE_gse <- DE_gene_set_enrichment(gene_set, DE_data) |>
        mutate(gene_set = name)
    
    DE_gse |> filter(Pval < 0.1)
    DE_gse |> filter(NumSig > 0) |> arrange(Pval)
    
    write_csv(DE_gse, file = here(data_dir, sprintf("DE_gene_enrichment_%s_%s.csv", opt$datatype, name)))
    
    ## enrichment heatmap
    pdf_height = 5 + length(cluster_levels)/4
    
    pdf(here(plot_dir, sprintf("DE_gene_enrichment_%s_%s.pdf", opt$datatype, name)), height = pdf_height)
    print(gene_set_enrichment_plot(DE_gse))
    dev.off()
    
    return(map_int(gene_set, length) )
    
})


# slurmjobs::job_loop(loops = list(datatype = c("sn_broad","sn_fine","Visium")),
#                     create_shell = TRUE,
#                     name = "12_DE_enrichment",
#                     create_script = FALSE)


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


