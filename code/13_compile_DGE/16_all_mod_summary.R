## Louise Huuki-Myers, April 2026
## Create summary plots all DE models

#### Set up ####
library("tidyverse")
library("here")
library("SummarizedExperiment")

## set plot dir
data_dir <- here("processed-data", "13_compile_DGE", "16_all_mod_summary")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

## set plot dir
plot_dir <- here("plots", "13_compile_DGE", "16_all_mod_summary")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)


#### Get DGE summary of n_gene + nFDR05 ####

datatypes <- c("Visium", "sn_broad", "sn_fine")
names(datatypes) <- datatypes

mod_dir <- c(carrier = "01_compile_DGE",
             carrier_anc = "05_compile_DGE_ancestry")

# list.files(here("processed-data", "12_voomLmFit", "01_Clusterwise_voomLmFit", paste0("vlmf_","sn_fine")), pattern = "FDR05_summary")

summary_files <- map2(mod_dir, names(mod_dir), function(mod, mod_name){
    unlist(map(c("Visium", "sn_broad", "sn_fine"), ~list.files(here("processed-data", "13_compile_DGE", mod, .x),
                                            pattern = "summary",
                                            full.names = TRUE))
           ) })

# head(read.csv(summary_files$carrier[[1]]))

summary_tbl <- map2_dfr(summary_files, names(summary_files), 
                                    function(files, name){map(files, ~read.csv(.x, row.names = 1) |> 
                                                                  mutate(model = name,
                                                                         file = basename(.x),
                                                                         datatype = gsub(".*_summary_(.*)\\.csv", "\\1", file))
                                    )}
) |>
    filter(mod == "carrier" | is.na(mod)) |>
    dplyr::rename(n_gene = n_genes)

head(summary_tbl)

summary_tbl |> dplyr::count(model, datatype)

summary_interaction_files <- map(datatypes, ~list.files(here("processed-data", "13_compile_DGE", "05_compile_DGE_interaction", .x),
                                              pattern = "summary",
                                              full.names = TRUE))

summary_interaction_tbl <- map2_dfr(summary_interaction_files, names(summary_interaction_files), 
                                function(files, name){map(files, ~read.csv(.x, row.names = 1) |> 
                                                              mutate(datatype = name, 
                                                                     file = basename(.x),
                                                                     model = gsub(".*_summary_(interaction_[^_]+)_.*", "\\1", file)
                                                              )
                                                          )}
                                ) |>
    select(-n_FDR10, -n_FDR20)

summary_tbl_all <- summary_tbl  |>
    bind_rows(summary_interaction_tbl)

any(is.na(summary_tbl_all$n_gene))
    
summary_tbl_all |> dplyr::count(model, datatype)

#### Get n_donors ####

n_donor_all <- map_dfr(datatypes, function(datatype){
    
    if(datatype == "Visium"){
        pb_fn <- here("processed-data", "09_pseudoBulkDGE_Visium", "01_pseudobulk_data_Visium", "spe_pseudo_DGE.RDS")
    } else if(datatype == "sn_broad"){
        pb_fn <- here("processed-data", "08_pseudoBulkDGE_sn", "01_pseudobulk_data_sn","sce_pseudo_DGE-cell_type_broad.RDS")
    }else if(datatype == "sn_fine"){
        pb_fn <- here("processed-data", "08_pseudoBulkDGE_sn", "01_pseudobulk_data_sn","sce_pseudo_DGE-cell_type_anno.RDS")
    } else {
        stop("non-valid datatype")
    }
    
    message(Sys.time(), sprintf(" - Datatype = %s, loading '%s'", datatype, basename(pb_fn)))
    sce_pb <- readRDS(pb_fn)
    sce_pb <- sce_pb[,sce_pb$BrNum != 'Br1289'] ## Drop sample Br1289

    pb <- as.data.frame(colData(sce_pb))
    
    # carrier & interaction
    n_donor <- pb |> 
        group_by(cluster = registration_variable) |>
        summarise(n_donor = n()) |>
        mutate(datatype = datatype)
    
    # carrier w/ ancestry contrast
    n_donor_contrast <- pb |> 
        group_by(cluster = registration_variable,
                 contrast = Ancestry) |>
        summarise(n_donor = n()) |>
        mutate(datatype = datatype,
               contrast = paste0("carrier_", contrast))
    
    
    n_donor <- bind_rows(n_donor, n_donor_contrast) |>
        mutate(cluster = gsub("_Sp", "~Sp", cluster))
    return(n_donor)
        
})

n_donor_all |> filter(datatype == "Visium")

n_donor_all |> group_by(datatype, contrast) |> summarise(min_donor = min(n_donor), max_donor = max(n_donor))

# datatype contrast   min_donor max_donor
# <chr>    <chr>          <int>     <int>
# 1 Visium   carrier_AA        16        16
# 2 Visium   carrier_EA        14        14
# 3 Visium   NA                30        30
# 4 sn_broad carrier_AA         8        16
# 5 sn_broad carrier_EA        10        14
# 6 sn_broad NA                18        30
# 7 sn_fine  carrier_AA         4        16
# 8 sn_fine  carrier_EA         5        14
# 9 sn_fine  NA                 9        30

summary_tbl_all <- summary_tbl_all |> left_join(n_donor_all)

write_csv(summary_tbl_all, here(data_dir, "DGE_summary_all_model.csv"))

n_donor_all <- n_donor_all |>
    mutate(group = fct_relevel(ifelse(is.na(contrast), "All", gsub("carrier_","", contrast)),"All"),
           datatype = factor(datatype, levels = datatypes))

write_csv(n_donor_all, here(data_dir, "DGE_n_donor.csv"))

#### plots ####
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

## boxplot of n donors
donor_boxplot <- n_donor_all |>
    ggplot(aes(y = datatype, x = n_donor, color = group)) +
    geom_boxplot() +
    theme_bw() +
    scale_color_manual(values = c(ancestry_colors, All = "black")) +
    scale_y_discrete(limits = rev) +
    labs(x = "number of donors") +
    theme(legend.position = "left")

ggsave(donor_boxplot, filename = here(plot_dir, "n_donor_boxplot.png"), height = 5, width = 4)

sum_FDR05_tbl <- summary_tbl_all |> 
    mutate(model = ifelse(model == "carrier_anc", 
                  paste0(model, gsub("carrier_", ":", contrast)),
                  model),
           datatype = factor(datatype, levels = datatypes)) |>
    group_by(model, datatype) |>
    summarize(sum_FDR05 = sum(n_FDR05))


model_FDR05_tileplot <- sum_FDR05_tbl |>
    ggplot(aes(y = datatype, x = model)) +
    geom_tile(aes(fill = sum_FDR05), color = "grey50") +
    geom_text(aes(label = sum_FDR05, color = sum_FDR05 > 500)) + 
    scale_y_discrete(limits = rev) +
    scale_x_discrete(limits = rev) +
    scale_color_manual(values = c("TRUE" = "white", "FALSE" = "black"), guide = "none") +
    scale_fill_distiller(palette = "Blues", direction = 1) +
    theme_bw() +
    theme(legend.position = "top", 
          legend.text = element_text(angle = 45, hjust = 1),
          axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(model_FDR05_tileplot, filename = here(plot_dir, "model_FDR05_tileplot.png"), height = 5, width = 6)

summary_tbl_all |> 
    ggplot(aes(x = n_gene, y = n_FDR05, color = model)) +
    geom_point()

gene_v_donor_FDRpoint <- summary_tbl_all |> 
    ggplot(aes(x = n_gene, y = n_donor, size = n_FDR05, color = n_FDR05)) +
    geom_point() +
    facet_grid(model~datatype) +
    theme_bw()

ggsave(gene_v_donor_FDRpoint, filename = here(plot_dir, "gene_v_donor_FDRpoint.png"))

model_dotplot <- summary_tbl_all |> 
    ggplot(aes(x = cluster, y = n_donor, size = n_FDR05, color = n_gene)) +
    geom_point() +
    facet_grid(model~datatype, scales = "free_x") +
    theme_bw()

ggsave(model_dotplot, filename = here(plot_dir, "model_dotplot.png"))






