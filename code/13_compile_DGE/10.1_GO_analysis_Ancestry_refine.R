## Louise Huuki-Myers, June 2025
## Compile and plot all DGE data

#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")

plot_dir <- here("plots", "13_compile_DGE", "10.1_GO_analysis_Ancestry_refine")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### load data ####

# list.files(here("processed-data", "13_compile_DGE", "10_GO_analysis_contrast", "GO_ancestry"))
compare_clus <- readRDS(here("processed-data", "13_compile_DGE", "10_GO_analysis_contrast", "GO_ancestry", "GO_compare_clus_ancestry_sn_fine.rds"))

compare_clus |> filter(grepl("Oligo.3", Cluster)) |> count(Cluster)

compare_clus |> select(Cluster, Description, geneID) |> as_tibble()

## validate data

valid_DEGs <- read_csv(here("processed-data", "13_compile_DGE", "19_validate_summary", "DGE_Xenium_validated_All.csv"))

valid_DEGs_Oligo3 <-  valid_DEGs |> filter(cell_type_anno == "Oligo.3") |> pull(gene_name) |> unique()


compare_clus |> count(ONTOLOGY)

compare_clus_valid <- compare_clus |> 
    select(Cluster, Description, geneID) |> 
    as_tibble() |>
    mutate(gene_list = str_split(geneID, "/"),
           n_valid = map_int(gene_list, ~sum(.x %in% valid_DEGs_Oligo3)))

compare_clus_valid |> count(n_valid)
compare_clus_valid |> filter()


