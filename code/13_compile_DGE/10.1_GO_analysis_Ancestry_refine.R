## Louise Huuki-Myers, June 2025
## Compile and plot all DGE data

#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")

plot_dir <- here("plots", "13_compile_DGE", "10.1_GO_analysis_Ancestry_refine")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

data_dir <- here("processed-data", "13_compile_DGE", "10.1_GO_analysis_Ancestry_refine")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### load data ####

# list.files(here("processed-data", "13_compile_DGE", "10_GO_analysis_contrast", "GO_ancestry"))
compare_clus <- readRDS(here("processed-data", "13_compile_DGE", "10_GO_analysis_contrast", "GO_ancestry", "GO_compare_clus_ancestry_sn_fine.rds"))

compare_clus |> filter(grepl("Oligo.3", Cluster)) |> count(Cluster)

compare_clus |> select(Cluster, Description, geneID) |> as_tibble()

## validate data

valid_DEGs <- read_csv(here("processed-data", "13_compile_DGE", "19_validate_summary", "DGE_Xenium_validated_All.csv"))

valid_DEGs_Oligo3 <-  valid_DEGs |> filter(cell_type_anno == "Oligo.3") |> pull(gene_name) |> unique()

DGE_results_ancestry <- readRDS(here("processed-data", "13_compile_DGE", "05_compile_DGE_ancestry", "Xenium", "DGE_results_ancestry_Xenium_wSN.Rds"))

valid_genes_anc <- DGE_results_ancestry |> filter(validate) |> pull(gene_name)

DGE_results_ancestry |> filter(validate) |>select(contrast, gene_name, vlmf_sn_t) |> mutate(valid_carrier = gene_name %in% valid_DEGs_Oligo3)
# contrast   gene_name vlmf_sn_t valid_carrier
# <chr>      <chr>         <dbl> <lgl>        
# 1 carrier_AA SLC24A2       -3.32 TRUE         
# 2 carrier_AA SNX19         -4.05 FALSE        
# 3 carrier_EA RABL2B        -3.66 FALSE        
# 4 carrier_EA RFTN1          4.23 FALSE        
# 5 carrier_EA NTRK3          3.81 TRUE  

intersect(valid_genes_anc, valid_DEGs_Oligo3)
# [1] "SLC24A2" "NTRK3"

compare_clus |> count(ONTOLOGY)

compare_clus_valid <- compare_clus |> 
    filter(grepl("Oligo.3", Cluster)) |>
    select(Cluster, Description, geneID, p.adjust, GeneRatio, ONTOLOGY) |> 
    as_tibble() |>
    mutate(matched = map(str_split(geneID, "/"), ~intersect(.x, valid_DEGs_Oligo3)),
           n_valid = map_int(matched, length),
           valid_genes = map_chr(matched, paste, collapse = "/"),
           matched_anc = map(str_split(geneID, "/"), ~intersect(.x, valid_genes_anc)),
           n_valid_anc = map_int(matched_anc, length),
           valid_genes_anc = map_chr(matched_anc, paste, collapse = "/"),
           GeneRatio_num = as.numeric(sub("/.*", "", GeneRatio)) / as.numeric(sub(".*/", "", GeneRatio))
           )

write_csv(compare_clus_valid |> select(-matched, -matched_anc),
          file = here(data_dir, "GO_compare_clus_ancestry_sn_fine_refine.csv"))


compare_clus_valid |> count(n_valid)
compare_clus_valid |> filter(n_valid_anc > 0) |> count(valid_genes_anc)

# valid_genes_anc     n
# <chr>           <int>
# 1 NTRK3              13
# 2 RFTN1               4

compare_clus_valid |> filter(n_valid > 1 | valid_genes_anc > 0)
compare_clus_valid |> filter(valid_genes_anc > 0)


plot_terms <- c(
    # ## Down
    # "myelination", "ensheathment of neurons", "oligodendrocyte differentiation", "myelin sheath", 
    # AA / synaptic-vesicle-calcium axis
    "calcium ion transport", "neurotransmitter transport", "synaptic vesicle cycle",
    "modulation of chemical synaptic transmission", "regulation of trans-synaptic signaling",
    "GABA-ergic synapse",
    "regulation of neuron projection development",
    # EA / ECM-growth factor-caveola axis
    "caveola", "extracellular matrix organization", "muscle tissue development",
    "growth factor binding"
    # # membrane microdomain thread (RFTN1-linked, ancestry-validated, ALL-only)
    # "membrane raft", "receptor complex",
)

cluster_order <- c("Oligo-3_down_ALL", "Oligo-3_up_ALL", "Oligo-3_up_AA", "Oligo-3_up_EA")

compare_clus_select <- compare_clus_valid |>
    filter(Description %in% plot_terms) |>
    mutate(Description = factor(str_wrap(Description, width = 35), levels = rev(str_wrap(plot_terms, width = 35))),
           Cluster = factor(Cluster, levels = cluster_order),
           validated = case_when(n_valid_anc > 0 ~ "ancestry-validated",
                                 n_valid > 0 ~ "carrier-validated",
                                 TRUE ~ "not validated"))

compare_clus_select |> 
    filter(validated == "ancestry-validated") |> 
    select(Cluster, Description, geneID,valid_genes, valid_genes_anc) |>
    arrange(valid_genes_anc)

# Cluster        Description                                    geneID                     valid_genes valid_genes_anc
# <fct>          <fct>                                          <chr>                      <chr>       <chr>          
# 1 Oligo-3_up_ALL "receptor complex"                             CALCRL/RAMP1/MERTK/CSF2RA… NTRK3       NTRK3          
# 2 Oligo-3_up_ALL "regulation of neuron projection\ndevelopment" UST/BAIAP2/LRRK2/CNR1/CRM… LRRK2/NTRK… NTRK3          
# 3 Oligo-3_up_ALL "growth factor binding"                        COL4A1/COL6A1/A2M/IGFBP5/… IGFBP5/NTR… NTRK3          
# 4 Oligo-3_up_EA  "growth factor binding"                        COL1A2/COL4A1/IGFBP7/NTRK… NTRK3       NTRK3          
# 5 Oligo-3_up_ALL "membrane raft"                                SORBS1/ADCY8/LRRK2/CNR1/S… LRRK2       RFTN1 

## Custom GO plot 

oligo3_go_select_plot <- compare_clus_select |>
    ggplot(aes(x = Cluster, y = Description)) +
    geom_point(aes(size = GeneRatio_num, fill = p.adjust, color = validated), shape = 21, stroke = 1.2) +
    facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradient(low = "red", high = "blue", trans = "log10") +
    scale_color_manual(values = c("ancestry-validated" = "black", "carrier-validated" = "grey40", "not validated" = "grey80")) +
    theme_bw() +
    labs(title = "Oligo.3 ancestry GO terms",
         x = NULL, y = NULL, size = "Gene Ratio", fill = "p.adjust", color = "Xenium validation") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 

ggsave(oligo3_go_select_plot, filename = here(plot_dir, "Oligo3_ancestry_GO_select_dotplot.png"), height = 5, width = 5)
