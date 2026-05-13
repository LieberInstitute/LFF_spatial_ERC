## Louise Huuki-Myers, May 2026
## Explore RCTD results

#### Set Up ####
library("SpatialExperiment")
library("qs2")
library("here")
library("sessioninfo")
library("tidyverse")
library("spatialLIBD")
library("DeconvoBuddies")


data_dir <- here("processed-data", "21_Xenium", "10_xenium_cell_types")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "21_Xenium", "10_xenium_cell_types")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE) 
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

#### Load data ####
message(Sys.time(), "- Load xenium data")
spe <- qs_read(here("processed-data", "21_Xenium", "08_xenium_QC_normalize","spe_xenium_QC.qs2"))
spe$sample_id <- spe$BrNum

#### Add RCTD cell types to spe ####
message(Sys.time(), "- Load rctd data")
rctd_data <- qs_read(here("processed-data", "21_Xenium", "09_xenium_label_transfer_RCTD","rctd_results_xenium.qs2"))

names(rctd_data@results)

head(rctd_data@results$results_df)

ncol(spe) == nrow(rctd_data@results$results_df)
identical(colnames(spe), rownames(rctd_data@results$results_df))

spe$cell_id <- colnames(spe)

colData(spe) <- cbind(colData(spe), rctd_data@results$results_df)

#### Explore RCTD results ####

table(is.na(spe$spot_class), spe$sum_gex < 100)

## from before sum_gex < 100 filter
#         FALSE   TRUE
# FALSE 450307      0
# TRUE       0  17209


table(rctd_data@results$results_df$spot_class)
# reject           singlet   doublet_certain doublet_uncertain 
# 5609            304401            122264             18033

table(rctd_data@results$results_df$first_type)

# Astro.1          Astro.2          Astro.3          Astro.4          Astro.5            Macro          Micro.1 
# 15125            18158              945             8954            42636             6285              291 
# Micro.2          Micro.3          Micro.4          Micro.5            OPC.1            OPC.2            OPC.3 
# 14780             2412              345            18178              120            11871               34 
# OPC.4            OPC.5          Oligo.1          Oligo.2          Oligo.3          Oligo.4          Oligo.5 
# 3814            11465            16888            11961            63536            37989             3599 
# Vasc.Endo          Vasc.PC        Vasc.VLMC         Excit.L2     Excit.L2_5.1     Excit.L2_5.2       Excit.L5.1 
# 36177            27006             8918            22124            10358             6601             8213 
# Excit.L5.2    Excit.L5_6_NP      Excit.L6_CT        Excit.L6b Inhib.Chandelier Inhib.Lamp5_Lhx6       Inhib.Pax6 
# 2341              983             5039             4974             1447             5284             1713 
# Inhib.Pvalb        Inhib.Sst        Inhib.Vip 
# 5247             4651             9845 

## spot class summary
rcdt_results_class_summary <- rctd_data@results$results_df |> 
    rownames_to_column("barcode") |>
    mutate(sample_id = gsub("_.*", "", barcode)) |>
    count(sample_id, spot_class) |>
    group_by(sample_id) |>
    mutate(prop = n/sum(n))

## spot class summary plots
rctd_class_bar_plot <- rcdt_results_class_summary |>
    ggplot(aes(x = sample_id, y = n, fill = spot_class)) +
    geom_col() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(rctd_class_bar_plot, filename = here(plot_dir, "xenium_rctd_class_bar_plot.png"))

rctd_class_prop_bar_plot <- rcdt_results_class_summary |>
    ggplot(aes(x = sample_id, y = prop, fill = spot_class)) +
    geom_col() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(rctd_class_prop_bar_plot, filename = here(plot_dir, "xenium_rctd_class_prop_bar_plot.png"))

## cell type summary 
rcdt_results_singlets_summary <- rctd_data@results$results_df |> 
    rownames_to_column("barcode") |>
    filter(spot_class == "singlet") |>
    mutate(sample_id = gsub("_.*", "", barcode)) |>
    group_by(sample_id, spot_class, cell_type_anno = first_type) |>
    summarise(xenium_n = n()) |>
    group_by(sample_id) |>
    mutate(xenium_singlet_prop = xenium_n/sum(xenium_n),
           type = "first_type")

rcdt_results_doublets <- rctd_data@results$results_df |> 
    rownames_to_column("barcode") |>
    filter(spot_class == "doublet_certain") |>
    mutate(sample_id = gsub("_.*", "", barcode)) |>
    select(barcode, spot_class, sample_id, first_type, second_type) |>
    pivot_longer(!c(barcode, sample_id, spot_class), names_to = "type", values_to = "cell_type_anno") |> 
    bind_rows(rctd_data@results$results_df |> 
                  rownames_to_column("barcode") |>
                  filter(spot_class == "doublet_uncertain") |>
                  mutate(sample_id = gsub("_.*", "", barcode),
                         type = "first_type") |>
                  select(barcode, spot_class, sample_id, type,  cell_type_anno = first_type))

rcdt_results_doublets |> count(spot_class, type)


rcdt_results_doublets |> 
    filter(spot_class == "doublet_certain") |>
    count(first_type, second_type)

doublet_combos <- rcdt_results_doublets |> 
    filter(spot_class == "doublet_certain") |>
    group_by(barcode, sample_id) |> 
    summarise(doublet = paste0(sort(cell_type_anno), collapse = "-")) |>
    ungroup() |>
    count(doublet)

doublet_combos |> arrange(-n) |> print(n = 20)

doublet_combos |> filter(grepl("Oligo.3", doublet)) |> arrange(-n)

# doublet                  n
# <chr>                <int>
# 1 Astro.5-Oligo.3       1714
# 2 Oligo.3-Excit.L6_CT   1615
# 3 Oligo.3-Excit.L6b     1608
# 4 Oligo.3-Excit.L2_5.1  1546
# 5 Oligo.3-Excit.L2_5.2  1022

rcdt_results_doublets_summary <- rcdt_results_doublets |>    
    group_by(sample_id, spot_class, type, cell_type_anno) |>
    summarise(xenium_n = n())  |>
    group_by(sample_id) |>
    mutate(xenium_doublet_prop = xenium_n/sum(xenium_n))

rcdt_results_summary <- rcdt_results_singlets_summary |> 
    bind_rows(rcdt_results_doublets_summary)  |>
    group_by(sample_id, cell_type_anno) |>
    summarise(xenium_n = sum(xenium_n)) |>
    group_by(sample_id) |>
    mutate(xenium_prop = xenium_n/sum(xenium_n))


cell_type_class_summary <- rcdt_results_singlets_summary |> 
    bind_rows(rcdt_results_doublets_summary) |>
    group_by(cell_type_anno, spot_class, type) |>
    summarise(xenium_n = sum(xenium_n)) |>
    group_by(cell_type_anno) |> 
    mutate(prop = xenium_n/sum(xenium_n),
           rctd_class = ifelse(spot_class == "doublet_certain", paste0("doublet_c_", gsub("_type", "", type)) , as.character(spot_class)))

rctd_class_ct_bar_plot <- cell_type_class_summary |>
    ggplot(aes(x = cell_type_anno, y = xenium_n, fill = rctd_class)) +
    geom_col() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(rctd_class_ct_bar_plot, filename = here(plot_dir, "xenium_rctd_class_ct_bar_plot.png"), width = 9)

rctd_class_ct_prop_bar_plot <- cell_type_class_summary |>
    ggplot(aes(x = cell_type_anno, y = prop, fill = rctd_class)) +
    geom_col() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(rctd_class_ct_prop_bar_plot, filename = here(plot_dir, "xenium_rctd_class_ct_prop_bar_plot.png"), width = 9)

write_csv(cell_type_class_summary, file = here(data_dir, "RCTD_cell_type_class_summary.csv"))

#### quality metrics vs. RCTD class ####

rctd_qc_summary <- as.data.frame(colData(spe)) |>
    group_by(spot_class) |>
    summarise(across(
        c(sum_gex, detected_gex, cell_area),
        list(
            min    = ~min(.,    na.rm = TRUE),
            Q1     = ~quantile(., 0.25, na.rm = TRUE),
            median = ~median(., na.rm = TRUE),
            mean   = ~mean(.,   na.rm = TRUE),
            Q3     = ~quantile(., 0.75, na.rm = TRUE),
            max    = ~max(.,    na.rm = TRUE)
        ),
        .names = "{.fn}_{.col}"
    )) |>
    pivot_longer(!spot_class, names_to = "stat", values_to = "value") |>
    separate(stat, into = c("summary", "stat"), sep = "_", extra = "merge")


rctd_qc_summary |> filter(stat == "cell_area", spot_class %in% c("singlet", "doublet_certain"))

gg_QC_plot_out <- GGally::ggpairs(as.data.frame(colData(spe)), columns = c("sum_gex", "detected_gex", "cell_area"), aes(colour = spot_class)) + theme_bw()
ggsave(gg_QC_plot_out, filename = here(plot_dir, "xenium_QC_metrics_ggpairs_plot_RCTD_spot_class.png"), height = 12, width = 12)


#### load SingleR results ####

singleR_results <- qs_read(here("processed-data", "21_Xenium", "09_xenium_label_transfer_singleR", "SingleR_results_xenium.qs2"))

rctd_test <- rctd_data@results$results_df[1:1000,]

identical(rownames(singleR_results), rownames(rctd_test))

table(singleR_results$labels, rctd_test$first_type)


#### Save SPE with additional data ####
message(Sys.time(), " - Saving cleaned SPE object")
qs2::qs_save(spe, here(data_dir, "spe_xenium_cell_types.qs2"))

# spe <- qs_read(here("processed-data", "21_Xenium", "10_xenium_cell_types","spe_xenium_cell_types.qs2"))

source(here("code", "utils", "my_plot_reduced_dim.R"))

#### Plot RCTD cell types in Reduced Dims ####
## categorical
walk2(c("first_type", "second_type"), list(cell_type_colors$anno, cell_type_colors$anno), ~my_plot_reduced_dim(spe, prefix = "ERC_xenium", dimred = "TSNE", my_var = .x, var_type = "cat", color_pal = .y))
walk2(c("first_type", "second_type"), list(cell_type_colors$anno, cell_type_colors$anno), ~my_plot_reduced_dim(spe, prefix = "ERC_xenium", dimred = "UMAP", my_var = .x, var_type = "cat", color_pal = .y))

walk2(c("first_type", "second_type"), 
      list(cell_type_colors$anno, cell_type_colors$anno), 
      ~my_plot_reduced_dim(spe, prefix = "ERC_xenium", 
                           dimred = "UMAP", 
                           my_var = .x, 
                           var_type = "cat", 
                           color_pal = .y, 
                           facet = "spot_class"))

walk(c("spot_class"), ~my_plot_reduced_dim(spe, prefix = "ERC_xenium", dimred = "TSNE", my_var = .x, var_type = "cat"))
walk(c("spot_class"), ~my_plot_reduced_dim(spe, prefix = "ERC_xenium", dimred = "TSNE", my_var = .x, var_type = "cat", facet = TRUE))

walk(c("spot_class"), ~my_plot_reduced_dim(spe, prefix = "ERC_xenium", dimred = "UMAP", my_var = .x, var_type = "cat"))
walk(c("spot_class"), ~my_plot_reduced_dim(spe, prefix = "ERC_xenium", dimred = "UMAP", my_var = .x, var_type = "cat", facet = TRUE))


#### Explore doublet pairings ####

rctd_tb_doublet_counts <- rctd_data@results$results_df |>
    filter(spot_class == "doublet_certain") |>
    count(first_type, second_type) |>
    as_tibble()

rctd_tb_doublet_counts |> arrange(-n) |> print(n=20)

rctd_tb_doublet_counts |>
    mutate(astro5_involved = str_detect(first_type, "Astro.5") | str_detect(second_type, "Astro.5")) |>
    group_by(astro5_involved) |>
    summarize(total_doub = sum(n)) |>
    ungroup() |>
    mutate(pct = total_doub/sum(total_doub))

# astro5_involved total_doub   pct
# <lgl>                <int> <dbl>
# 1 FALSE                89456 0.732
# 2 TRUE                 32808 0.268

rctd_tb_doublet_counts |> filter(first_type == "Oligo.3") |> arrange(-n) 
    

rctd_tb_doublet_tile <- rctd_tb_doublet_counts |>
    ggplot(aes(first_type, second_type, fill = n)) +
    geom_tile() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave(rctd_tb_doublet_tile, filename = here(plot_dir, "xenium_rctd_doublet_tile.png"))

#### Spatial plots ####

spe$Oligo <- droplevels(spe$first_type)
spe$Oligo[!grepl("Oligo", spe$Oligo)] <- NA

# load(here("processed-data","00_project_prep","Oligo_OPC_colors.Rdata"), verbose = TRUE)

Oligo_colors <- c(Oligo.1  = "#57DBCE", 
                  Oligo.4  = "#68DE6E",
                  # Oligo.4  = "#20978C",
                  Oligo.2  = "#0072E5",  
                  Oligo.3 = "#FF700A", 
                  Oligo.5  = "#FF5EFA")

walk(unique(spe$BrNum), function(samp){
    
    # vis_clus_class <- vis_clus(spe, 
    #                            sampleid = samp,
    #                            clustervar = "spot_class",
    #                            datatype = "Xenium", 
    #                            point_size = 1.5,
    #                            # alpha = 0.5,
    #                            colors = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF"),
    #                            guide_point_size = 2)
    # 
    # ggsave(vis_clus_class, filename = here(plot_dir, sprintf("Xenium_vis_spot_class_%s.png", samp)), width = 12)
    # 
    
    # vis_clus_first_type <- vis_clus(spe, 
    #                                 sampleid = "Br1556",
    #                                 clustervar = "first_type",
    #                                 datatype = "Xenium", 
    #                                 point_size = 1.5,
    #                                 colors = cell_type_colors$anno,
    #                                 guide_point_size = 2)
    # 
    # ggsave(vis_clus_first_type, filename = here(plot_dir, "Xenium_vis_first_type_Br1556.png"), width = 12)
    # 
    
    # vis_clus_first_type_singlet <- vis_clus(spe[, spe$spot_class == "singlet"], 
    #                                         sampleid = samp,
    #                                         clustervar = "first_type",
    #                                         datatype = "Xenium", 
    #                                         point_size = 1.5,
    #                                         colors = cell_type_colors$anno,
    #                                         guide_point_size = 2)
    # 
    # ggsave(vis_clus_first_type_singlet, filename = here(plot_dir, sprintf("Xenium_vis_first_type_singlet_%s.png", samp)), width = 12)
    # 
    
    # vis_clus_first_type_Oligo.3 <- vis_clus(spe[, spe$spot_class == "singlet" & spe$first_type == "Oligo.3"], 
    #                                         sampleid = "Br1556",
    #                                         clustervar = "first_type",
    #                                         datatype = "Xenium", 
    #                                         point_size = 1.5,
    #                                         colors = cell_type_colors$anno,
    #                                         guide_point_size = 2)
    # 
    # ggsave(vis_clus_first_type_Oligo.3, filename = here(plot_dir, "Xenium_vis_first_type_Oligo.3_Br1556.png"), width = 12)
    # 
    vis_clus_first_type_Oligo <- vis_clus(spe[, spe$spot_class == "singlet"], 
                                            sampleid = samp,
                                            clustervar = "Oligo",
                                            datatype = "Xenium", 
                                            point_size = 1.5,
                                            colors = Oligo_colors,
                                            guide_point_size = 2)
    
    ggsave(vis_clus_first_type_Oligo, filename = here(plot_dir, sprintf("Xenium_vis_first_type_Oligo_%s.png", samp)), width = 12)
    
    
})



#### Check marker gene expression ####

probes_long <- readxl::read_xlsx(here("processed-data", "21_Xenium", "02_xenium_compile_custom", "ERC_Xenium_ALL_probes_long.xlsx"))

probes_markers <- probes_long |>
    filter(grepl("Marker|Base", goal), !grepl("~|_Sp|Proliferation|Glioblastoma|Broad", target) ) |>
    mutate(target_broad = gsub('\\..*$',"",target)) |>
    group_by(target_broad, gene_name) |>
    summarise(n_target = n(),
              target = list(sort(unique(target))),
              targets = paste0(unlist(target), collapse = ", "),
              goal = list(sort(unique(goal)))
    ) |>
    group_by(gene_name) |>
    mutate(n_broad = n(),
           unique_target = n_broad ==1 & n_target == 1) |>
    ungroup()

probes_markers |> count(target_broad)
probes_markers |> count(n_broad)
probes_markers |> count(unique_target)

probes_markers |> filter(n_broad == 1) |> count(target_broad, unique_target, n_target)

probes_markers |> filter(gene_name == 'DNER')

probes_markers |> ungroup() |> count(gene_name) |> arrange(-n)

my_genes <- c(Astro = "AQP4", 
              # Micro = "CD86",
              Micro = "CTSH",
              Oligo = "MBP",
              Oligo.M = "OPALIN",
              Oligo.3 = "LINGO2",
              OPC = "PDGFRA",
              Vasc = "PECAM1",
              Excit = "SLC17A7",
              Inhib = "GAD1")

broad_markers <- list(
    Astro = c("AQP4", "SOX9", "GJA1", "FGFR3", "TGFB2", "SPON1"),
    Oligo = c("MBP", "MOBP", "MOG", "MAG", "OPALIN","LINGO2"),
    OPC   = c("PDGFRA", "OLIG2", "VCAN", "SOX10", "PTPRZ1", "BRINP3"),
    Micro = c("P2RY12", "SPI1", "CTSH", "TREM2", "P2RY13", "ITGAM"),
    Macro = c("CD163", "LYVE1", "PTPRC", "ITGAM", "ITGAX", "TGFB1"),
    Vasc  = c("PECAM1", "FLT1", "NRP1", "ABCC9", "CSPG4", "NR2F2"),
    Excit = c("SLC17A7", "SLC17A6", "CUX2", "RORB", "THEMIS", "CRYM"),
    Inhib = c("GAD1", "GAD2", "SST", "PVALB", "VIP", "LAMP5")
)

# map(broad_markers, ~all(.x %in% rownames(spe)))

plot_marker_express_List(
    spe[, spe$spot_class == "singlet"],
    gene_list = broad_markers,
    pdf_fn = here(plot_dir, "xenium_rcdt_singlets_broad_marker_expres_violin.pdf"),
    cellType_col = "first_type",
    gene_name_col = "Symbol",
    color_pal = cell_type_colors$anno
)

plot_marker_express_List(
    spe[, spe$spot_class == "doublet_certain"],
    gene_list = broad_markers,
    pdf_fn = here(plot_dir, "xenium_rcdt_doubletc_broad_marker_expres_violin.pdf"),
    cellType_col = "first_type",
    gene_name_col = "Symbol",
    color_pal = cell_type_colors$anno
)

#### scDot plots ####

library("scDotPlot")

# slurmjobs::job_single('10_xenium_cell_types', create_shell = TRUE, memory = '100G', command = "Rscript 10_xenium_cell_types.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

