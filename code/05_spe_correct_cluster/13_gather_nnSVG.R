## October 2024, Louise Huuki-Myers
## Load sample-wise nnSVG data, average the ranks of the SVGs
## Examine gene set similarities & differences 

library("SpatialExperiment")
library("HDF5Array")
library("spatialLIBD")
library("tidyverse")
library("UpSetR")
library("here")
library("sessioninfo")

#### load spe data ####
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC"))

rowData(spe)

#### define dirs ####
data_dir <- here("processed-data", "05_spe_correct_cluster", "13_gather_nnSVG")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "05_spe_correct_cluster", "13_gather_nnSVG")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### load nnSVG outputs ####

nnSVG_files <- list.files(here("processed-data", "05_spe_correct_cluster", "12_nnSVG"), full.names = TRUE)
names(nnSVG_files) <- gsub("nnSVG_(Br[0-9]+).csv", "\\1", basename(nnSVG_files))

nnSVG_data <- do.call("rbind", map2(nnSVG_files, names(nnSVG_files), ~read.csv(.x, row.names = 1) |> mutate(sample_id = .y)))

head(nnSVG_data)

nnSVG_data |> count(sample_id) |> arrange(n)
#    sample_id    n
# 1     Br5276  304
# 2     Br5529  848
# 3     Br6263 1159
# 4     Br6321 1555
# 5     Br6161 1602

dlpfc_layers_top100 <-  read_csv(here("processed-data", "00_project_prep", "06_marker_genes", "dlpfc_layers_top100.csv"))
length(unique(dlpfc_layers_top100$ensembl)) # [1] 1129

## other gene sets
load(here("processed-data", "02_build_spe","01_preprocess_spe", "top_hvgs.Rdata"), verbose = TRUE)
load(here("processed-data", "05_spe_correct_cluster", "03_GLM_Harmony", "top_hdgs.rdata"), verbose = TRUE)

# hdg_vs_hvg_overlap <- intersect(top.hvgs$p1, top.hdgs$`2k`)
# length(hdg_vs_hvg_overlap) #[1] 405
# 
# SVG_vs_hvg_overlap <- intersect(top.hvgs$p1, top.svg)
# SVG_vs_hdg_overlap <- intersect(top.hdgs$`1k`, top.svg)
# 
# length(SVG_vs_hvg_overlap) #[1] 176
# length(SVG_vs_hdg_overlap) #[1] 464

nnSVG_avg <- nnSVG_data |>
    group_by(gene_id, gene_name) |>
    summarize(nnsvg_avg_rank = mean(rank),
              n = n()) |>
    arrange(nnsvg_avg_rank) |>
    ungroup() |>
    mutate(top.svg = (n > 5 & nnsvg_avg_rank < 1500),
           DLPFC_layer_marker = gene_id %in% dlpfc_layers_top100$ensembl,
           top.hvg = gene_id %in% top.hvgs$p1,
           top.hdg = gene_id %in% top.hdgs$`1k`) 

nnSVG_avg |> count(top.svg)

nnSVG_avg |> count(top.svg, DLPFC_layer_marker)

nnSVG_avg |> count(top.svg, DLPFC_layer_marker)
nnSVG_avg |> count(top.svg, top.hvg, top.hdg, DLPFC_layer_marker)

write_csv(nnSVG_avg, file = here(data_dir, "nnSVG_avg.csv"))

# gene_id         gene_name nnsvg_avg_rank     n
# <chr>           <chr>              <dbl> <int>
# 1 ENSG00000197971 MBP                 2       31
# 2 ENSG00000131095 GFAP                9.29    31
# 3 ENSG00000106809 OGN                10        1
# 4 ENSG00000123560 PLP1               10.5     31
# 5 ENSG00000118271 TTR                12.7     18
# 6 ENSG00000132639 SNAP25             19.6     31
# 7 ENSG00000154146 NRGN               26.5     31
# 8 ENSG00000142173 COL6A2             29        1
# 9 ENSG00000165507 DEPP1              29.5      2
# 10 ENSG00000171617 ENC1               29.6     31

nnSVG_n_vs_avg_rank <- nnSVG_avg |>
    ggplot(aes(x = n, y = nnsvg_avg_rank, color = top.svg)) +
    geom_point()

ggsave(nnSVG_n_vs_avg_rank, filename = here(plot_dir, "nnSVG_n_vs_avg_rank.png"))

nnSVG_n_vs_avg_rank_DLPFC <- nnSVG_avg |>
    ggplot(aes(x = n, y = nnsvg_avg_rank, color = DLPFC_layer_marker)) +
    geom_point()

ggsave(nnSVG_n_vs_avg_rank_DLPFC, filename = here(plot_dir, "nnSVG_n_vs_avg_rank-DLPFC_markers.png"))

nnSVG_n_histogram <- nnSVG_avg |>
    ggplot(aes(x = n)) +
    geom_histogram(binwidth = 1)

ggsave(nnSVG_n_histogram, filename = here(plot_dir, "nnSVG_n_histogram.png"))

## genes were rarely evaluated in all samples
nnSVG_avg |> ungroup() |> count(n == 31)
# `n == 31`     n
# <lgl>     <int>
# 1 FALSE      5614
# 2 TRUE        265

nnSVG_avg |> ungroup() |> count(n >= 20)


top.svg <- nnSVG_avg |> 
    ungroup() |>
    filter(n > 5,
           nnsvg_avg_rank < 1500) |>
    pull(gene_id)

length(top.svg)
# [1] 1454

save(top.svg, file = here(data_dir, "top_svg.Rdata"))

#### upset plots ####

gene_list = list(HVG_10p = top.hvgs$p1,
                 HDG_1k = top.hdgs$`1k`,
                 SVG = top.svg,
                 spatialDLPFC = unique(dlpfc_layers_top100$ensembl))

## Selected gene list only
pdf(here(plot_dir, "gene_set_upset.pdf"))
upset(fromList(gene_list[1:3]), 
      order.by = "freq", 
      sets = names(gene_list[1:3]), 
      keep.order = TRUE
)
dev.off()

## compare to 
pdf(here(plot_dir, "gene_set_upset_markers.pdf"))
upset(fromList(gene_list), 
      order.by = "freq", 
      sets = names(gene_list), 
      keep.order = TRUE
)

dev.off()


 
#### Compute the gene set enrichment results ####
modeling_results <- fetch_data(type = "spatialDLPFC_Visium_modeling_results")
gene_list_enrichment <- gene_set_enrichment(
    gene_list = gene_list[1:3],
    modeling_results = modeling_results,
    model_type = "enrichment"
)

dlpfc_anno <- read.csv(here("processed-data", "00_project_prep", "spatialDLPFC_Data", "bayesSpace_layer_annotations.csv")) |>
    dplyr::filter(bayesSpace == "k09") |>
    select(SpD = cluster,
           layer = layer_annotation,
           layer_combo)

gene_list_enrichment$test <- dlpfc_anno$layer_combo[match(gene_list_enrichment$test, dlpfc_anno$SpD)]

pdf(here(plot_dir, "gene_set_enrichment.pdf"))
gene_set_enrichment_plot(
    gene_list_enrichment,
    xlabs = unique(gene_list_enrichment$ID),
    )
dev.off()

#### Check top marker genes for each layer ####

dlpfc_layers_top100_sets <- dlpfc_layers_top100 |>
    mutate(top.svg = ensembl %in% top.svg,
           top.hvg = ensembl %in% top.hvgs$p1,
           top.hdg = ensembl %in% top.hdgs$`1k`) |>
    mutate(layer = ifelse(is.na(SpD), layer, paste(layer,"~", SpD)))

dlpfc_layers_top100_sets |>
    pivot_longer(!c("ensembl", "gene", "dataset", "layer", "top", "marker_anno", "SpD", "layer_combo"),
                 names_to = "gene_set", values_to = "member") |>
    filter(top <= 25) |>
    count(dataset, gene_set, member) |>
    pivot_wider(names_from = "member", values_from = "n") |>
    mutate(precent_TRUE = 100*`TRUE`/(`TRUE` + `FALSE`))

# dataset      gene_set `FALSE` `TRUE` precent_TRUE
# <chr>        <chr>      <int>  <int>        <dbl>
# 1 HumanPilot   top.hdg      139     36         20.6
# 2 HumanPilot   top.hvg      115     60         34.3
# 3 HumanPilot   top.svg       93     82         46.9
# 4 spatialDLPFC top.hdg      111    114         50.7
# 5 spatialDLPFC top.hvg       95    130         57.8
# 6 spatialDLPFC top.svg       66    159         70.7

top_marker_set_tile <-  dlpfc_layers_top100_sets |>
    # filter(top <= 25) |>
    pivot_longer(!c("ensembl", "gene", "dataset", "layer", "top", "marker_anno", "SpD", "layer_combo"),
                 names_to = "gene_set", values_to = "member") |>
    group_by(dataset, layer, gene_set) |>
    count(member) |>
    ggplot(aes(x = member, y = layer, fill = n)) +
    geom_tile() +
    geom_text(aes(label = n)) +
    facet_grid(dataset~gene_set, scales = "free") +
    theme_bw()

ggsave(top_marker_set_tile, filename = here(plot_dir, "marker_set_tile_top100.png"))
    

# slurmjobs::job_single('13_gather_nnSVG', create_shell = TRUE, memory = '10G', command = "Rscript 13_gather_nnSVG.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()s

